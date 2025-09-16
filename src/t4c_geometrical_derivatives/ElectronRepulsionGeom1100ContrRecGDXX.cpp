#include "ElectronRepulsionGeom1100ContrRecGDXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom11_hrr_electron_repulsion_gdxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_gdxx,
                                            const size_t idx_geom_01_fdxx,
                                            const size_t idx_geom_10_fdxx,
                                            const size_t idx_geom_11_fdxx,
                                            const size_t idx_geom_11_ffxx,
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

            const auto fd_geom_01_off = idx_geom_01_fdxx + i * dcomps + j;

            auto g_0_x_xxx_xx = cbuffer.data(fd_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxx_xy = cbuffer.data(fd_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxx_xz = cbuffer.data(fd_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxx_yy = cbuffer.data(fd_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxx_yz = cbuffer.data(fd_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxx_zz = cbuffer.data(fd_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxy_xx = cbuffer.data(fd_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxy_xy = cbuffer.data(fd_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxy_xz = cbuffer.data(fd_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxy_yy = cbuffer.data(fd_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxy_yz = cbuffer.data(fd_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxy_zz = cbuffer.data(fd_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxz_xx = cbuffer.data(fd_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxz_xy = cbuffer.data(fd_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxz_xz = cbuffer.data(fd_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxz_yy = cbuffer.data(fd_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxz_yz = cbuffer.data(fd_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxz_zz = cbuffer.data(fd_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xyy_xx = cbuffer.data(fd_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xyy_xy = cbuffer.data(fd_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xyy_xz = cbuffer.data(fd_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xyy_yy = cbuffer.data(fd_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xyy_yz = cbuffer.data(fd_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xyy_zz = cbuffer.data(fd_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xyz_xx = cbuffer.data(fd_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xyz_xy = cbuffer.data(fd_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xyz_xz = cbuffer.data(fd_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xyz_yy = cbuffer.data(fd_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xyz_yz = cbuffer.data(fd_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xyz_zz = cbuffer.data(fd_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xzz_xx = cbuffer.data(fd_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xzz_xy = cbuffer.data(fd_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xzz_xz = cbuffer.data(fd_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xzz_yy = cbuffer.data(fd_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xzz_yz = cbuffer.data(fd_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xzz_zz = cbuffer.data(fd_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_yyy_xx = cbuffer.data(fd_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_yyy_xy = cbuffer.data(fd_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_yyy_xz = cbuffer.data(fd_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_yyy_yy = cbuffer.data(fd_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_yyy_yz = cbuffer.data(fd_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_yyy_zz = cbuffer.data(fd_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_yyz_xx = cbuffer.data(fd_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_yyz_xy = cbuffer.data(fd_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_yyz_xz = cbuffer.data(fd_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_yyz_yy = cbuffer.data(fd_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_yyz_yz = cbuffer.data(fd_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_yyz_zz = cbuffer.data(fd_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_yzz_xx = cbuffer.data(fd_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_yzz_xy = cbuffer.data(fd_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_yzz_xz = cbuffer.data(fd_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_yzz_yy = cbuffer.data(fd_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_yzz_yz = cbuffer.data(fd_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_yzz_zz = cbuffer.data(fd_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_zzz_xx = cbuffer.data(fd_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_zzz_xy = cbuffer.data(fd_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_zzz_xz = cbuffer.data(fd_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_zzz_yy = cbuffer.data(fd_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_zzz_yz = cbuffer.data(fd_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_zzz_zz = cbuffer.data(fd_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_y_xxx_xx = cbuffer.data(fd_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_y_xxx_xy = cbuffer.data(fd_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_y_xxx_xz = cbuffer.data(fd_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_y_xxx_yy = cbuffer.data(fd_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_y_xxx_yz = cbuffer.data(fd_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_y_xxx_zz = cbuffer.data(fd_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_y_xxy_xx = cbuffer.data(fd_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_y_xxy_xy = cbuffer.data(fd_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_y_xxy_xz = cbuffer.data(fd_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_y_xxy_yy = cbuffer.data(fd_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_y_xxy_yz = cbuffer.data(fd_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_y_xxy_zz = cbuffer.data(fd_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_y_xxz_xx = cbuffer.data(fd_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_y_xxz_xy = cbuffer.data(fd_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_y_xxz_xz = cbuffer.data(fd_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_y_xxz_yy = cbuffer.data(fd_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_y_xxz_yz = cbuffer.data(fd_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_y_xxz_zz = cbuffer.data(fd_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_y_xyy_xx = cbuffer.data(fd_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_y_xyy_xy = cbuffer.data(fd_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_y_xyy_xz = cbuffer.data(fd_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_y_xyy_yy = cbuffer.data(fd_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_y_xyy_yz = cbuffer.data(fd_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_y_xyy_zz = cbuffer.data(fd_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_y_xyz_xx = cbuffer.data(fd_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_y_xyz_xy = cbuffer.data(fd_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_y_xyz_xz = cbuffer.data(fd_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_y_xyz_yy = cbuffer.data(fd_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_y_xyz_yz = cbuffer.data(fd_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_y_xyz_zz = cbuffer.data(fd_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_y_xzz_xx = cbuffer.data(fd_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_y_xzz_xy = cbuffer.data(fd_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_y_xzz_xz = cbuffer.data(fd_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_y_xzz_yy = cbuffer.data(fd_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_y_xzz_yz = cbuffer.data(fd_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_y_xzz_zz = cbuffer.data(fd_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_y_yyy_xx = cbuffer.data(fd_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_y_yyy_xy = cbuffer.data(fd_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_y_yyy_xz = cbuffer.data(fd_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_y_yyy_yy = cbuffer.data(fd_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_y_yyy_yz = cbuffer.data(fd_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_y_yyy_zz = cbuffer.data(fd_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_y_yyz_xx = cbuffer.data(fd_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_y_yyz_xy = cbuffer.data(fd_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_y_yyz_xz = cbuffer.data(fd_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_y_yyz_yy = cbuffer.data(fd_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_y_yyz_yz = cbuffer.data(fd_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_y_yyz_zz = cbuffer.data(fd_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_y_yzz_xx = cbuffer.data(fd_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_yzz_xy = cbuffer.data(fd_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_yzz_xz = cbuffer.data(fd_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_y_yzz_yy = cbuffer.data(fd_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_y_yzz_yz = cbuffer.data(fd_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_yzz_zz = cbuffer.data(fd_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_y_zzz_xx = cbuffer.data(fd_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_zzz_xy = cbuffer.data(fd_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_zzz_xz = cbuffer.data(fd_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_y_zzz_yy = cbuffer.data(fd_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_zzz_yz = cbuffer.data(fd_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_zzz_zz = cbuffer.data(fd_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_z_xxx_xx = cbuffer.data(fd_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_z_xxx_xy = cbuffer.data(fd_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_z_xxx_xz = cbuffer.data(fd_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_z_xxx_yy = cbuffer.data(fd_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_z_xxx_yz = cbuffer.data(fd_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_z_xxx_zz = cbuffer.data(fd_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_z_xxy_xx = cbuffer.data(fd_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_z_xxy_xy = cbuffer.data(fd_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_z_xxy_xz = cbuffer.data(fd_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_z_xxy_yy = cbuffer.data(fd_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_z_xxy_yz = cbuffer.data(fd_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_z_xxy_zz = cbuffer.data(fd_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_z_xxz_xx = cbuffer.data(fd_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_z_xxz_xy = cbuffer.data(fd_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_z_xxz_xz = cbuffer.data(fd_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_z_xxz_yy = cbuffer.data(fd_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_z_xxz_yz = cbuffer.data(fd_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_z_xxz_zz = cbuffer.data(fd_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_z_xyy_xx = cbuffer.data(fd_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_z_xyy_xy = cbuffer.data(fd_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_z_xyy_xz = cbuffer.data(fd_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_z_xyy_yy = cbuffer.data(fd_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_z_xyy_yz = cbuffer.data(fd_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_z_xyy_zz = cbuffer.data(fd_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_z_xyz_xx = cbuffer.data(fd_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_z_xyz_xy = cbuffer.data(fd_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_z_xyz_xz = cbuffer.data(fd_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_z_xyz_yy = cbuffer.data(fd_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_z_xyz_yz = cbuffer.data(fd_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_z_xyz_zz = cbuffer.data(fd_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_z_xzz_xx = cbuffer.data(fd_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_z_xzz_xy = cbuffer.data(fd_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_z_xzz_xz = cbuffer.data(fd_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_z_xzz_yy = cbuffer.data(fd_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_z_xzz_yz = cbuffer.data(fd_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_z_xzz_zz = cbuffer.data(fd_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_z_yyy_xx = cbuffer.data(fd_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_z_yyy_xy = cbuffer.data(fd_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_z_yyy_xz = cbuffer.data(fd_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_z_yyy_yy = cbuffer.data(fd_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_z_yyy_yz = cbuffer.data(fd_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_z_yyy_zz = cbuffer.data(fd_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_z_yyz_xx = cbuffer.data(fd_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_z_yyz_xy = cbuffer.data(fd_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_z_yyz_xz = cbuffer.data(fd_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_z_yyz_yy = cbuffer.data(fd_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_z_yyz_yz = cbuffer.data(fd_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_z_yyz_zz = cbuffer.data(fd_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_z_yzz_xx = cbuffer.data(fd_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_z_yzz_xy = cbuffer.data(fd_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_z_yzz_xz = cbuffer.data(fd_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_z_yzz_yy = cbuffer.data(fd_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_z_yzz_yz = cbuffer.data(fd_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_z_yzz_zz = cbuffer.data(fd_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_z_zzz_xx = cbuffer.data(fd_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_z_zzz_xy = cbuffer.data(fd_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_z_zzz_xz = cbuffer.data(fd_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_z_zzz_yy = cbuffer.data(fd_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_z_zzz_yz = cbuffer.data(fd_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_z_zzz_zz = cbuffer.data(fd_geom_01_off + 179 * ccomps * dcomps);

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

            const auto fd_geom_11_off = idx_geom_11_fdxx + i * dcomps + j;

            auto g_x_x_xxx_xx = cbuffer.data(fd_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xxx_xy = cbuffer.data(fd_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xxx_xz = cbuffer.data(fd_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_xxx_yy = cbuffer.data(fd_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xxx_yz = cbuffer.data(fd_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xxx_zz = cbuffer.data(fd_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_xxy_xx = cbuffer.data(fd_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xxy_xy = cbuffer.data(fd_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xxy_xz = cbuffer.data(fd_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_xxy_yy = cbuffer.data(fd_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_xxy_yz = cbuffer.data(fd_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_xxy_zz = cbuffer.data(fd_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_xxz_xx = cbuffer.data(fd_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_xxz_xy = cbuffer.data(fd_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_xxz_xz = cbuffer.data(fd_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_xxz_yy = cbuffer.data(fd_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_xxz_yz = cbuffer.data(fd_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_xxz_zz = cbuffer.data(fd_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_xyy_xx = cbuffer.data(fd_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_xyy_xy = cbuffer.data(fd_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_xyy_xz = cbuffer.data(fd_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_xyy_yy = cbuffer.data(fd_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_xyy_yz = cbuffer.data(fd_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_xyy_zz = cbuffer.data(fd_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_xyz_xx = cbuffer.data(fd_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_xyz_xy = cbuffer.data(fd_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_xyz_xz = cbuffer.data(fd_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_xyz_yy = cbuffer.data(fd_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_xyz_yz = cbuffer.data(fd_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_xyz_zz = cbuffer.data(fd_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_x_xzz_xx = cbuffer.data(fd_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_xzz_xy = cbuffer.data(fd_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_xzz_xz = cbuffer.data(fd_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_xzz_yy = cbuffer.data(fd_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_xzz_yz = cbuffer.data(fd_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_xzz_zz = cbuffer.data(fd_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_x_yyy_xx = cbuffer.data(fd_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_yyy_xy = cbuffer.data(fd_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_yyy_xz = cbuffer.data(fd_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_yyy_yy = cbuffer.data(fd_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_yyy_yz = cbuffer.data(fd_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_yyy_zz = cbuffer.data(fd_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_x_yyz_xx = cbuffer.data(fd_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_yyz_xy = cbuffer.data(fd_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_yyz_xz = cbuffer.data(fd_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_x_yyz_yy = cbuffer.data(fd_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_x_yyz_yz = cbuffer.data(fd_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_x_yyz_zz = cbuffer.data(fd_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_x_yzz_xx = cbuffer.data(fd_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_x_yzz_xy = cbuffer.data(fd_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_x_yzz_xz = cbuffer.data(fd_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_x_yzz_yy = cbuffer.data(fd_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_x_yzz_yz = cbuffer.data(fd_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_x_yzz_zz = cbuffer.data(fd_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_x_zzz_xx = cbuffer.data(fd_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_x_zzz_xy = cbuffer.data(fd_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_x_zzz_xz = cbuffer.data(fd_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_x_zzz_yy = cbuffer.data(fd_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_x_zzz_yz = cbuffer.data(fd_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_x_zzz_zz = cbuffer.data(fd_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_y_xxx_xx = cbuffer.data(fd_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_y_xxx_xy = cbuffer.data(fd_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_y_xxx_xz = cbuffer.data(fd_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_y_xxx_yy = cbuffer.data(fd_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_y_xxx_yz = cbuffer.data(fd_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_y_xxx_zz = cbuffer.data(fd_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_y_xxy_xx = cbuffer.data(fd_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_y_xxy_xy = cbuffer.data(fd_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_y_xxy_xz = cbuffer.data(fd_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_y_xxy_yy = cbuffer.data(fd_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_y_xxy_yz = cbuffer.data(fd_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_y_xxy_zz = cbuffer.data(fd_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_y_xxz_xx = cbuffer.data(fd_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_y_xxz_xy = cbuffer.data(fd_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_y_xxz_xz = cbuffer.data(fd_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_y_xxz_yy = cbuffer.data(fd_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_y_xxz_yz = cbuffer.data(fd_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_y_xxz_zz = cbuffer.data(fd_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_y_xyy_xx = cbuffer.data(fd_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_y_xyy_xy = cbuffer.data(fd_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_y_xyy_xz = cbuffer.data(fd_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_y_xyy_yy = cbuffer.data(fd_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_y_xyy_yz = cbuffer.data(fd_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_y_xyy_zz = cbuffer.data(fd_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_y_xyz_xx = cbuffer.data(fd_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_y_xyz_xy = cbuffer.data(fd_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_y_xyz_xz = cbuffer.data(fd_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_y_xyz_yy = cbuffer.data(fd_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_y_xyz_yz = cbuffer.data(fd_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_y_xyz_zz = cbuffer.data(fd_geom_11_off + 89 * ccomps * dcomps);

            auto g_x_y_xzz_xx = cbuffer.data(fd_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_y_xzz_xy = cbuffer.data(fd_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_y_xzz_xz = cbuffer.data(fd_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_y_xzz_yy = cbuffer.data(fd_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_y_xzz_yz = cbuffer.data(fd_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_y_xzz_zz = cbuffer.data(fd_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_y_yyy_xx = cbuffer.data(fd_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_y_yyy_xy = cbuffer.data(fd_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_y_yyy_xz = cbuffer.data(fd_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_y_yyy_yy = cbuffer.data(fd_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_y_yyy_yz = cbuffer.data(fd_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_y_yyy_zz = cbuffer.data(fd_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_y_yyz_xx = cbuffer.data(fd_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_y_yyz_xy = cbuffer.data(fd_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_y_yyz_xz = cbuffer.data(fd_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_y_yyz_yy = cbuffer.data(fd_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_y_yyz_yz = cbuffer.data(fd_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_y_yyz_zz = cbuffer.data(fd_geom_11_off + 107 * ccomps * dcomps);

            auto g_x_y_yzz_xx = cbuffer.data(fd_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_y_yzz_xy = cbuffer.data(fd_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_y_yzz_xz = cbuffer.data(fd_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_y_yzz_yy = cbuffer.data(fd_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_y_yzz_yz = cbuffer.data(fd_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_y_yzz_zz = cbuffer.data(fd_geom_11_off + 113 * ccomps * dcomps);

            auto g_x_y_zzz_xx = cbuffer.data(fd_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_y_zzz_xy = cbuffer.data(fd_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_y_zzz_xz = cbuffer.data(fd_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_y_zzz_yy = cbuffer.data(fd_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_y_zzz_yz = cbuffer.data(fd_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_y_zzz_zz = cbuffer.data(fd_geom_11_off + 119 * ccomps * dcomps);

            auto g_x_z_xxx_xx = cbuffer.data(fd_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_z_xxx_xy = cbuffer.data(fd_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_z_xxx_xz = cbuffer.data(fd_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_z_xxx_yy = cbuffer.data(fd_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_z_xxx_yz = cbuffer.data(fd_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_z_xxx_zz = cbuffer.data(fd_geom_11_off + 125 * ccomps * dcomps);

            auto g_x_z_xxy_xx = cbuffer.data(fd_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_z_xxy_xy = cbuffer.data(fd_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_z_xxy_xz = cbuffer.data(fd_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_z_xxy_yy = cbuffer.data(fd_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_z_xxy_yz = cbuffer.data(fd_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_z_xxy_zz = cbuffer.data(fd_geom_11_off + 131 * ccomps * dcomps);

            auto g_x_z_xxz_xx = cbuffer.data(fd_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_z_xxz_xy = cbuffer.data(fd_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_z_xxz_xz = cbuffer.data(fd_geom_11_off + 134 * ccomps * dcomps);

            auto g_x_z_xxz_yy = cbuffer.data(fd_geom_11_off + 135 * ccomps * dcomps);

            auto g_x_z_xxz_yz = cbuffer.data(fd_geom_11_off + 136 * ccomps * dcomps);

            auto g_x_z_xxz_zz = cbuffer.data(fd_geom_11_off + 137 * ccomps * dcomps);

            auto g_x_z_xyy_xx = cbuffer.data(fd_geom_11_off + 138 * ccomps * dcomps);

            auto g_x_z_xyy_xy = cbuffer.data(fd_geom_11_off + 139 * ccomps * dcomps);

            auto g_x_z_xyy_xz = cbuffer.data(fd_geom_11_off + 140 * ccomps * dcomps);

            auto g_x_z_xyy_yy = cbuffer.data(fd_geom_11_off + 141 * ccomps * dcomps);

            auto g_x_z_xyy_yz = cbuffer.data(fd_geom_11_off + 142 * ccomps * dcomps);

            auto g_x_z_xyy_zz = cbuffer.data(fd_geom_11_off + 143 * ccomps * dcomps);

            auto g_x_z_xyz_xx = cbuffer.data(fd_geom_11_off + 144 * ccomps * dcomps);

            auto g_x_z_xyz_xy = cbuffer.data(fd_geom_11_off + 145 * ccomps * dcomps);

            auto g_x_z_xyz_xz = cbuffer.data(fd_geom_11_off + 146 * ccomps * dcomps);

            auto g_x_z_xyz_yy = cbuffer.data(fd_geom_11_off + 147 * ccomps * dcomps);

            auto g_x_z_xyz_yz = cbuffer.data(fd_geom_11_off + 148 * ccomps * dcomps);

            auto g_x_z_xyz_zz = cbuffer.data(fd_geom_11_off + 149 * ccomps * dcomps);

            auto g_x_z_xzz_xx = cbuffer.data(fd_geom_11_off + 150 * ccomps * dcomps);

            auto g_x_z_xzz_xy = cbuffer.data(fd_geom_11_off + 151 * ccomps * dcomps);

            auto g_x_z_xzz_xz = cbuffer.data(fd_geom_11_off + 152 * ccomps * dcomps);

            auto g_x_z_xzz_yy = cbuffer.data(fd_geom_11_off + 153 * ccomps * dcomps);

            auto g_x_z_xzz_yz = cbuffer.data(fd_geom_11_off + 154 * ccomps * dcomps);

            auto g_x_z_xzz_zz = cbuffer.data(fd_geom_11_off + 155 * ccomps * dcomps);

            auto g_x_z_yyy_xx = cbuffer.data(fd_geom_11_off + 156 * ccomps * dcomps);

            auto g_x_z_yyy_xy = cbuffer.data(fd_geom_11_off + 157 * ccomps * dcomps);

            auto g_x_z_yyy_xz = cbuffer.data(fd_geom_11_off + 158 * ccomps * dcomps);

            auto g_x_z_yyy_yy = cbuffer.data(fd_geom_11_off + 159 * ccomps * dcomps);

            auto g_x_z_yyy_yz = cbuffer.data(fd_geom_11_off + 160 * ccomps * dcomps);

            auto g_x_z_yyy_zz = cbuffer.data(fd_geom_11_off + 161 * ccomps * dcomps);

            auto g_x_z_yyz_xx = cbuffer.data(fd_geom_11_off + 162 * ccomps * dcomps);

            auto g_x_z_yyz_xy = cbuffer.data(fd_geom_11_off + 163 * ccomps * dcomps);

            auto g_x_z_yyz_xz = cbuffer.data(fd_geom_11_off + 164 * ccomps * dcomps);

            auto g_x_z_yyz_yy = cbuffer.data(fd_geom_11_off + 165 * ccomps * dcomps);

            auto g_x_z_yyz_yz = cbuffer.data(fd_geom_11_off + 166 * ccomps * dcomps);

            auto g_x_z_yyz_zz = cbuffer.data(fd_geom_11_off + 167 * ccomps * dcomps);

            auto g_x_z_yzz_xx = cbuffer.data(fd_geom_11_off + 168 * ccomps * dcomps);

            auto g_x_z_yzz_xy = cbuffer.data(fd_geom_11_off + 169 * ccomps * dcomps);

            auto g_x_z_yzz_xz = cbuffer.data(fd_geom_11_off + 170 * ccomps * dcomps);

            auto g_x_z_yzz_yy = cbuffer.data(fd_geom_11_off + 171 * ccomps * dcomps);

            auto g_x_z_yzz_yz = cbuffer.data(fd_geom_11_off + 172 * ccomps * dcomps);

            auto g_x_z_yzz_zz = cbuffer.data(fd_geom_11_off + 173 * ccomps * dcomps);

            auto g_x_z_zzz_xx = cbuffer.data(fd_geom_11_off + 174 * ccomps * dcomps);

            auto g_x_z_zzz_xy = cbuffer.data(fd_geom_11_off + 175 * ccomps * dcomps);

            auto g_x_z_zzz_xz = cbuffer.data(fd_geom_11_off + 176 * ccomps * dcomps);

            auto g_x_z_zzz_yy = cbuffer.data(fd_geom_11_off + 177 * ccomps * dcomps);

            auto g_x_z_zzz_yz = cbuffer.data(fd_geom_11_off + 178 * ccomps * dcomps);

            auto g_x_z_zzz_zz = cbuffer.data(fd_geom_11_off + 179 * ccomps * dcomps);

            auto g_y_x_xxx_xx = cbuffer.data(fd_geom_11_off + 180 * ccomps * dcomps);

            auto g_y_x_xxx_xy = cbuffer.data(fd_geom_11_off + 181 * ccomps * dcomps);

            auto g_y_x_xxx_xz = cbuffer.data(fd_geom_11_off + 182 * ccomps * dcomps);

            auto g_y_x_xxx_yy = cbuffer.data(fd_geom_11_off + 183 * ccomps * dcomps);

            auto g_y_x_xxx_yz = cbuffer.data(fd_geom_11_off + 184 * ccomps * dcomps);

            auto g_y_x_xxx_zz = cbuffer.data(fd_geom_11_off + 185 * ccomps * dcomps);

            auto g_y_x_xxy_xx = cbuffer.data(fd_geom_11_off + 186 * ccomps * dcomps);

            auto g_y_x_xxy_xy = cbuffer.data(fd_geom_11_off + 187 * ccomps * dcomps);

            auto g_y_x_xxy_xz = cbuffer.data(fd_geom_11_off + 188 * ccomps * dcomps);

            auto g_y_x_xxy_yy = cbuffer.data(fd_geom_11_off + 189 * ccomps * dcomps);

            auto g_y_x_xxy_yz = cbuffer.data(fd_geom_11_off + 190 * ccomps * dcomps);

            auto g_y_x_xxy_zz = cbuffer.data(fd_geom_11_off + 191 * ccomps * dcomps);

            auto g_y_x_xxz_xx = cbuffer.data(fd_geom_11_off + 192 * ccomps * dcomps);

            auto g_y_x_xxz_xy = cbuffer.data(fd_geom_11_off + 193 * ccomps * dcomps);

            auto g_y_x_xxz_xz = cbuffer.data(fd_geom_11_off + 194 * ccomps * dcomps);

            auto g_y_x_xxz_yy = cbuffer.data(fd_geom_11_off + 195 * ccomps * dcomps);

            auto g_y_x_xxz_yz = cbuffer.data(fd_geom_11_off + 196 * ccomps * dcomps);

            auto g_y_x_xxz_zz = cbuffer.data(fd_geom_11_off + 197 * ccomps * dcomps);

            auto g_y_x_xyy_xx = cbuffer.data(fd_geom_11_off + 198 * ccomps * dcomps);

            auto g_y_x_xyy_xy = cbuffer.data(fd_geom_11_off + 199 * ccomps * dcomps);

            auto g_y_x_xyy_xz = cbuffer.data(fd_geom_11_off + 200 * ccomps * dcomps);

            auto g_y_x_xyy_yy = cbuffer.data(fd_geom_11_off + 201 * ccomps * dcomps);

            auto g_y_x_xyy_yz = cbuffer.data(fd_geom_11_off + 202 * ccomps * dcomps);

            auto g_y_x_xyy_zz = cbuffer.data(fd_geom_11_off + 203 * ccomps * dcomps);

            auto g_y_x_xyz_xx = cbuffer.data(fd_geom_11_off + 204 * ccomps * dcomps);

            auto g_y_x_xyz_xy = cbuffer.data(fd_geom_11_off + 205 * ccomps * dcomps);

            auto g_y_x_xyz_xz = cbuffer.data(fd_geom_11_off + 206 * ccomps * dcomps);

            auto g_y_x_xyz_yy = cbuffer.data(fd_geom_11_off + 207 * ccomps * dcomps);

            auto g_y_x_xyz_yz = cbuffer.data(fd_geom_11_off + 208 * ccomps * dcomps);

            auto g_y_x_xyz_zz = cbuffer.data(fd_geom_11_off + 209 * ccomps * dcomps);

            auto g_y_x_xzz_xx = cbuffer.data(fd_geom_11_off + 210 * ccomps * dcomps);

            auto g_y_x_xzz_xy = cbuffer.data(fd_geom_11_off + 211 * ccomps * dcomps);

            auto g_y_x_xzz_xz = cbuffer.data(fd_geom_11_off + 212 * ccomps * dcomps);

            auto g_y_x_xzz_yy = cbuffer.data(fd_geom_11_off + 213 * ccomps * dcomps);

            auto g_y_x_xzz_yz = cbuffer.data(fd_geom_11_off + 214 * ccomps * dcomps);

            auto g_y_x_xzz_zz = cbuffer.data(fd_geom_11_off + 215 * ccomps * dcomps);

            auto g_y_x_yyy_xx = cbuffer.data(fd_geom_11_off + 216 * ccomps * dcomps);

            auto g_y_x_yyy_xy = cbuffer.data(fd_geom_11_off + 217 * ccomps * dcomps);

            auto g_y_x_yyy_xz = cbuffer.data(fd_geom_11_off + 218 * ccomps * dcomps);

            auto g_y_x_yyy_yy = cbuffer.data(fd_geom_11_off + 219 * ccomps * dcomps);

            auto g_y_x_yyy_yz = cbuffer.data(fd_geom_11_off + 220 * ccomps * dcomps);

            auto g_y_x_yyy_zz = cbuffer.data(fd_geom_11_off + 221 * ccomps * dcomps);

            auto g_y_x_yyz_xx = cbuffer.data(fd_geom_11_off + 222 * ccomps * dcomps);

            auto g_y_x_yyz_xy = cbuffer.data(fd_geom_11_off + 223 * ccomps * dcomps);

            auto g_y_x_yyz_xz = cbuffer.data(fd_geom_11_off + 224 * ccomps * dcomps);

            auto g_y_x_yyz_yy = cbuffer.data(fd_geom_11_off + 225 * ccomps * dcomps);

            auto g_y_x_yyz_yz = cbuffer.data(fd_geom_11_off + 226 * ccomps * dcomps);

            auto g_y_x_yyz_zz = cbuffer.data(fd_geom_11_off + 227 * ccomps * dcomps);

            auto g_y_x_yzz_xx = cbuffer.data(fd_geom_11_off + 228 * ccomps * dcomps);

            auto g_y_x_yzz_xy = cbuffer.data(fd_geom_11_off + 229 * ccomps * dcomps);

            auto g_y_x_yzz_xz = cbuffer.data(fd_geom_11_off + 230 * ccomps * dcomps);

            auto g_y_x_yzz_yy = cbuffer.data(fd_geom_11_off + 231 * ccomps * dcomps);

            auto g_y_x_yzz_yz = cbuffer.data(fd_geom_11_off + 232 * ccomps * dcomps);

            auto g_y_x_yzz_zz = cbuffer.data(fd_geom_11_off + 233 * ccomps * dcomps);

            auto g_y_x_zzz_xx = cbuffer.data(fd_geom_11_off + 234 * ccomps * dcomps);

            auto g_y_x_zzz_xy = cbuffer.data(fd_geom_11_off + 235 * ccomps * dcomps);

            auto g_y_x_zzz_xz = cbuffer.data(fd_geom_11_off + 236 * ccomps * dcomps);

            auto g_y_x_zzz_yy = cbuffer.data(fd_geom_11_off + 237 * ccomps * dcomps);

            auto g_y_x_zzz_yz = cbuffer.data(fd_geom_11_off + 238 * ccomps * dcomps);

            auto g_y_x_zzz_zz = cbuffer.data(fd_geom_11_off + 239 * ccomps * dcomps);

            auto g_y_y_xxx_xx = cbuffer.data(fd_geom_11_off + 240 * ccomps * dcomps);

            auto g_y_y_xxx_xy = cbuffer.data(fd_geom_11_off + 241 * ccomps * dcomps);

            auto g_y_y_xxx_xz = cbuffer.data(fd_geom_11_off + 242 * ccomps * dcomps);

            auto g_y_y_xxx_yy = cbuffer.data(fd_geom_11_off + 243 * ccomps * dcomps);

            auto g_y_y_xxx_yz = cbuffer.data(fd_geom_11_off + 244 * ccomps * dcomps);

            auto g_y_y_xxx_zz = cbuffer.data(fd_geom_11_off + 245 * ccomps * dcomps);

            auto g_y_y_xxy_xx = cbuffer.data(fd_geom_11_off + 246 * ccomps * dcomps);

            auto g_y_y_xxy_xy = cbuffer.data(fd_geom_11_off + 247 * ccomps * dcomps);

            auto g_y_y_xxy_xz = cbuffer.data(fd_geom_11_off + 248 * ccomps * dcomps);

            auto g_y_y_xxy_yy = cbuffer.data(fd_geom_11_off + 249 * ccomps * dcomps);

            auto g_y_y_xxy_yz = cbuffer.data(fd_geom_11_off + 250 * ccomps * dcomps);

            auto g_y_y_xxy_zz = cbuffer.data(fd_geom_11_off + 251 * ccomps * dcomps);

            auto g_y_y_xxz_xx = cbuffer.data(fd_geom_11_off + 252 * ccomps * dcomps);

            auto g_y_y_xxz_xy = cbuffer.data(fd_geom_11_off + 253 * ccomps * dcomps);

            auto g_y_y_xxz_xz = cbuffer.data(fd_geom_11_off + 254 * ccomps * dcomps);

            auto g_y_y_xxz_yy = cbuffer.data(fd_geom_11_off + 255 * ccomps * dcomps);

            auto g_y_y_xxz_yz = cbuffer.data(fd_geom_11_off + 256 * ccomps * dcomps);

            auto g_y_y_xxz_zz = cbuffer.data(fd_geom_11_off + 257 * ccomps * dcomps);

            auto g_y_y_xyy_xx = cbuffer.data(fd_geom_11_off + 258 * ccomps * dcomps);

            auto g_y_y_xyy_xy = cbuffer.data(fd_geom_11_off + 259 * ccomps * dcomps);

            auto g_y_y_xyy_xz = cbuffer.data(fd_geom_11_off + 260 * ccomps * dcomps);

            auto g_y_y_xyy_yy = cbuffer.data(fd_geom_11_off + 261 * ccomps * dcomps);

            auto g_y_y_xyy_yz = cbuffer.data(fd_geom_11_off + 262 * ccomps * dcomps);

            auto g_y_y_xyy_zz = cbuffer.data(fd_geom_11_off + 263 * ccomps * dcomps);

            auto g_y_y_xyz_xx = cbuffer.data(fd_geom_11_off + 264 * ccomps * dcomps);

            auto g_y_y_xyz_xy = cbuffer.data(fd_geom_11_off + 265 * ccomps * dcomps);

            auto g_y_y_xyz_xz = cbuffer.data(fd_geom_11_off + 266 * ccomps * dcomps);

            auto g_y_y_xyz_yy = cbuffer.data(fd_geom_11_off + 267 * ccomps * dcomps);

            auto g_y_y_xyz_yz = cbuffer.data(fd_geom_11_off + 268 * ccomps * dcomps);

            auto g_y_y_xyz_zz = cbuffer.data(fd_geom_11_off + 269 * ccomps * dcomps);

            auto g_y_y_xzz_xx = cbuffer.data(fd_geom_11_off + 270 * ccomps * dcomps);

            auto g_y_y_xzz_xy = cbuffer.data(fd_geom_11_off + 271 * ccomps * dcomps);

            auto g_y_y_xzz_xz = cbuffer.data(fd_geom_11_off + 272 * ccomps * dcomps);

            auto g_y_y_xzz_yy = cbuffer.data(fd_geom_11_off + 273 * ccomps * dcomps);

            auto g_y_y_xzz_yz = cbuffer.data(fd_geom_11_off + 274 * ccomps * dcomps);

            auto g_y_y_xzz_zz = cbuffer.data(fd_geom_11_off + 275 * ccomps * dcomps);

            auto g_y_y_yyy_xx = cbuffer.data(fd_geom_11_off + 276 * ccomps * dcomps);

            auto g_y_y_yyy_xy = cbuffer.data(fd_geom_11_off + 277 * ccomps * dcomps);

            auto g_y_y_yyy_xz = cbuffer.data(fd_geom_11_off + 278 * ccomps * dcomps);

            auto g_y_y_yyy_yy = cbuffer.data(fd_geom_11_off + 279 * ccomps * dcomps);

            auto g_y_y_yyy_yz = cbuffer.data(fd_geom_11_off + 280 * ccomps * dcomps);

            auto g_y_y_yyy_zz = cbuffer.data(fd_geom_11_off + 281 * ccomps * dcomps);

            auto g_y_y_yyz_xx = cbuffer.data(fd_geom_11_off + 282 * ccomps * dcomps);

            auto g_y_y_yyz_xy = cbuffer.data(fd_geom_11_off + 283 * ccomps * dcomps);

            auto g_y_y_yyz_xz = cbuffer.data(fd_geom_11_off + 284 * ccomps * dcomps);

            auto g_y_y_yyz_yy = cbuffer.data(fd_geom_11_off + 285 * ccomps * dcomps);

            auto g_y_y_yyz_yz = cbuffer.data(fd_geom_11_off + 286 * ccomps * dcomps);

            auto g_y_y_yyz_zz = cbuffer.data(fd_geom_11_off + 287 * ccomps * dcomps);

            auto g_y_y_yzz_xx = cbuffer.data(fd_geom_11_off + 288 * ccomps * dcomps);

            auto g_y_y_yzz_xy = cbuffer.data(fd_geom_11_off + 289 * ccomps * dcomps);

            auto g_y_y_yzz_xz = cbuffer.data(fd_geom_11_off + 290 * ccomps * dcomps);

            auto g_y_y_yzz_yy = cbuffer.data(fd_geom_11_off + 291 * ccomps * dcomps);

            auto g_y_y_yzz_yz = cbuffer.data(fd_geom_11_off + 292 * ccomps * dcomps);

            auto g_y_y_yzz_zz = cbuffer.data(fd_geom_11_off + 293 * ccomps * dcomps);

            auto g_y_y_zzz_xx = cbuffer.data(fd_geom_11_off + 294 * ccomps * dcomps);

            auto g_y_y_zzz_xy = cbuffer.data(fd_geom_11_off + 295 * ccomps * dcomps);

            auto g_y_y_zzz_xz = cbuffer.data(fd_geom_11_off + 296 * ccomps * dcomps);

            auto g_y_y_zzz_yy = cbuffer.data(fd_geom_11_off + 297 * ccomps * dcomps);

            auto g_y_y_zzz_yz = cbuffer.data(fd_geom_11_off + 298 * ccomps * dcomps);

            auto g_y_y_zzz_zz = cbuffer.data(fd_geom_11_off + 299 * ccomps * dcomps);

            auto g_y_z_xxx_xx = cbuffer.data(fd_geom_11_off + 300 * ccomps * dcomps);

            auto g_y_z_xxx_xy = cbuffer.data(fd_geom_11_off + 301 * ccomps * dcomps);

            auto g_y_z_xxx_xz = cbuffer.data(fd_geom_11_off + 302 * ccomps * dcomps);

            auto g_y_z_xxx_yy = cbuffer.data(fd_geom_11_off + 303 * ccomps * dcomps);

            auto g_y_z_xxx_yz = cbuffer.data(fd_geom_11_off + 304 * ccomps * dcomps);

            auto g_y_z_xxx_zz = cbuffer.data(fd_geom_11_off + 305 * ccomps * dcomps);

            auto g_y_z_xxy_xx = cbuffer.data(fd_geom_11_off + 306 * ccomps * dcomps);

            auto g_y_z_xxy_xy = cbuffer.data(fd_geom_11_off + 307 * ccomps * dcomps);

            auto g_y_z_xxy_xz = cbuffer.data(fd_geom_11_off + 308 * ccomps * dcomps);

            auto g_y_z_xxy_yy = cbuffer.data(fd_geom_11_off + 309 * ccomps * dcomps);

            auto g_y_z_xxy_yz = cbuffer.data(fd_geom_11_off + 310 * ccomps * dcomps);

            auto g_y_z_xxy_zz = cbuffer.data(fd_geom_11_off + 311 * ccomps * dcomps);

            auto g_y_z_xxz_xx = cbuffer.data(fd_geom_11_off + 312 * ccomps * dcomps);

            auto g_y_z_xxz_xy = cbuffer.data(fd_geom_11_off + 313 * ccomps * dcomps);

            auto g_y_z_xxz_xz = cbuffer.data(fd_geom_11_off + 314 * ccomps * dcomps);

            auto g_y_z_xxz_yy = cbuffer.data(fd_geom_11_off + 315 * ccomps * dcomps);

            auto g_y_z_xxz_yz = cbuffer.data(fd_geom_11_off + 316 * ccomps * dcomps);

            auto g_y_z_xxz_zz = cbuffer.data(fd_geom_11_off + 317 * ccomps * dcomps);

            auto g_y_z_xyy_xx = cbuffer.data(fd_geom_11_off + 318 * ccomps * dcomps);

            auto g_y_z_xyy_xy = cbuffer.data(fd_geom_11_off + 319 * ccomps * dcomps);

            auto g_y_z_xyy_xz = cbuffer.data(fd_geom_11_off + 320 * ccomps * dcomps);

            auto g_y_z_xyy_yy = cbuffer.data(fd_geom_11_off + 321 * ccomps * dcomps);

            auto g_y_z_xyy_yz = cbuffer.data(fd_geom_11_off + 322 * ccomps * dcomps);

            auto g_y_z_xyy_zz = cbuffer.data(fd_geom_11_off + 323 * ccomps * dcomps);

            auto g_y_z_xyz_xx = cbuffer.data(fd_geom_11_off + 324 * ccomps * dcomps);

            auto g_y_z_xyz_xy = cbuffer.data(fd_geom_11_off + 325 * ccomps * dcomps);

            auto g_y_z_xyz_xz = cbuffer.data(fd_geom_11_off + 326 * ccomps * dcomps);

            auto g_y_z_xyz_yy = cbuffer.data(fd_geom_11_off + 327 * ccomps * dcomps);

            auto g_y_z_xyz_yz = cbuffer.data(fd_geom_11_off + 328 * ccomps * dcomps);

            auto g_y_z_xyz_zz = cbuffer.data(fd_geom_11_off + 329 * ccomps * dcomps);

            auto g_y_z_xzz_xx = cbuffer.data(fd_geom_11_off + 330 * ccomps * dcomps);

            auto g_y_z_xzz_xy = cbuffer.data(fd_geom_11_off + 331 * ccomps * dcomps);

            auto g_y_z_xzz_xz = cbuffer.data(fd_geom_11_off + 332 * ccomps * dcomps);

            auto g_y_z_xzz_yy = cbuffer.data(fd_geom_11_off + 333 * ccomps * dcomps);

            auto g_y_z_xzz_yz = cbuffer.data(fd_geom_11_off + 334 * ccomps * dcomps);

            auto g_y_z_xzz_zz = cbuffer.data(fd_geom_11_off + 335 * ccomps * dcomps);

            auto g_y_z_yyy_xx = cbuffer.data(fd_geom_11_off + 336 * ccomps * dcomps);

            auto g_y_z_yyy_xy = cbuffer.data(fd_geom_11_off + 337 * ccomps * dcomps);

            auto g_y_z_yyy_xz = cbuffer.data(fd_geom_11_off + 338 * ccomps * dcomps);

            auto g_y_z_yyy_yy = cbuffer.data(fd_geom_11_off + 339 * ccomps * dcomps);

            auto g_y_z_yyy_yz = cbuffer.data(fd_geom_11_off + 340 * ccomps * dcomps);

            auto g_y_z_yyy_zz = cbuffer.data(fd_geom_11_off + 341 * ccomps * dcomps);

            auto g_y_z_yyz_xx = cbuffer.data(fd_geom_11_off + 342 * ccomps * dcomps);

            auto g_y_z_yyz_xy = cbuffer.data(fd_geom_11_off + 343 * ccomps * dcomps);

            auto g_y_z_yyz_xz = cbuffer.data(fd_geom_11_off + 344 * ccomps * dcomps);

            auto g_y_z_yyz_yy = cbuffer.data(fd_geom_11_off + 345 * ccomps * dcomps);

            auto g_y_z_yyz_yz = cbuffer.data(fd_geom_11_off + 346 * ccomps * dcomps);

            auto g_y_z_yyz_zz = cbuffer.data(fd_geom_11_off + 347 * ccomps * dcomps);

            auto g_y_z_yzz_xx = cbuffer.data(fd_geom_11_off + 348 * ccomps * dcomps);

            auto g_y_z_yzz_xy = cbuffer.data(fd_geom_11_off + 349 * ccomps * dcomps);

            auto g_y_z_yzz_xz = cbuffer.data(fd_geom_11_off + 350 * ccomps * dcomps);

            auto g_y_z_yzz_yy = cbuffer.data(fd_geom_11_off + 351 * ccomps * dcomps);

            auto g_y_z_yzz_yz = cbuffer.data(fd_geom_11_off + 352 * ccomps * dcomps);

            auto g_y_z_yzz_zz = cbuffer.data(fd_geom_11_off + 353 * ccomps * dcomps);

            auto g_y_z_zzz_xx = cbuffer.data(fd_geom_11_off + 354 * ccomps * dcomps);

            auto g_y_z_zzz_xy = cbuffer.data(fd_geom_11_off + 355 * ccomps * dcomps);

            auto g_y_z_zzz_xz = cbuffer.data(fd_geom_11_off + 356 * ccomps * dcomps);

            auto g_y_z_zzz_yy = cbuffer.data(fd_geom_11_off + 357 * ccomps * dcomps);

            auto g_y_z_zzz_yz = cbuffer.data(fd_geom_11_off + 358 * ccomps * dcomps);

            auto g_y_z_zzz_zz = cbuffer.data(fd_geom_11_off + 359 * ccomps * dcomps);

            auto g_z_x_xxx_xx = cbuffer.data(fd_geom_11_off + 360 * ccomps * dcomps);

            auto g_z_x_xxx_xy = cbuffer.data(fd_geom_11_off + 361 * ccomps * dcomps);

            auto g_z_x_xxx_xz = cbuffer.data(fd_geom_11_off + 362 * ccomps * dcomps);

            auto g_z_x_xxx_yy = cbuffer.data(fd_geom_11_off + 363 * ccomps * dcomps);

            auto g_z_x_xxx_yz = cbuffer.data(fd_geom_11_off + 364 * ccomps * dcomps);

            auto g_z_x_xxx_zz = cbuffer.data(fd_geom_11_off + 365 * ccomps * dcomps);

            auto g_z_x_xxy_xx = cbuffer.data(fd_geom_11_off + 366 * ccomps * dcomps);

            auto g_z_x_xxy_xy = cbuffer.data(fd_geom_11_off + 367 * ccomps * dcomps);

            auto g_z_x_xxy_xz = cbuffer.data(fd_geom_11_off + 368 * ccomps * dcomps);

            auto g_z_x_xxy_yy = cbuffer.data(fd_geom_11_off + 369 * ccomps * dcomps);

            auto g_z_x_xxy_yz = cbuffer.data(fd_geom_11_off + 370 * ccomps * dcomps);

            auto g_z_x_xxy_zz = cbuffer.data(fd_geom_11_off + 371 * ccomps * dcomps);

            auto g_z_x_xxz_xx = cbuffer.data(fd_geom_11_off + 372 * ccomps * dcomps);

            auto g_z_x_xxz_xy = cbuffer.data(fd_geom_11_off + 373 * ccomps * dcomps);

            auto g_z_x_xxz_xz = cbuffer.data(fd_geom_11_off + 374 * ccomps * dcomps);

            auto g_z_x_xxz_yy = cbuffer.data(fd_geom_11_off + 375 * ccomps * dcomps);

            auto g_z_x_xxz_yz = cbuffer.data(fd_geom_11_off + 376 * ccomps * dcomps);

            auto g_z_x_xxz_zz = cbuffer.data(fd_geom_11_off + 377 * ccomps * dcomps);

            auto g_z_x_xyy_xx = cbuffer.data(fd_geom_11_off + 378 * ccomps * dcomps);

            auto g_z_x_xyy_xy = cbuffer.data(fd_geom_11_off + 379 * ccomps * dcomps);

            auto g_z_x_xyy_xz = cbuffer.data(fd_geom_11_off + 380 * ccomps * dcomps);

            auto g_z_x_xyy_yy = cbuffer.data(fd_geom_11_off + 381 * ccomps * dcomps);

            auto g_z_x_xyy_yz = cbuffer.data(fd_geom_11_off + 382 * ccomps * dcomps);

            auto g_z_x_xyy_zz = cbuffer.data(fd_geom_11_off + 383 * ccomps * dcomps);

            auto g_z_x_xyz_xx = cbuffer.data(fd_geom_11_off + 384 * ccomps * dcomps);

            auto g_z_x_xyz_xy = cbuffer.data(fd_geom_11_off + 385 * ccomps * dcomps);

            auto g_z_x_xyz_xz = cbuffer.data(fd_geom_11_off + 386 * ccomps * dcomps);

            auto g_z_x_xyz_yy = cbuffer.data(fd_geom_11_off + 387 * ccomps * dcomps);

            auto g_z_x_xyz_yz = cbuffer.data(fd_geom_11_off + 388 * ccomps * dcomps);

            auto g_z_x_xyz_zz = cbuffer.data(fd_geom_11_off + 389 * ccomps * dcomps);

            auto g_z_x_xzz_xx = cbuffer.data(fd_geom_11_off + 390 * ccomps * dcomps);

            auto g_z_x_xzz_xy = cbuffer.data(fd_geom_11_off + 391 * ccomps * dcomps);

            auto g_z_x_xzz_xz = cbuffer.data(fd_geom_11_off + 392 * ccomps * dcomps);

            auto g_z_x_xzz_yy = cbuffer.data(fd_geom_11_off + 393 * ccomps * dcomps);

            auto g_z_x_xzz_yz = cbuffer.data(fd_geom_11_off + 394 * ccomps * dcomps);

            auto g_z_x_xzz_zz = cbuffer.data(fd_geom_11_off + 395 * ccomps * dcomps);

            auto g_z_x_yyy_xx = cbuffer.data(fd_geom_11_off + 396 * ccomps * dcomps);

            auto g_z_x_yyy_xy = cbuffer.data(fd_geom_11_off + 397 * ccomps * dcomps);

            auto g_z_x_yyy_xz = cbuffer.data(fd_geom_11_off + 398 * ccomps * dcomps);

            auto g_z_x_yyy_yy = cbuffer.data(fd_geom_11_off + 399 * ccomps * dcomps);

            auto g_z_x_yyy_yz = cbuffer.data(fd_geom_11_off + 400 * ccomps * dcomps);

            auto g_z_x_yyy_zz = cbuffer.data(fd_geom_11_off + 401 * ccomps * dcomps);

            auto g_z_x_yyz_xx = cbuffer.data(fd_geom_11_off + 402 * ccomps * dcomps);

            auto g_z_x_yyz_xy = cbuffer.data(fd_geom_11_off + 403 * ccomps * dcomps);

            auto g_z_x_yyz_xz = cbuffer.data(fd_geom_11_off + 404 * ccomps * dcomps);

            auto g_z_x_yyz_yy = cbuffer.data(fd_geom_11_off + 405 * ccomps * dcomps);

            auto g_z_x_yyz_yz = cbuffer.data(fd_geom_11_off + 406 * ccomps * dcomps);

            auto g_z_x_yyz_zz = cbuffer.data(fd_geom_11_off + 407 * ccomps * dcomps);

            auto g_z_x_yzz_xx = cbuffer.data(fd_geom_11_off + 408 * ccomps * dcomps);

            auto g_z_x_yzz_xy = cbuffer.data(fd_geom_11_off + 409 * ccomps * dcomps);

            auto g_z_x_yzz_xz = cbuffer.data(fd_geom_11_off + 410 * ccomps * dcomps);

            auto g_z_x_yzz_yy = cbuffer.data(fd_geom_11_off + 411 * ccomps * dcomps);

            auto g_z_x_yzz_yz = cbuffer.data(fd_geom_11_off + 412 * ccomps * dcomps);

            auto g_z_x_yzz_zz = cbuffer.data(fd_geom_11_off + 413 * ccomps * dcomps);

            auto g_z_x_zzz_xx = cbuffer.data(fd_geom_11_off + 414 * ccomps * dcomps);

            auto g_z_x_zzz_xy = cbuffer.data(fd_geom_11_off + 415 * ccomps * dcomps);

            auto g_z_x_zzz_xz = cbuffer.data(fd_geom_11_off + 416 * ccomps * dcomps);

            auto g_z_x_zzz_yy = cbuffer.data(fd_geom_11_off + 417 * ccomps * dcomps);

            auto g_z_x_zzz_yz = cbuffer.data(fd_geom_11_off + 418 * ccomps * dcomps);

            auto g_z_x_zzz_zz = cbuffer.data(fd_geom_11_off + 419 * ccomps * dcomps);

            auto g_z_y_xxx_xx = cbuffer.data(fd_geom_11_off + 420 * ccomps * dcomps);

            auto g_z_y_xxx_xy = cbuffer.data(fd_geom_11_off + 421 * ccomps * dcomps);

            auto g_z_y_xxx_xz = cbuffer.data(fd_geom_11_off + 422 * ccomps * dcomps);

            auto g_z_y_xxx_yy = cbuffer.data(fd_geom_11_off + 423 * ccomps * dcomps);

            auto g_z_y_xxx_yz = cbuffer.data(fd_geom_11_off + 424 * ccomps * dcomps);

            auto g_z_y_xxx_zz = cbuffer.data(fd_geom_11_off + 425 * ccomps * dcomps);

            auto g_z_y_xxy_xx = cbuffer.data(fd_geom_11_off + 426 * ccomps * dcomps);

            auto g_z_y_xxy_xy = cbuffer.data(fd_geom_11_off + 427 * ccomps * dcomps);

            auto g_z_y_xxy_xz = cbuffer.data(fd_geom_11_off + 428 * ccomps * dcomps);

            auto g_z_y_xxy_yy = cbuffer.data(fd_geom_11_off + 429 * ccomps * dcomps);

            auto g_z_y_xxy_yz = cbuffer.data(fd_geom_11_off + 430 * ccomps * dcomps);

            auto g_z_y_xxy_zz = cbuffer.data(fd_geom_11_off + 431 * ccomps * dcomps);

            auto g_z_y_xxz_xx = cbuffer.data(fd_geom_11_off + 432 * ccomps * dcomps);

            auto g_z_y_xxz_xy = cbuffer.data(fd_geom_11_off + 433 * ccomps * dcomps);

            auto g_z_y_xxz_xz = cbuffer.data(fd_geom_11_off + 434 * ccomps * dcomps);

            auto g_z_y_xxz_yy = cbuffer.data(fd_geom_11_off + 435 * ccomps * dcomps);

            auto g_z_y_xxz_yz = cbuffer.data(fd_geom_11_off + 436 * ccomps * dcomps);

            auto g_z_y_xxz_zz = cbuffer.data(fd_geom_11_off + 437 * ccomps * dcomps);

            auto g_z_y_xyy_xx = cbuffer.data(fd_geom_11_off + 438 * ccomps * dcomps);

            auto g_z_y_xyy_xy = cbuffer.data(fd_geom_11_off + 439 * ccomps * dcomps);

            auto g_z_y_xyy_xz = cbuffer.data(fd_geom_11_off + 440 * ccomps * dcomps);

            auto g_z_y_xyy_yy = cbuffer.data(fd_geom_11_off + 441 * ccomps * dcomps);

            auto g_z_y_xyy_yz = cbuffer.data(fd_geom_11_off + 442 * ccomps * dcomps);

            auto g_z_y_xyy_zz = cbuffer.data(fd_geom_11_off + 443 * ccomps * dcomps);

            auto g_z_y_xyz_xx = cbuffer.data(fd_geom_11_off + 444 * ccomps * dcomps);

            auto g_z_y_xyz_xy = cbuffer.data(fd_geom_11_off + 445 * ccomps * dcomps);

            auto g_z_y_xyz_xz = cbuffer.data(fd_geom_11_off + 446 * ccomps * dcomps);

            auto g_z_y_xyz_yy = cbuffer.data(fd_geom_11_off + 447 * ccomps * dcomps);

            auto g_z_y_xyz_yz = cbuffer.data(fd_geom_11_off + 448 * ccomps * dcomps);

            auto g_z_y_xyz_zz = cbuffer.data(fd_geom_11_off + 449 * ccomps * dcomps);

            auto g_z_y_xzz_xx = cbuffer.data(fd_geom_11_off + 450 * ccomps * dcomps);

            auto g_z_y_xzz_xy = cbuffer.data(fd_geom_11_off + 451 * ccomps * dcomps);

            auto g_z_y_xzz_xz = cbuffer.data(fd_geom_11_off + 452 * ccomps * dcomps);

            auto g_z_y_xzz_yy = cbuffer.data(fd_geom_11_off + 453 * ccomps * dcomps);

            auto g_z_y_xzz_yz = cbuffer.data(fd_geom_11_off + 454 * ccomps * dcomps);

            auto g_z_y_xzz_zz = cbuffer.data(fd_geom_11_off + 455 * ccomps * dcomps);

            auto g_z_y_yyy_xx = cbuffer.data(fd_geom_11_off + 456 * ccomps * dcomps);

            auto g_z_y_yyy_xy = cbuffer.data(fd_geom_11_off + 457 * ccomps * dcomps);

            auto g_z_y_yyy_xz = cbuffer.data(fd_geom_11_off + 458 * ccomps * dcomps);

            auto g_z_y_yyy_yy = cbuffer.data(fd_geom_11_off + 459 * ccomps * dcomps);

            auto g_z_y_yyy_yz = cbuffer.data(fd_geom_11_off + 460 * ccomps * dcomps);

            auto g_z_y_yyy_zz = cbuffer.data(fd_geom_11_off + 461 * ccomps * dcomps);

            auto g_z_y_yyz_xx = cbuffer.data(fd_geom_11_off + 462 * ccomps * dcomps);

            auto g_z_y_yyz_xy = cbuffer.data(fd_geom_11_off + 463 * ccomps * dcomps);

            auto g_z_y_yyz_xz = cbuffer.data(fd_geom_11_off + 464 * ccomps * dcomps);

            auto g_z_y_yyz_yy = cbuffer.data(fd_geom_11_off + 465 * ccomps * dcomps);

            auto g_z_y_yyz_yz = cbuffer.data(fd_geom_11_off + 466 * ccomps * dcomps);

            auto g_z_y_yyz_zz = cbuffer.data(fd_geom_11_off + 467 * ccomps * dcomps);

            auto g_z_y_yzz_xx = cbuffer.data(fd_geom_11_off + 468 * ccomps * dcomps);

            auto g_z_y_yzz_xy = cbuffer.data(fd_geom_11_off + 469 * ccomps * dcomps);

            auto g_z_y_yzz_xz = cbuffer.data(fd_geom_11_off + 470 * ccomps * dcomps);

            auto g_z_y_yzz_yy = cbuffer.data(fd_geom_11_off + 471 * ccomps * dcomps);

            auto g_z_y_yzz_yz = cbuffer.data(fd_geom_11_off + 472 * ccomps * dcomps);

            auto g_z_y_yzz_zz = cbuffer.data(fd_geom_11_off + 473 * ccomps * dcomps);

            auto g_z_y_zzz_xx = cbuffer.data(fd_geom_11_off + 474 * ccomps * dcomps);

            auto g_z_y_zzz_xy = cbuffer.data(fd_geom_11_off + 475 * ccomps * dcomps);

            auto g_z_y_zzz_xz = cbuffer.data(fd_geom_11_off + 476 * ccomps * dcomps);

            auto g_z_y_zzz_yy = cbuffer.data(fd_geom_11_off + 477 * ccomps * dcomps);

            auto g_z_y_zzz_yz = cbuffer.data(fd_geom_11_off + 478 * ccomps * dcomps);

            auto g_z_y_zzz_zz = cbuffer.data(fd_geom_11_off + 479 * ccomps * dcomps);

            auto g_z_z_xxx_xx = cbuffer.data(fd_geom_11_off + 480 * ccomps * dcomps);

            auto g_z_z_xxx_xy = cbuffer.data(fd_geom_11_off + 481 * ccomps * dcomps);

            auto g_z_z_xxx_xz = cbuffer.data(fd_geom_11_off + 482 * ccomps * dcomps);

            auto g_z_z_xxx_yy = cbuffer.data(fd_geom_11_off + 483 * ccomps * dcomps);

            auto g_z_z_xxx_yz = cbuffer.data(fd_geom_11_off + 484 * ccomps * dcomps);

            auto g_z_z_xxx_zz = cbuffer.data(fd_geom_11_off + 485 * ccomps * dcomps);

            auto g_z_z_xxy_xx = cbuffer.data(fd_geom_11_off + 486 * ccomps * dcomps);

            auto g_z_z_xxy_xy = cbuffer.data(fd_geom_11_off + 487 * ccomps * dcomps);

            auto g_z_z_xxy_xz = cbuffer.data(fd_geom_11_off + 488 * ccomps * dcomps);

            auto g_z_z_xxy_yy = cbuffer.data(fd_geom_11_off + 489 * ccomps * dcomps);

            auto g_z_z_xxy_yz = cbuffer.data(fd_geom_11_off + 490 * ccomps * dcomps);

            auto g_z_z_xxy_zz = cbuffer.data(fd_geom_11_off + 491 * ccomps * dcomps);

            auto g_z_z_xxz_xx = cbuffer.data(fd_geom_11_off + 492 * ccomps * dcomps);

            auto g_z_z_xxz_xy = cbuffer.data(fd_geom_11_off + 493 * ccomps * dcomps);

            auto g_z_z_xxz_xz = cbuffer.data(fd_geom_11_off + 494 * ccomps * dcomps);

            auto g_z_z_xxz_yy = cbuffer.data(fd_geom_11_off + 495 * ccomps * dcomps);

            auto g_z_z_xxz_yz = cbuffer.data(fd_geom_11_off + 496 * ccomps * dcomps);

            auto g_z_z_xxz_zz = cbuffer.data(fd_geom_11_off + 497 * ccomps * dcomps);

            auto g_z_z_xyy_xx = cbuffer.data(fd_geom_11_off + 498 * ccomps * dcomps);

            auto g_z_z_xyy_xy = cbuffer.data(fd_geom_11_off + 499 * ccomps * dcomps);

            auto g_z_z_xyy_xz = cbuffer.data(fd_geom_11_off + 500 * ccomps * dcomps);

            auto g_z_z_xyy_yy = cbuffer.data(fd_geom_11_off + 501 * ccomps * dcomps);

            auto g_z_z_xyy_yz = cbuffer.data(fd_geom_11_off + 502 * ccomps * dcomps);

            auto g_z_z_xyy_zz = cbuffer.data(fd_geom_11_off + 503 * ccomps * dcomps);

            auto g_z_z_xyz_xx = cbuffer.data(fd_geom_11_off + 504 * ccomps * dcomps);

            auto g_z_z_xyz_xy = cbuffer.data(fd_geom_11_off + 505 * ccomps * dcomps);

            auto g_z_z_xyz_xz = cbuffer.data(fd_geom_11_off + 506 * ccomps * dcomps);

            auto g_z_z_xyz_yy = cbuffer.data(fd_geom_11_off + 507 * ccomps * dcomps);

            auto g_z_z_xyz_yz = cbuffer.data(fd_geom_11_off + 508 * ccomps * dcomps);

            auto g_z_z_xyz_zz = cbuffer.data(fd_geom_11_off + 509 * ccomps * dcomps);

            auto g_z_z_xzz_xx = cbuffer.data(fd_geom_11_off + 510 * ccomps * dcomps);

            auto g_z_z_xzz_xy = cbuffer.data(fd_geom_11_off + 511 * ccomps * dcomps);

            auto g_z_z_xzz_xz = cbuffer.data(fd_geom_11_off + 512 * ccomps * dcomps);

            auto g_z_z_xzz_yy = cbuffer.data(fd_geom_11_off + 513 * ccomps * dcomps);

            auto g_z_z_xzz_yz = cbuffer.data(fd_geom_11_off + 514 * ccomps * dcomps);

            auto g_z_z_xzz_zz = cbuffer.data(fd_geom_11_off + 515 * ccomps * dcomps);

            auto g_z_z_yyy_xx = cbuffer.data(fd_geom_11_off + 516 * ccomps * dcomps);

            auto g_z_z_yyy_xy = cbuffer.data(fd_geom_11_off + 517 * ccomps * dcomps);

            auto g_z_z_yyy_xz = cbuffer.data(fd_geom_11_off + 518 * ccomps * dcomps);

            auto g_z_z_yyy_yy = cbuffer.data(fd_geom_11_off + 519 * ccomps * dcomps);

            auto g_z_z_yyy_yz = cbuffer.data(fd_geom_11_off + 520 * ccomps * dcomps);

            auto g_z_z_yyy_zz = cbuffer.data(fd_geom_11_off + 521 * ccomps * dcomps);

            auto g_z_z_yyz_xx = cbuffer.data(fd_geom_11_off + 522 * ccomps * dcomps);

            auto g_z_z_yyz_xy = cbuffer.data(fd_geom_11_off + 523 * ccomps * dcomps);

            auto g_z_z_yyz_xz = cbuffer.data(fd_geom_11_off + 524 * ccomps * dcomps);

            auto g_z_z_yyz_yy = cbuffer.data(fd_geom_11_off + 525 * ccomps * dcomps);

            auto g_z_z_yyz_yz = cbuffer.data(fd_geom_11_off + 526 * ccomps * dcomps);

            auto g_z_z_yyz_zz = cbuffer.data(fd_geom_11_off + 527 * ccomps * dcomps);

            auto g_z_z_yzz_xx = cbuffer.data(fd_geom_11_off + 528 * ccomps * dcomps);

            auto g_z_z_yzz_xy = cbuffer.data(fd_geom_11_off + 529 * ccomps * dcomps);

            auto g_z_z_yzz_xz = cbuffer.data(fd_geom_11_off + 530 * ccomps * dcomps);

            auto g_z_z_yzz_yy = cbuffer.data(fd_geom_11_off + 531 * ccomps * dcomps);

            auto g_z_z_yzz_yz = cbuffer.data(fd_geom_11_off + 532 * ccomps * dcomps);

            auto g_z_z_yzz_zz = cbuffer.data(fd_geom_11_off + 533 * ccomps * dcomps);

            auto g_z_z_zzz_xx = cbuffer.data(fd_geom_11_off + 534 * ccomps * dcomps);

            auto g_z_z_zzz_xy = cbuffer.data(fd_geom_11_off + 535 * ccomps * dcomps);

            auto g_z_z_zzz_xz = cbuffer.data(fd_geom_11_off + 536 * ccomps * dcomps);

            auto g_z_z_zzz_yy = cbuffer.data(fd_geom_11_off + 537 * ccomps * dcomps);

            auto g_z_z_zzz_yz = cbuffer.data(fd_geom_11_off + 538 * ccomps * dcomps);

            auto g_z_z_zzz_zz = cbuffer.data(fd_geom_11_off + 539 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FFSS

            const auto ff_geom_11_off = idx_geom_11_ffxx + i * dcomps + j;

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

            /// set up bra offset for contr_buffer_gdxx

            const auto gd_geom_11_off = idx_geom_11_gdxx + i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxxx_xx = cbuffer.data(gd_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xxxx_xy = cbuffer.data(gd_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xxxx_xz = cbuffer.data(gd_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_xxxx_yy = cbuffer.data(gd_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xxxx_yz = cbuffer.data(gd_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xxxx_zz = cbuffer.data(gd_geom_11_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxx_xx, g_0_x_xxx_xy, g_0_x_xxx_xz, g_0_x_xxx_yy, g_0_x_xxx_yz, g_0_x_xxx_zz, g_x_0_xxx_xx, g_x_0_xxx_xy, g_x_0_xxx_xz, g_x_0_xxx_yy, g_x_0_xxx_yz, g_x_0_xxx_zz, g_x_x_xxx_xx, g_x_x_xxx_xxx, g_x_x_xxx_xxy, g_x_x_xxx_xxz, g_x_x_xxx_xy, g_x_x_xxx_xyy, g_x_x_xxx_xyz, g_x_x_xxx_xz, g_x_x_xxx_xzz, g_x_x_xxx_yy, g_x_x_xxx_yz, g_x_x_xxx_zz, g_x_x_xxxx_xx, g_x_x_xxxx_xy, g_x_x_xxxx_xz, g_x_x_xxxx_yy, g_x_x_xxxx_yz, g_x_x_xxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxxx_xx[k] = -g_0_x_xxx_xx[k] + g_x_0_xxx_xx[k] - g_x_x_xxx_xx[k] * ab_x + g_x_x_xxx_xxx[k];

                g_x_x_xxxx_xy[k] = -g_0_x_xxx_xy[k] + g_x_0_xxx_xy[k] - g_x_x_xxx_xy[k] * ab_x + g_x_x_xxx_xxy[k];

                g_x_x_xxxx_xz[k] = -g_0_x_xxx_xz[k] + g_x_0_xxx_xz[k] - g_x_x_xxx_xz[k] * ab_x + g_x_x_xxx_xxz[k];

                g_x_x_xxxx_yy[k] = -g_0_x_xxx_yy[k] + g_x_0_xxx_yy[k] - g_x_x_xxx_yy[k] * ab_x + g_x_x_xxx_xyy[k];

                g_x_x_xxxx_yz[k] = -g_0_x_xxx_yz[k] + g_x_0_xxx_yz[k] - g_x_x_xxx_yz[k] * ab_x + g_x_x_xxx_xyz[k];

                g_x_x_xxxx_zz[k] = -g_0_x_xxx_zz[k] + g_x_0_xxx_zz[k] - g_x_x_xxx_zz[k] * ab_x + g_x_x_xxx_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxxy_xx = cbuffer.data(gd_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xxxy_xy = cbuffer.data(gd_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xxxy_xz = cbuffer.data(gd_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_xxxy_yy = cbuffer.data(gd_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_xxxy_yz = cbuffer.data(gd_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_xxxy_zz = cbuffer.data(gd_geom_11_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxx_xx, g_x_x_xxx_xxy, g_x_x_xxx_xy, g_x_x_xxx_xyy, g_x_x_xxx_xyz, g_x_x_xxx_xz, g_x_x_xxx_yy, g_x_x_xxx_yyy, g_x_x_xxx_yyz, g_x_x_xxx_yz, g_x_x_xxx_yzz, g_x_x_xxx_zz, g_x_x_xxxy_xx, g_x_x_xxxy_xy, g_x_x_xxxy_xz, g_x_x_xxxy_yy, g_x_x_xxxy_yz, g_x_x_xxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxxy_xx[k] = -g_x_x_xxx_xx[k] * ab_y + g_x_x_xxx_xxy[k];

                g_x_x_xxxy_xy[k] = -g_x_x_xxx_xy[k] * ab_y + g_x_x_xxx_xyy[k];

                g_x_x_xxxy_xz[k] = -g_x_x_xxx_xz[k] * ab_y + g_x_x_xxx_xyz[k];

                g_x_x_xxxy_yy[k] = -g_x_x_xxx_yy[k] * ab_y + g_x_x_xxx_yyy[k];

                g_x_x_xxxy_yz[k] = -g_x_x_xxx_yz[k] * ab_y + g_x_x_xxx_yyz[k];

                g_x_x_xxxy_zz[k] = -g_x_x_xxx_zz[k] * ab_y + g_x_x_xxx_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxxz_xx = cbuffer.data(gd_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_xxxz_xy = cbuffer.data(gd_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_xxxz_xz = cbuffer.data(gd_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_xxxz_yy = cbuffer.data(gd_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_xxxz_yz = cbuffer.data(gd_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_xxxz_zz = cbuffer.data(gd_geom_11_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxx_xx, g_x_x_xxx_xxz, g_x_x_xxx_xy, g_x_x_xxx_xyz, g_x_x_xxx_xz, g_x_x_xxx_xzz, g_x_x_xxx_yy, g_x_x_xxx_yyz, g_x_x_xxx_yz, g_x_x_xxx_yzz, g_x_x_xxx_zz, g_x_x_xxx_zzz, g_x_x_xxxz_xx, g_x_x_xxxz_xy, g_x_x_xxxz_xz, g_x_x_xxxz_yy, g_x_x_xxxz_yz, g_x_x_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxxz_xx[k] = -g_x_x_xxx_xx[k] * ab_z + g_x_x_xxx_xxz[k];

                g_x_x_xxxz_xy[k] = -g_x_x_xxx_xy[k] * ab_z + g_x_x_xxx_xyz[k];

                g_x_x_xxxz_xz[k] = -g_x_x_xxx_xz[k] * ab_z + g_x_x_xxx_xzz[k];

                g_x_x_xxxz_yy[k] = -g_x_x_xxx_yy[k] * ab_z + g_x_x_xxx_yyz[k];

                g_x_x_xxxz_yz[k] = -g_x_x_xxx_yz[k] * ab_z + g_x_x_xxx_yzz[k];

                g_x_x_xxxz_zz[k] = -g_x_x_xxx_zz[k] * ab_z + g_x_x_xxx_zzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxyy_xx = cbuffer.data(gd_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_xxyy_xy = cbuffer.data(gd_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_xxyy_xz = cbuffer.data(gd_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_xxyy_yy = cbuffer.data(gd_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_xxyy_yz = cbuffer.data(gd_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_xxyy_zz = cbuffer.data(gd_geom_11_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxy_xx, g_x_x_xxy_xxy, g_x_x_xxy_xy, g_x_x_xxy_xyy, g_x_x_xxy_xyz, g_x_x_xxy_xz, g_x_x_xxy_yy, g_x_x_xxy_yyy, g_x_x_xxy_yyz, g_x_x_xxy_yz, g_x_x_xxy_yzz, g_x_x_xxy_zz, g_x_x_xxyy_xx, g_x_x_xxyy_xy, g_x_x_xxyy_xz, g_x_x_xxyy_yy, g_x_x_xxyy_yz, g_x_x_xxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxyy_xx[k] = -g_x_x_xxy_xx[k] * ab_y + g_x_x_xxy_xxy[k];

                g_x_x_xxyy_xy[k] = -g_x_x_xxy_xy[k] * ab_y + g_x_x_xxy_xyy[k];

                g_x_x_xxyy_xz[k] = -g_x_x_xxy_xz[k] * ab_y + g_x_x_xxy_xyz[k];

                g_x_x_xxyy_yy[k] = -g_x_x_xxy_yy[k] * ab_y + g_x_x_xxy_yyy[k];

                g_x_x_xxyy_yz[k] = -g_x_x_xxy_yz[k] * ab_y + g_x_x_xxy_yyz[k];

                g_x_x_xxyy_zz[k] = -g_x_x_xxy_zz[k] * ab_y + g_x_x_xxy_yzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxyz_xx = cbuffer.data(gd_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_xxyz_xy = cbuffer.data(gd_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_xxyz_xz = cbuffer.data(gd_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_xxyz_yy = cbuffer.data(gd_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_xxyz_yz = cbuffer.data(gd_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_xxyz_zz = cbuffer.data(gd_geom_11_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxyz_xx, g_x_x_xxyz_xy, g_x_x_xxyz_xz, g_x_x_xxyz_yy, g_x_x_xxyz_yz, g_x_x_xxyz_zz, g_x_x_xxz_xx, g_x_x_xxz_xxy, g_x_x_xxz_xy, g_x_x_xxz_xyy, g_x_x_xxz_xyz, g_x_x_xxz_xz, g_x_x_xxz_yy, g_x_x_xxz_yyy, g_x_x_xxz_yyz, g_x_x_xxz_yz, g_x_x_xxz_yzz, g_x_x_xxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxyz_xx[k] = -g_x_x_xxz_xx[k] * ab_y + g_x_x_xxz_xxy[k];

                g_x_x_xxyz_xy[k] = -g_x_x_xxz_xy[k] * ab_y + g_x_x_xxz_xyy[k];

                g_x_x_xxyz_xz[k] = -g_x_x_xxz_xz[k] * ab_y + g_x_x_xxz_xyz[k];

                g_x_x_xxyz_yy[k] = -g_x_x_xxz_yy[k] * ab_y + g_x_x_xxz_yyy[k];

                g_x_x_xxyz_yz[k] = -g_x_x_xxz_yz[k] * ab_y + g_x_x_xxz_yyz[k];

                g_x_x_xxyz_zz[k] = -g_x_x_xxz_zz[k] * ab_y + g_x_x_xxz_yzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxzz_xx = cbuffer.data(gd_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_xxzz_xy = cbuffer.data(gd_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_xxzz_xz = cbuffer.data(gd_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_xxzz_yy = cbuffer.data(gd_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_xxzz_yz = cbuffer.data(gd_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_xxzz_zz = cbuffer.data(gd_geom_11_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxz_xx, g_x_x_xxz_xxz, g_x_x_xxz_xy, g_x_x_xxz_xyz, g_x_x_xxz_xz, g_x_x_xxz_xzz, g_x_x_xxz_yy, g_x_x_xxz_yyz, g_x_x_xxz_yz, g_x_x_xxz_yzz, g_x_x_xxz_zz, g_x_x_xxz_zzz, g_x_x_xxzz_xx, g_x_x_xxzz_xy, g_x_x_xxzz_xz, g_x_x_xxzz_yy, g_x_x_xxzz_yz, g_x_x_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxzz_xx[k] = -g_x_x_xxz_xx[k] * ab_z + g_x_x_xxz_xxz[k];

                g_x_x_xxzz_xy[k] = -g_x_x_xxz_xy[k] * ab_z + g_x_x_xxz_xyz[k];

                g_x_x_xxzz_xz[k] = -g_x_x_xxz_xz[k] * ab_z + g_x_x_xxz_xzz[k];

                g_x_x_xxzz_yy[k] = -g_x_x_xxz_yy[k] * ab_z + g_x_x_xxz_yyz[k];

                g_x_x_xxzz_yz[k] = -g_x_x_xxz_yz[k] * ab_z + g_x_x_xxz_yzz[k];

                g_x_x_xxzz_zz[k] = -g_x_x_xxz_zz[k] * ab_z + g_x_x_xxz_zzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyyy_xx = cbuffer.data(gd_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_xyyy_xy = cbuffer.data(gd_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_xyyy_xz = cbuffer.data(gd_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_xyyy_yy = cbuffer.data(gd_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_xyyy_yz = cbuffer.data(gd_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_xyyy_zz = cbuffer.data(gd_geom_11_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xyy_xx, g_x_x_xyy_xxy, g_x_x_xyy_xy, g_x_x_xyy_xyy, g_x_x_xyy_xyz, g_x_x_xyy_xz, g_x_x_xyy_yy, g_x_x_xyy_yyy, g_x_x_xyy_yyz, g_x_x_xyy_yz, g_x_x_xyy_yzz, g_x_x_xyy_zz, g_x_x_xyyy_xx, g_x_x_xyyy_xy, g_x_x_xyyy_xz, g_x_x_xyyy_yy, g_x_x_xyyy_yz, g_x_x_xyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyyy_xx[k] = -g_x_x_xyy_xx[k] * ab_y + g_x_x_xyy_xxy[k];

                g_x_x_xyyy_xy[k] = -g_x_x_xyy_xy[k] * ab_y + g_x_x_xyy_xyy[k];

                g_x_x_xyyy_xz[k] = -g_x_x_xyy_xz[k] * ab_y + g_x_x_xyy_xyz[k];

                g_x_x_xyyy_yy[k] = -g_x_x_xyy_yy[k] * ab_y + g_x_x_xyy_yyy[k];

                g_x_x_xyyy_yz[k] = -g_x_x_xyy_yz[k] * ab_y + g_x_x_xyy_yyz[k];

                g_x_x_xyyy_zz[k] = -g_x_x_xyy_zz[k] * ab_y + g_x_x_xyy_yzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyyz_xx = cbuffer.data(gd_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_xyyz_xy = cbuffer.data(gd_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_xyyz_xz = cbuffer.data(gd_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_x_xyyz_yy = cbuffer.data(gd_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_x_xyyz_yz = cbuffer.data(gd_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_x_xyyz_zz = cbuffer.data(gd_geom_11_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xyyz_xx, g_x_x_xyyz_xy, g_x_x_xyyz_xz, g_x_x_xyyz_yy, g_x_x_xyyz_yz, g_x_x_xyyz_zz, g_x_x_xyz_xx, g_x_x_xyz_xxy, g_x_x_xyz_xy, g_x_x_xyz_xyy, g_x_x_xyz_xyz, g_x_x_xyz_xz, g_x_x_xyz_yy, g_x_x_xyz_yyy, g_x_x_xyz_yyz, g_x_x_xyz_yz, g_x_x_xyz_yzz, g_x_x_xyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyyz_xx[k] = -g_x_x_xyz_xx[k] * ab_y + g_x_x_xyz_xxy[k];

                g_x_x_xyyz_xy[k] = -g_x_x_xyz_xy[k] * ab_y + g_x_x_xyz_xyy[k];

                g_x_x_xyyz_xz[k] = -g_x_x_xyz_xz[k] * ab_y + g_x_x_xyz_xyz[k];

                g_x_x_xyyz_yy[k] = -g_x_x_xyz_yy[k] * ab_y + g_x_x_xyz_yyy[k];

                g_x_x_xyyz_yz[k] = -g_x_x_xyz_yz[k] * ab_y + g_x_x_xyz_yyz[k];

                g_x_x_xyyz_zz[k] = -g_x_x_xyz_zz[k] * ab_y + g_x_x_xyz_yzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyzz_xx = cbuffer.data(gd_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_x_xyzz_xy = cbuffer.data(gd_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_x_xyzz_xz = cbuffer.data(gd_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_x_xyzz_yy = cbuffer.data(gd_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_x_xyzz_yz = cbuffer.data(gd_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_x_xyzz_zz = cbuffer.data(gd_geom_11_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xyzz_xx, g_x_x_xyzz_xy, g_x_x_xyzz_xz, g_x_x_xyzz_yy, g_x_x_xyzz_yz, g_x_x_xyzz_zz, g_x_x_xzz_xx, g_x_x_xzz_xxy, g_x_x_xzz_xy, g_x_x_xzz_xyy, g_x_x_xzz_xyz, g_x_x_xzz_xz, g_x_x_xzz_yy, g_x_x_xzz_yyy, g_x_x_xzz_yyz, g_x_x_xzz_yz, g_x_x_xzz_yzz, g_x_x_xzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyzz_xx[k] = -g_x_x_xzz_xx[k] * ab_y + g_x_x_xzz_xxy[k];

                g_x_x_xyzz_xy[k] = -g_x_x_xzz_xy[k] * ab_y + g_x_x_xzz_xyy[k];

                g_x_x_xyzz_xz[k] = -g_x_x_xzz_xz[k] * ab_y + g_x_x_xzz_xyz[k];

                g_x_x_xyzz_yy[k] = -g_x_x_xzz_yy[k] * ab_y + g_x_x_xzz_yyy[k];

                g_x_x_xyzz_yz[k] = -g_x_x_xzz_yz[k] * ab_y + g_x_x_xzz_yyz[k];

                g_x_x_xyzz_zz[k] = -g_x_x_xzz_zz[k] * ab_y + g_x_x_xzz_yzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_x_x_xzzz_xx = cbuffer.data(gd_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_x_xzzz_xy = cbuffer.data(gd_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_x_xzzz_xz = cbuffer.data(gd_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_x_xzzz_yy = cbuffer.data(gd_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_x_xzzz_yz = cbuffer.data(gd_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_x_xzzz_zz = cbuffer.data(gd_geom_11_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xzz_xx, g_x_x_xzz_xxz, g_x_x_xzz_xy, g_x_x_xzz_xyz, g_x_x_xzz_xz, g_x_x_xzz_xzz, g_x_x_xzz_yy, g_x_x_xzz_yyz, g_x_x_xzz_yz, g_x_x_xzz_yzz, g_x_x_xzz_zz, g_x_x_xzz_zzz, g_x_x_xzzz_xx, g_x_x_xzzz_xy, g_x_x_xzzz_xz, g_x_x_xzzz_yy, g_x_x_xzzz_yz, g_x_x_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xzzz_xx[k] = -g_x_x_xzz_xx[k] * ab_z + g_x_x_xzz_xxz[k];

                g_x_x_xzzz_xy[k] = -g_x_x_xzz_xy[k] * ab_z + g_x_x_xzz_xyz[k];

                g_x_x_xzzz_xz[k] = -g_x_x_xzz_xz[k] * ab_z + g_x_x_xzz_xzz[k];

                g_x_x_xzzz_yy[k] = -g_x_x_xzz_yy[k] * ab_z + g_x_x_xzz_yyz[k];

                g_x_x_xzzz_yz[k] = -g_x_x_xzz_yz[k] * ab_z + g_x_x_xzz_yzz[k];

                g_x_x_xzzz_zz[k] = -g_x_x_xzz_zz[k] * ab_z + g_x_x_xzz_zzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyyy_xx = cbuffer.data(gd_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_x_yyyy_xy = cbuffer.data(gd_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_x_yyyy_xz = cbuffer.data(gd_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_x_yyyy_yy = cbuffer.data(gd_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_x_yyyy_yz = cbuffer.data(gd_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_x_yyyy_zz = cbuffer.data(gd_geom_11_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yyy_xx, g_x_x_yyy_xxy, g_x_x_yyy_xy, g_x_x_yyy_xyy, g_x_x_yyy_xyz, g_x_x_yyy_xz, g_x_x_yyy_yy, g_x_x_yyy_yyy, g_x_x_yyy_yyz, g_x_x_yyy_yz, g_x_x_yyy_yzz, g_x_x_yyy_zz, g_x_x_yyyy_xx, g_x_x_yyyy_xy, g_x_x_yyyy_xz, g_x_x_yyyy_yy, g_x_x_yyyy_yz, g_x_x_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyyy_xx[k] = -g_x_x_yyy_xx[k] * ab_y + g_x_x_yyy_xxy[k];

                g_x_x_yyyy_xy[k] = -g_x_x_yyy_xy[k] * ab_y + g_x_x_yyy_xyy[k];

                g_x_x_yyyy_xz[k] = -g_x_x_yyy_xz[k] * ab_y + g_x_x_yyy_xyz[k];

                g_x_x_yyyy_yy[k] = -g_x_x_yyy_yy[k] * ab_y + g_x_x_yyy_yyy[k];

                g_x_x_yyyy_yz[k] = -g_x_x_yyy_yz[k] * ab_y + g_x_x_yyy_yyz[k];

                g_x_x_yyyy_zz[k] = -g_x_x_yyy_zz[k] * ab_y + g_x_x_yyy_yzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyyz_xx = cbuffer.data(gd_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_x_yyyz_xy = cbuffer.data(gd_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_x_yyyz_xz = cbuffer.data(gd_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_x_yyyz_yy = cbuffer.data(gd_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_x_yyyz_yz = cbuffer.data(gd_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_x_yyyz_zz = cbuffer.data(gd_geom_11_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yyyz_xx, g_x_x_yyyz_xy, g_x_x_yyyz_xz, g_x_x_yyyz_yy, g_x_x_yyyz_yz, g_x_x_yyyz_zz, g_x_x_yyz_xx, g_x_x_yyz_xxy, g_x_x_yyz_xy, g_x_x_yyz_xyy, g_x_x_yyz_xyz, g_x_x_yyz_xz, g_x_x_yyz_yy, g_x_x_yyz_yyy, g_x_x_yyz_yyz, g_x_x_yyz_yz, g_x_x_yyz_yzz, g_x_x_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyyz_xx[k] = -g_x_x_yyz_xx[k] * ab_y + g_x_x_yyz_xxy[k];

                g_x_x_yyyz_xy[k] = -g_x_x_yyz_xy[k] * ab_y + g_x_x_yyz_xyy[k];

                g_x_x_yyyz_xz[k] = -g_x_x_yyz_xz[k] * ab_y + g_x_x_yyz_xyz[k];

                g_x_x_yyyz_yy[k] = -g_x_x_yyz_yy[k] * ab_y + g_x_x_yyz_yyy[k];

                g_x_x_yyyz_yz[k] = -g_x_x_yyz_yz[k] * ab_y + g_x_x_yyz_yyz[k];

                g_x_x_yyyz_zz[k] = -g_x_x_yyz_zz[k] * ab_y + g_x_x_yyz_yzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyzz_xx = cbuffer.data(gd_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_x_yyzz_xy = cbuffer.data(gd_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_x_yyzz_xz = cbuffer.data(gd_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_x_yyzz_yy = cbuffer.data(gd_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_x_yyzz_yz = cbuffer.data(gd_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_x_yyzz_zz = cbuffer.data(gd_geom_11_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yyzz_xx, g_x_x_yyzz_xy, g_x_x_yyzz_xz, g_x_x_yyzz_yy, g_x_x_yyzz_yz, g_x_x_yyzz_zz, g_x_x_yzz_xx, g_x_x_yzz_xxy, g_x_x_yzz_xy, g_x_x_yzz_xyy, g_x_x_yzz_xyz, g_x_x_yzz_xz, g_x_x_yzz_yy, g_x_x_yzz_yyy, g_x_x_yzz_yyz, g_x_x_yzz_yz, g_x_x_yzz_yzz, g_x_x_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyzz_xx[k] = -g_x_x_yzz_xx[k] * ab_y + g_x_x_yzz_xxy[k];

                g_x_x_yyzz_xy[k] = -g_x_x_yzz_xy[k] * ab_y + g_x_x_yzz_xyy[k];

                g_x_x_yyzz_xz[k] = -g_x_x_yzz_xz[k] * ab_y + g_x_x_yzz_xyz[k];

                g_x_x_yyzz_yy[k] = -g_x_x_yzz_yy[k] * ab_y + g_x_x_yzz_yyy[k];

                g_x_x_yyzz_yz[k] = -g_x_x_yzz_yz[k] * ab_y + g_x_x_yzz_yyz[k];

                g_x_x_yyzz_zz[k] = -g_x_x_yzz_zz[k] * ab_y + g_x_x_yzz_yzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_x_x_yzzz_xx = cbuffer.data(gd_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_x_yzzz_xy = cbuffer.data(gd_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_x_yzzz_xz = cbuffer.data(gd_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_x_yzzz_yy = cbuffer.data(gd_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_x_yzzz_yz = cbuffer.data(gd_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_x_yzzz_zz = cbuffer.data(gd_geom_11_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yzzz_xx, g_x_x_yzzz_xy, g_x_x_yzzz_xz, g_x_x_yzzz_yy, g_x_x_yzzz_yz, g_x_x_yzzz_zz, g_x_x_zzz_xx, g_x_x_zzz_xxy, g_x_x_zzz_xy, g_x_x_zzz_xyy, g_x_x_zzz_xyz, g_x_x_zzz_xz, g_x_x_zzz_yy, g_x_x_zzz_yyy, g_x_x_zzz_yyz, g_x_x_zzz_yz, g_x_x_zzz_yzz, g_x_x_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yzzz_xx[k] = -g_x_x_zzz_xx[k] * ab_y + g_x_x_zzz_xxy[k];

                g_x_x_yzzz_xy[k] = -g_x_x_zzz_xy[k] * ab_y + g_x_x_zzz_xyy[k];

                g_x_x_yzzz_xz[k] = -g_x_x_zzz_xz[k] * ab_y + g_x_x_zzz_xyz[k];

                g_x_x_yzzz_yy[k] = -g_x_x_zzz_yy[k] * ab_y + g_x_x_zzz_yyy[k];

                g_x_x_yzzz_yz[k] = -g_x_x_zzz_yz[k] * ab_y + g_x_x_zzz_yyz[k];

                g_x_x_yzzz_zz[k] = -g_x_x_zzz_zz[k] * ab_y + g_x_x_zzz_yzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_x_x_zzzz_xx = cbuffer.data(gd_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_x_zzzz_xy = cbuffer.data(gd_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_x_zzzz_xz = cbuffer.data(gd_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_x_zzzz_yy = cbuffer.data(gd_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_x_zzzz_yz = cbuffer.data(gd_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_x_zzzz_zz = cbuffer.data(gd_geom_11_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_zzz_xx, g_x_x_zzz_xxz, g_x_x_zzz_xy, g_x_x_zzz_xyz, g_x_x_zzz_xz, g_x_x_zzz_xzz, g_x_x_zzz_yy, g_x_x_zzz_yyz, g_x_x_zzz_yz, g_x_x_zzz_yzz, g_x_x_zzz_zz, g_x_x_zzz_zzz, g_x_x_zzzz_xx, g_x_x_zzzz_xy, g_x_x_zzzz_xz, g_x_x_zzzz_yy, g_x_x_zzzz_yz, g_x_x_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_zzzz_xx[k] = -g_x_x_zzz_xx[k] * ab_z + g_x_x_zzz_xxz[k];

                g_x_x_zzzz_xy[k] = -g_x_x_zzz_xy[k] * ab_z + g_x_x_zzz_xyz[k];

                g_x_x_zzzz_xz[k] = -g_x_x_zzz_xz[k] * ab_z + g_x_x_zzz_xzz[k];

                g_x_x_zzzz_yy[k] = -g_x_x_zzz_yy[k] * ab_z + g_x_x_zzz_yyz[k];

                g_x_x_zzzz_yz[k] = -g_x_x_zzz_yz[k] * ab_z + g_x_x_zzz_yzz[k];

                g_x_x_zzzz_zz[k] = -g_x_x_zzz_zz[k] * ab_z + g_x_x_zzz_zzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxxx_xx = cbuffer.data(gd_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_y_xxxx_xy = cbuffer.data(gd_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_y_xxxx_xz = cbuffer.data(gd_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_y_xxxx_yy = cbuffer.data(gd_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_y_xxxx_yz = cbuffer.data(gd_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_y_xxxx_zz = cbuffer.data(gd_geom_11_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxx_xx, g_0_y_xxx_xy, g_0_y_xxx_xz, g_0_y_xxx_yy, g_0_y_xxx_yz, g_0_y_xxx_zz, g_x_y_xxx_xx, g_x_y_xxx_xxx, g_x_y_xxx_xxy, g_x_y_xxx_xxz, g_x_y_xxx_xy, g_x_y_xxx_xyy, g_x_y_xxx_xyz, g_x_y_xxx_xz, g_x_y_xxx_xzz, g_x_y_xxx_yy, g_x_y_xxx_yz, g_x_y_xxx_zz, g_x_y_xxxx_xx, g_x_y_xxxx_xy, g_x_y_xxxx_xz, g_x_y_xxxx_yy, g_x_y_xxxx_yz, g_x_y_xxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxxx_xx[k] = -g_0_y_xxx_xx[k] - g_x_y_xxx_xx[k] * ab_x + g_x_y_xxx_xxx[k];

                g_x_y_xxxx_xy[k] = -g_0_y_xxx_xy[k] - g_x_y_xxx_xy[k] * ab_x + g_x_y_xxx_xxy[k];

                g_x_y_xxxx_xz[k] = -g_0_y_xxx_xz[k] - g_x_y_xxx_xz[k] * ab_x + g_x_y_xxx_xxz[k];

                g_x_y_xxxx_yy[k] = -g_0_y_xxx_yy[k] - g_x_y_xxx_yy[k] * ab_x + g_x_y_xxx_xyy[k];

                g_x_y_xxxx_yz[k] = -g_0_y_xxx_yz[k] - g_x_y_xxx_yz[k] * ab_x + g_x_y_xxx_xyz[k];

                g_x_y_xxxx_zz[k] = -g_0_y_xxx_zz[k] - g_x_y_xxx_zz[k] * ab_x + g_x_y_xxx_xzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxxy_xx = cbuffer.data(gd_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_y_xxxy_xy = cbuffer.data(gd_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_y_xxxy_xz = cbuffer.data(gd_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_y_xxxy_yy = cbuffer.data(gd_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_y_xxxy_yz = cbuffer.data(gd_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_y_xxxy_zz = cbuffer.data(gd_geom_11_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxy_xx, g_0_y_xxy_xy, g_0_y_xxy_xz, g_0_y_xxy_yy, g_0_y_xxy_yz, g_0_y_xxy_zz, g_x_y_xxxy_xx, g_x_y_xxxy_xy, g_x_y_xxxy_xz, g_x_y_xxxy_yy, g_x_y_xxxy_yz, g_x_y_xxxy_zz, g_x_y_xxy_xx, g_x_y_xxy_xxx, g_x_y_xxy_xxy, g_x_y_xxy_xxz, g_x_y_xxy_xy, g_x_y_xxy_xyy, g_x_y_xxy_xyz, g_x_y_xxy_xz, g_x_y_xxy_xzz, g_x_y_xxy_yy, g_x_y_xxy_yz, g_x_y_xxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxxy_xx[k] = -g_0_y_xxy_xx[k] - g_x_y_xxy_xx[k] * ab_x + g_x_y_xxy_xxx[k];

                g_x_y_xxxy_xy[k] = -g_0_y_xxy_xy[k] - g_x_y_xxy_xy[k] * ab_x + g_x_y_xxy_xxy[k];

                g_x_y_xxxy_xz[k] = -g_0_y_xxy_xz[k] - g_x_y_xxy_xz[k] * ab_x + g_x_y_xxy_xxz[k];

                g_x_y_xxxy_yy[k] = -g_0_y_xxy_yy[k] - g_x_y_xxy_yy[k] * ab_x + g_x_y_xxy_xyy[k];

                g_x_y_xxxy_yz[k] = -g_0_y_xxy_yz[k] - g_x_y_xxy_yz[k] * ab_x + g_x_y_xxy_xyz[k];

                g_x_y_xxxy_zz[k] = -g_0_y_xxy_zz[k] - g_x_y_xxy_zz[k] * ab_x + g_x_y_xxy_xzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxxz_xx = cbuffer.data(gd_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_y_xxxz_xy = cbuffer.data(gd_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_y_xxxz_xz = cbuffer.data(gd_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_y_xxxz_yy = cbuffer.data(gd_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_y_xxxz_yz = cbuffer.data(gd_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_y_xxxz_zz = cbuffer.data(gd_geom_11_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xxx_xx, g_x_y_xxx_xxz, g_x_y_xxx_xy, g_x_y_xxx_xyz, g_x_y_xxx_xz, g_x_y_xxx_xzz, g_x_y_xxx_yy, g_x_y_xxx_yyz, g_x_y_xxx_yz, g_x_y_xxx_yzz, g_x_y_xxx_zz, g_x_y_xxx_zzz, g_x_y_xxxz_xx, g_x_y_xxxz_xy, g_x_y_xxxz_xz, g_x_y_xxxz_yy, g_x_y_xxxz_yz, g_x_y_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxxz_xx[k] = -g_x_y_xxx_xx[k] * ab_z + g_x_y_xxx_xxz[k];

                g_x_y_xxxz_xy[k] = -g_x_y_xxx_xy[k] * ab_z + g_x_y_xxx_xyz[k];

                g_x_y_xxxz_xz[k] = -g_x_y_xxx_xz[k] * ab_z + g_x_y_xxx_xzz[k];

                g_x_y_xxxz_yy[k] = -g_x_y_xxx_yy[k] * ab_z + g_x_y_xxx_yyz[k];

                g_x_y_xxxz_yz[k] = -g_x_y_xxx_yz[k] * ab_z + g_x_y_xxx_yzz[k];

                g_x_y_xxxz_zz[k] = -g_x_y_xxx_zz[k] * ab_z + g_x_y_xxx_zzz[k];
            }

            /// Set up 108-114 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxyy_xx = cbuffer.data(gd_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_y_xxyy_xy = cbuffer.data(gd_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_y_xxyy_xz = cbuffer.data(gd_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_y_xxyy_yy = cbuffer.data(gd_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_y_xxyy_yz = cbuffer.data(gd_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_y_xxyy_zz = cbuffer.data(gd_geom_11_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyy_xx, g_0_y_xyy_xy, g_0_y_xyy_xz, g_0_y_xyy_yy, g_0_y_xyy_yz, g_0_y_xyy_zz, g_x_y_xxyy_xx, g_x_y_xxyy_xy, g_x_y_xxyy_xz, g_x_y_xxyy_yy, g_x_y_xxyy_yz, g_x_y_xxyy_zz, g_x_y_xyy_xx, g_x_y_xyy_xxx, g_x_y_xyy_xxy, g_x_y_xyy_xxz, g_x_y_xyy_xy, g_x_y_xyy_xyy, g_x_y_xyy_xyz, g_x_y_xyy_xz, g_x_y_xyy_xzz, g_x_y_xyy_yy, g_x_y_xyy_yz, g_x_y_xyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxyy_xx[k] = -g_0_y_xyy_xx[k] - g_x_y_xyy_xx[k] * ab_x + g_x_y_xyy_xxx[k];

                g_x_y_xxyy_xy[k] = -g_0_y_xyy_xy[k] - g_x_y_xyy_xy[k] * ab_x + g_x_y_xyy_xxy[k];

                g_x_y_xxyy_xz[k] = -g_0_y_xyy_xz[k] - g_x_y_xyy_xz[k] * ab_x + g_x_y_xyy_xxz[k];

                g_x_y_xxyy_yy[k] = -g_0_y_xyy_yy[k] - g_x_y_xyy_yy[k] * ab_x + g_x_y_xyy_xyy[k];

                g_x_y_xxyy_yz[k] = -g_0_y_xyy_yz[k] - g_x_y_xyy_yz[k] * ab_x + g_x_y_xyy_xyz[k];

                g_x_y_xxyy_zz[k] = -g_0_y_xyy_zz[k] - g_x_y_xyy_zz[k] * ab_x + g_x_y_xyy_xzz[k];
            }

            /// Set up 114-120 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxyz_xx = cbuffer.data(gd_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_y_xxyz_xy = cbuffer.data(gd_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_y_xxyz_xz = cbuffer.data(gd_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_y_xxyz_yy = cbuffer.data(gd_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_y_xxyz_yz = cbuffer.data(gd_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_y_xxyz_zz = cbuffer.data(gd_geom_11_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xxy_xx, g_x_y_xxy_xxz, g_x_y_xxy_xy, g_x_y_xxy_xyz, g_x_y_xxy_xz, g_x_y_xxy_xzz, g_x_y_xxy_yy, g_x_y_xxy_yyz, g_x_y_xxy_yz, g_x_y_xxy_yzz, g_x_y_xxy_zz, g_x_y_xxy_zzz, g_x_y_xxyz_xx, g_x_y_xxyz_xy, g_x_y_xxyz_xz, g_x_y_xxyz_yy, g_x_y_xxyz_yz, g_x_y_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxyz_xx[k] = -g_x_y_xxy_xx[k] * ab_z + g_x_y_xxy_xxz[k];

                g_x_y_xxyz_xy[k] = -g_x_y_xxy_xy[k] * ab_z + g_x_y_xxy_xyz[k];

                g_x_y_xxyz_xz[k] = -g_x_y_xxy_xz[k] * ab_z + g_x_y_xxy_xzz[k];

                g_x_y_xxyz_yy[k] = -g_x_y_xxy_yy[k] * ab_z + g_x_y_xxy_yyz[k];

                g_x_y_xxyz_yz[k] = -g_x_y_xxy_yz[k] * ab_z + g_x_y_xxy_yzz[k];

                g_x_y_xxyz_zz[k] = -g_x_y_xxy_zz[k] * ab_z + g_x_y_xxy_zzz[k];
            }

            /// Set up 120-126 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxzz_xx = cbuffer.data(gd_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_y_xxzz_xy = cbuffer.data(gd_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_y_xxzz_xz = cbuffer.data(gd_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_y_xxzz_yy = cbuffer.data(gd_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_y_xxzz_yz = cbuffer.data(gd_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_y_xxzz_zz = cbuffer.data(gd_geom_11_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xxz_xx, g_x_y_xxz_xxz, g_x_y_xxz_xy, g_x_y_xxz_xyz, g_x_y_xxz_xz, g_x_y_xxz_xzz, g_x_y_xxz_yy, g_x_y_xxz_yyz, g_x_y_xxz_yz, g_x_y_xxz_yzz, g_x_y_xxz_zz, g_x_y_xxz_zzz, g_x_y_xxzz_xx, g_x_y_xxzz_xy, g_x_y_xxzz_xz, g_x_y_xxzz_yy, g_x_y_xxzz_yz, g_x_y_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxzz_xx[k] = -g_x_y_xxz_xx[k] * ab_z + g_x_y_xxz_xxz[k];

                g_x_y_xxzz_xy[k] = -g_x_y_xxz_xy[k] * ab_z + g_x_y_xxz_xyz[k];

                g_x_y_xxzz_xz[k] = -g_x_y_xxz_xz[k] * ab_z + g_x_y_xxz_xzz[k];

                g_x_y_xxzz_yy[k] = -g_x_y_xxz_yy[k] * ab_z + g_x_y_xxz_yyz[k];

                g_x_y_xxzz_yz[k] = -g_x_y_xxz_yz[k] * ab_z + g_x_y_xxz_yzz[k];

                g_x_y_xxzz_zz[k] = -g_x_y_xxz_zz[k] * ab_z + g_x_y_xxz_zzz[k];
            }

            /// Set up 126-132 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyyy_xx = cbuffer.data(gd_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_y_xyyy_xy = cbuffer.data(gd_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_y_xyyy_xz = cbuffer.data(gd_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_y_xyyy_yy = cbuffer.data(gd_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_y_xyyy_yz = cbuffer.data(gd_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_y_xyyy_zz = cbuffer.data(gd_geom_11_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyy_xx, g_0_y_yyy_xy, g_0_y_yyy_xz, g_0_y_yyy_yy, g_0_y_yyy_yz, g_0_y_yyy_zz, g_x_y_xyyy_xx, g_x_y_xyyy_xy, g_x_y_xyyy_xz, g_x_y_xyyy_yy, g_x_y_xyyy_yz, g_x_y_xyyy_zz, g_x_y_yyy_xx, g_x_y_yyy_xxx, g_x_y_yyy_xxy, g_x_y_yyy_xxz, g_x_y_yyy_xy, g_x_y_yyy_xyy, g_x_y_yyy_xyz, g_x_y_yyy_xz, g_x_y_yyy_xzz, g_x_y_yyy_yy, g_x_y_yyy_yz, g_x_y_yyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyyy_xx[k] = -g_0_y_yyy_xx[k] - g_x_y_yyy_xx[k] * ab_x + g_x_y_yyy_xxx[k];

                g_x_y_xyyy_xy[k] = -g_0_y_yyy_xy[k] - g_x_y_yyy_xy[k] * ab_x + g_x_y_yyy_xxy[k];

                g_x_y_xyyy_xz[k] = -g_0_y_yyy_xz[k] - g_x_y_yyy_xz[k] * ab_x + g_x_y_yyy_xxz[k];

                g_x_y_xyyy_yy[k] = -g_0_y_yyy_yy[k] - g_x_y_yyy_yy[k] * ab_x + g_x_y_yyy_xyy[k];

                g_x_y_xyyy_yz[k] = -g_0_y_yyy_yz[k] - g_x_y_yyy_yz[k] * ab_x + g_x_y_yyy_xyz[k];

                g_x_y_xyyy_zz[k] = -g_0_y_yyy_zz[k] - g_x_y_yyy_zz[k] * ab_x + g_x_y_yyy_xzz[k];
            }

            /// Set up 132-138 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyyz_xx = cbuffer.data(gd_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_y_xyyz_xy = cbuffer.data(gd_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_y_xyyz_xz = cbuffer.data(gd_geom_11_off + 134 * ccomps * dcomps);

            auto g_x_y_xyyz_yy = cbuffer.data(gd_geom_11_off + 135 * ccomps * dcomps);

            auto g_x_y_xyyz_yz = cbuffer.data(gd_geom_11_off + 136 * ccomps * dcomps);

            auto g_x_y_xyyz_zz = cbuffer.data(gd_geom_11_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xyy_xx, g_x_y_xyy_xxz, g_x_y_xyy_xy, g_x_y_xyy_xyz, g_x_y_xyy_xz, g_x_y_xyy_xzz, g_x_y_xyy_yy, g_x_y_xyy_yyz, g_x_y_xyy_yz, g_x_y_xyy_yzz, g_x_y_xyy_zz, g_x_y_xyy_zzz, g_x_y_xyyz_xx, g_x_y_xyyz_xy, g_x_y_xyyz_xz, g_x_y_xyyz_yy, g_x_y_xyyz_yz, g_x_y_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyyz_xx[k] = -g_x_y_xyy_xx[k] * ab_z + g_x_y_xyy_xxz[k];

                g_x_y_xyyz_xy[k] = -g_x_y_xyy_xy[k] * ab_z + g_x_y_xyy_xyz[k];

                g_x_y_xyyz_xz[k] = -g_x_y_xyy_xz[k] * ab_z + g_x_y_xyy_xzz[k];

                g_x_y_xyyz_yy[k] = -g_x_y_xyy_yy[k] * ab_z + g_x_y_xyy_yyz[k];

                g_x_y_xyyz_yz[k] = -g_x_y_xyy_yz[k] * ab_z + g_x_y_xyy_yzz[k];

                g_x_y_xyyz_zz[k] = -g_x_y_xyy_zz[k] * ab_z + g_x_y_xyy_zzz[k];
            }

            /// Set up 138-144 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyzz_xx = cbuffer.data(gd_geom_11_off + 138 * ccomps * dcomps);

            auto g_x_y_xyzz_xy = cbuffer.data(gd_geom_11_off + 139 * ccomps * dcomps);

            auto g_x_y_xyzz_xz = cbuffer.data(gd_geom_11_off + 140 * ccomps * dcomps);

            auto g_x_y_xyzz_yy = cbuffer.data(gd_geom_11_off + 141 * ccomps * dcomps);

            auto g_x_y_xyzz_yz = cbuffer.data(gd_geom_11_off + 142 * ccomps * dcomps);

            auto g_x_y_xyzz_zz = cbuffer.data(gd_geom_11_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xyz_xx, g_x_y_xyz_xxz, g_x_y_xyz_xy, g_x_y_xyz_xyz, g_x_y_xyz_xz, g_x_y_xyz_xzz, g_x_y_xyz_yy, g_x_y_xyz_yyz, g_x_y_xyz_yz, g_x_y_xyz_yzz, g_x_y_xyz_zz, g_x_y_xyz_zzz, g_x_y_xyzz_xx, g_x_y_xyzz_xy, g_x_y_xyzz_xz, g_x_y_xyzz_yy, g_x_y_xyzz_yz, g_x_y_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyzz_xx[k] = -g_x_y_xyz_xx[k] * ab_z + g_x_y_xyz_xxz[k];

                g_x_y_xyzz_xy[k] = -g_x_y_xyz_xy[k] * ab_z + g_x_y_xyz_xyz[k];

                g_x_y_xyzz_xz[k] = -g_x_y_xyz_xz[k] * ab_z + g_x_y_xyz_xzz[k];

                g_x_y_xyzz_yy[k] = -g_x_y_xyz_yy[k] * ab_z + g_x_y_xyz_yyz[k];

                g_x_y_xyzz_yz[k] = -g_x_y_xyz_yz[k] * ab_z + g_x_y_xyz_yzz[k];

                g_x_y_xyzz_zz[k] = -g_x_y_xyz_zz[k] * ab_z + g_x_y_xyz_zzz[k];
            }

            /// Set up 144-150 components of targeted buffer : cbuffer.data(

            auto g_x_y_xzzz_xx = cbuffer.data(gd_geom_11_off + 144 * ccomps * dcomps);

            auto g_x_y_xzzz_xy = cbuffer.data(gd_geom_11_off + 145 * ccomps * dcomps);

            auto g_x_y_xzzz_xz = cbuffer.data(gd_geom_11_off + 146 * ccomps * dcomps);

            auto g_x_y_xzzz_yy = cbuffer.data(gd_geom_11_off + 147 * ccomps * dcomps);

            auto g_x_y_xzzz_yz = cbuffer.data(gd_geom_11_off + 148 * ccomps * dcomps);

            auto g_x_y_xzzz_zz = cbuffer.data(gd_geom_11_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xzz_xx, g_x_y_xzz_xxz, g_x_y_xzz_xy, g_x_y_xzz_xyz, g_x_y_xzz_xz, g_x_y_xzz_xzz, g_x_y_xzz_yy, g_x_y_xzz_yyz, g_x_y_xzz_yz, g_x_y_xzz_yzz, g_x_y_xzz_zz, g_x_y_xzz_zzz, g_x_y_xzzz_xx, g_x_y_xzzz_xy, g_x_y_xzzz_xz, g_x_y_xzzz_yy, g_x_y_xzzz_yz, g_x_y_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xzzz_xx[k] = -g_x_y_xzz_xx[k] * ab_z + g_x_y_xzz_xxz[k];

                g_x_y_xzzz_xy[k] = -g_x_y_xzz_xy[k] * ab_z + g_x_y_xzz_xyz[k];

                g_x_y_xzzz_xz[k] = -g_x_y_xzz_xz[k] * ab_z + g_x_y_xzz_xzz[k];

                g_x_y_xzzz_yy[k] = -g_x_y_xzz_yy[k] * ab_z + g_x_y_xzz_yyz[k];

                g_x_y_xzzz_yz[k] = -g_x_y_xzz_yz[k] * ab_z + g_x_y_xzz_yzz[k];

                g_x_y_xzzz_zz[k] = -g_x_y_xzz_zz[k] * ab_z + g_x_y_xzz_zzz[k];
            }

            /// Set up 150-156 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyyy_xx = cbuffer.data(gd_geom_11_off + 150 * ccomps * dcomps);

            auto g_x_y_yyyy_xy = cbuffer.data(gd_geom_11_off + 151 * ccomps * dcomps);

            auto g_x_y_yyyy_xz = cbuffer.data(gd_geom_11_off + 152 * ccomps * dcomps);

            auto g_x_y_yyyy_yy = cbuffer.data(gd_geom_11_off + 153 * ccomps * dcomps);

            auto g_x_y_yyyy_yz = cbuffer.data(gd_geom_11_off + 154 * ccomps * dcomps);

            auto g_x_y_yyyy_zz = cbuffer.data(gd_geom_11_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyy_xx, g_x_0_yyy_xy, g_x_0_yyy_xz, g_x_0_yyy_yy, g_x_0_yyy_yz, g_x_0_yyy_zz, g_x_y_yyy_xx, g_x_y_yyy_xxy, g_x_y_yyy_xy, g_x_y_yyy_xyy, g_x_y_yyy_xyz, g_x_y_yyy_xz, g_x_y_yyy_yy, g_x_y_yyy_yyy, g_x_y_yyy_yyz, g_x_y_yyy_yz, g_x_y_yyy_yzz, g_x_y_yyy_zz, g_x_y_yyyy_xx, g_x_y_yyyy_xy, g_x_y_yyyy_xz, g_x_y_yyyy_yy, g_x_y_yyyy_yz, g_x_y_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyyy_xx[k] = g_x_0_yyy_xx[k] - g_x_y_yyy_xx[k] * ab_y + g_x_y_yyy_xxy[k];

                g_x_y_yyyy_xy[k] = g_x_0_yyy_xy[k] - g_x_y_yyy_xy[k] * ab_y + g_x_y_yyy_xyy[k];

                g_x_y_yyyy_xz[k] = g_x_0_yyy_xz[k] - g_x_y_yyy_xz[k] * ab_y + g_x_y_yyy_xyz[k];

                g_x_y_yyyy_yy[k] = g_x_0_yyy_yy[k] - g_x_y_yyy_yy[k] * ab_y + g_x_y_yyy_yyy[k];

                g_x_y_yyyy_yz[k] = g_x_0_yyy_yz[k] - g_x_y_yyy_yz[k] * ab_y + g_x_y_yyy_yyz[k];

                g_x_y_yyyy_zz[k] = g_x_0_yyy_zz[k] - g_x_y_yyy_zz[k] * ab_y + g_x_y_yyy_yzz[k];
            }

            /// Set up 156-162 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyyz_xx = cbuffer.data(gd_geom_11_off + 156 * ccomps * dcomps);

            auto g_x_y_yyyz_xy = cbuffer.data(gd_geom_11_off + 157 * ccomps * dcomps);

            auto g_x_y_yyyz_xz = cbuffer.data(gd_geom_11_off + 158 * ccomps * dcomps);

            auto g_x_y_yyyz_yy = cbuffer.data(gd_geom_11_off + 159 * ccomps * dcomps);

            auto g_x_y_yyyz_yz = cbuffer.data(gd_geom_11_off + 160 * ccomps * dcomps);

            auto g_x_y_yyyz_zz = cbuffer.data(gd_geom_11_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yyy_xx, g_x_y_yyy_xxz, g_x_y_yyy_xy, g_x_y_yyy_xyz, g_x_y_yyy_xz, g_x_y_yyy_xzz, g_x_y_yyy_yy, g_x_y_yyy_yyz, g_x_y_yyy_yz, g_x_y_yyy_yzz, g_x_y_yyy_zz, g_x_y_yyy_zzz, g_x_y_yyyz_xx, g_x_y_yyyz_xy, g_x_y_yyyz_xz, g_x_y_yyyz_yy, g_x_y_yyyz_yz, g_x_y_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyyz_xx[k] = -g_x_y_yyy_xx[k] * ab_z + g_x_y_yyy_xxz[k];

                g_x_y_yyyz_xy[k] = -g_x_y_yyy_xy[k] * ab_z + g_x_y_yyy_xyz[k];

                g_x_y_yyyz_xz[k] = -g_x_y_yyy_xz[k] * ab_z + g_x_y_yyy_xzz[k];

                g_x_y_yyyz_yy[k] = -g_x_y_yyy_yy[k] * ab_z + g_x_y_yyy_yyz[k];

                g_x_y_yyyz_yz[k] = -g_x_y_yyy_yz[k] * ab_z + g_x_y_yyy_yzz[k];

                g_x_y_yyyz_zz[k] = -g_x_y_yyy_zz[k] * ab_z + g_x_y_yyy_zzz[k];
            }

            /// Set up 162-168 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyzz_xx = cbuffer.data(gd_geom_11_off + 162 * ccomps * dcomps);

            auto g_x_y_yyzz_xy = cbuffer.data(gd_geom_11_off + 163 * ccomps * dcomps);

            auto g_x_y_yyzz_xz = cbuffer.data(gd_geom_11_off + 164 * ccomps * dcomps);

            auto g_x_y_yyzz_yy = cbuffer.data(gd_geom_11_off + 165 * ccomps * dcomps);

            auto g_x_y_yyzz_yz = cbuffer.data(gd_geom_11_off + 166 * ccomps * dcomps);

            auto g_x_y_yyzz_zz = cbuffer.data(gd_geom_11_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yyz_xx, g_x_y_yyz_xxz, g_x_y_yyz_xy, g_x_y_yyz_xyz, g_x_y_yyz_xz, g_x_y_yyz_xzz, g_x_y_yyz_yy, g_x_y_yyz_yyz, g_x_y_yyz_yz, g_x_y_yyz_yzz, g_x_y_yyz_zz, g_x_y_yyz_zzz, g_x_y_yyzz_xx, g_x_y_yyzz_xy, g_x_y_yyzz_xz, g_x_y_yyzz_yy, g_x_y_yyzz_yz, g_x_y_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyzz_xx[k] = -g_x_y_yyz_xx[k] * ab_z + g_x_y_yyz_xxz[k];

                g_x_y_yyzz_xy[k] = -g_x_y_yyz_xy[k] * ab_z + g_x_y_yyz_xyz[k];

                g_x_y_yyzz_xz[k] = -g_x_y_yyz_xz[k] * ab_z + g_x_y_yyz_xzz[k];

                g_x_y_yyzz_yy[k] = -g_x_y_yyz_yy[k] * ab_z + g_x_y_yyz_yyz[k];

                g_x_y_yyzz_yz[k] = -g_x_y_yyz_yz[k] * ab_z + g_x_y_yyz_yzz[k];

                g_x_y_yyzz_zz[k] = -g_x_y_yyz_zz[k] * ab_z + g_x_y_yyz_zzz[k];
            }

            /// Set up 168-174 components of targeted buffer : cbuffer.data(

            auto g_x_y_yzzz_xx = cbuffer.data(gd_geom_11_off + 168 * ccomps * dcomps);

            auto g_x_y_yzzz_xy = cbuffer.data(gd_geom_11_off + 169 * ccomps * dcomps);

            auto g_x_y_yzzz_xz = cbuffer.data(gd_geom_11_off + 170 * ccomps * dcomps);

            auto g_x_y_yzzz_yy = cbuffer.data(gd_geom_11_off + 171 * ccomps * dcomps);

            auto g_x_y_yzzz_yz = cbuffer.data(gd_geom_11_off + 172 * ccomps * dcomps);

            auto g_x_y_yzzz_zz = cbuffer.data(gd_geom_11_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yzz_xx, g_x_y_yzz_xxz, g_x_y_yzz_xy, g_x_y_yzz_xyz, g_x_y_yzz_xz, g_x_y_yzz_xzz, g_x_y_yzz_yy, g_x_y_yzz_yyz, g_x_y_yzz_yz, g_x_y_yzz_yzz, g_x_y_yzz_zz, g_x_y_yzz_zzz, g_x_y_yzzz_xx, g_x_y_yzzz_xy, g_x_y_yzzz_xz, g_x_y_yzzz_yy, g_x_y_yzzz_yz, g_x_y_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yzzz_xx[k] = -g_x_y_yzz_xx[k] * ab_z + g_x_y_yzz_xxz[k];

                g_x_y_yzzz_xy[k] = -g_x_y_yzz_xy[k] * ab_z + g_x_y_yzz_xyz[k];

                g_x_y_yzzz_xz[k] = -g_x_y_yzz_xz[k] * ab_z + g_x_y_yzz_xzz[k];

                g_x_y_yzzz_yy[k] = -g_x_y_yzz_yy[k] * ab_z + g_x_y_yzz_yyz[k];

                g_x_y_yzzz_yz[k] = -g_x_y_yzz_yz[k] * ab_z + g_x_y_yzz_yzz[k];

                g_x_y_yzzz_zz[k] = -g_x_y_yzz_zz[k] * ab_z + g_x_y_yzz_zzz[k];
            }

            /// Set up 174-180 components of targeted buffer : cbuffer.data(

            auto g_x_y_zzzz_xx = cbuffer.data(gd_geom_11_off + 174 * ccomps * dcomps);

            auto g_x_y_zzzz_xy = cbuffer.data(gd_geom_11_off + 175 * ccomps * dcomps);

            auto g_x_y_zzzz_xz = cbuffer.data(gd_geom_11_off + 176 * ccomps * dcomps);

            auto g_x_y_zzzz_yy = cbuffer.data(gd_geom_11_off + 177 * ccomps * dcomps);

            auto g_x_y_zzzz_yz = cbuffer.data(gd_geom_11_off + 178 * ccomps * dcomps);

            auto g_x_y_zzzz_zz = cbuffer.data(gd_geom_11_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_zzz_xx, g_x_y_zzz_xxz, g_x_y_zzz_xy, g_x_y_zzz_xyz, g_x_y_zzz_xz, g_x_y_zzz_xzz, g_x_y_zzz_yy, g_x_y_zzz_yyz, g_x_y_zzz_yz, g_x_y_zzz_yzz, g_x_y_zzz_zz, g_x_y_zzz_zzz, g_x_y_zzzz_xx, g_x_y_zzzz_xy, g_x_y_zzzz_xz, g_x_y_zzzz_yy, g_x_y_zzzz_yz, g_x_y_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_zzzz_xx[k] = -g_x_y_zzz_xx[k] * ab_z + g_x_y_zzz_xxz[k];

                g_x_y_zzzz_xy[k] = -g_x_y_zzz_xy[k] * ab_z + g_x_y_zzz_xyz[k];

                g_x_y_zzzz_xz[k] = -g_x_y_zzz_xz[k] * ab_z + g_x_y_zzz_xzz[k];

                g_x_y_zzzz_yy[k] = -g_x_y_zzz_yy[k] * ab_z + g_x_y_zzz_yyz[k];

                g_x_y_zzzz_yz[k] = -g_x_y_zzz_yz[k] * ab_z + g_x_y_zzz_yzz[k];

                g_x_y_zzzz_zz[k] = -g_x_y_zzz_zz[k] * ab_z + g_x_y_zzz_zzz[k];
            }

            /// Set up 180-186 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxxx_xx = cbuffer.data(gd_geom_11_off + 180 * ccomps * dcomps);

            auto g_x_z_xxxx_xy = cbuffer.data(gd_geom_11_off + 181 * ccomps * dcomps);

            auto g_x_z_xxxx_xz = cbuffer.data(gd_geom_11_off + 182 * ccomps * dcomps);

            auto g_x_z_xxxx_yy = cbuffer.data(gd_geom_11_off + 183 * ccomps * dcomps);

            auto g_x_z_xxxx_yz = cbuffer.data(gd_geom_11_off + 184 * ccomps * dcomps);

            auto g_x_z_xxxx_zz = cbuffer.data(gd_geom_11_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxx_xx, g_0_z_xxx_xy, g_0_z_xxx_xz, g_0_z_xxx_yy, g_0_z_xxx_yz, g_0_z_xxx_zz, g_x_z_xxx_xx, g_x_z_xxx_xxx, g_x_z_xxx_xxy, g_x_z_xxx_xxz, g_x_z_xxx_xy, g_x_z_xxx_xyy, g_x_z_xxx_xyz, g_x_z_xxx_xz, g_x_z_xxx_xzz, g_x_z_xxx_yy, g_x_z_xxx_yz, g_x_z_xxx_zz, g_x_z_xxxx_xx, g_x_z_xxxx_xy, g_x_z_xxxx_xz, g_x_z_xxxx_yy, g_x_z_xxxx_yz, g_x_z_xxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxxx_xx[k] = -g_0_z_xxx_xx[k] - g_x_z_xxx_xx[k] * ab_x + g_x_z_xxx_xxx[k];

                g_x_z_xxxx_xy[k] = -g_0_z_xxx_xy[k] - g_x_z_xxx_xy[k] * ab_x + g_x_z_xxx_xxy[k];

                g_x_z_xxxx_xz[k] = -g_0_z_xxx_xz[k] - g_x_z_xxx_xz[k] * ab_x + g_x_z_xxx_xxz[k];

                g_x_z_xxxx_yy[k] = -g_0_z_xxx_yy[k] - g_x_z_xxx_yy[k] * ab_x + g_x_z_xxx_xyy[k];

                g_x_z_xxxx_yz[k] = -g_0_z_xxx_yz[k] - g_x_z_xxx_yz[k] * ab_x + g_x_z_xxx_xyz[k];

                g_x_z_xxxx_zz[k] = -g_0_z_xxx_zz[k] - g_x_z_xxx_zz[k] * ab_x + g_x_z_xxx_xzz[k];
            }

            /// Set up 186-192 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxxy_xx = cbuffer.data(gd_geom_11_off + 186 * ccomps * dcomps);

            auto g_x_z_xxxy_xy = cbuffer.data(gd_geom_11_off + 187 * ccomps * dcomps);

            auto g_x_z_xxxy_xz = cbuffer.data(gd_geom_11_off + 188 * ccomps * dcomps);

            auto g_x_z_xxxy_yy = cbuffer.data(gd_geom_11_off + 189 * ccomps * dcomps);

            auto g_x_z_xxxy_yz = cbuffer.data(gd_geom_11_off + 190 * ccomps * dcomps);

            auto g_x_z_xxxy_zz = cbuffer.data(gd_geom_11_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xxx_xx, g_x_z_xxx_xxy, g_x_z_xxx_xy, g_x_z_xxx_xyy, g_x_z_xxx_xyz, g_x_z_xxx_xz, g_x_z_xxx_yy, g_x_z_xxx_yyy, g_x_z_xxx_yyz, g_x_z_xxx_yz, g_x_z_xxx_yzz, g_x_z_xxx_zz, g_x_z_xxxy_xx, g_x_z_xxxy_xy, g_x_z_xxxy_xz, g_x_z_xxxy_yy, g_x_z_xxxy_yz, g_x_z_xxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxxy_xx[k] = -g_x_z_xxx_xx[k] * ab_y + g_x_z_xxx_xxy[k];

                g_x_z_xxxy_xy[k] = -g_x_z_xxx_xy[k] * ab_y + g_x_z_xxx_xyy[k];

                g_x_z_xxxy_xz[k] = -g_x_z_xxx_xz[k] * ab_y + g_x_z_xxx_xyz[k];

                g_x_z_xxxy_yy[k] = -g_x_z_xxx_yy[k] * ab_y + g_x_z_xxx_yyy[k];

                g_x_z_xxxy_yz[k] = -g_x_z_xxx_yz[k] * ab_y + g_x_z_xxx_yyz[k];

                g_x_z_xxxy_zz[k] = -g_x_z_xxx_zz[k] * ab_y + g_x_z_xxx_yzz[k];
            }

            /// Set up 192-198 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxxz_xx = cbuffer.data(gd_geom_11_off + 192 * ccomps * dcomps);

            auto g_x_z_xxxz_xy = cbuffer.data(gd_geom_11_off + 193 * ccomps * dcomps);

            auto g_x_z_xxxz_xz = cbuffer.data(gd_geom_11_off + 194 * ccomps * dcomps);

            auto g_x_z_xxxz_yy = cbuffer.data(gd_geom_11_off + 195 * ccomps * dcomps);

            auto g_x_z_xxxz_yz = cbuffer.data(gd_geom_11_off + 196 * ccomps * dcomps);

            auto g_x_z_xxxz_zz = cbuffer.data(gd_geom_11_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxz_xx, g_0_z_xxz_xy, g_0_z_xxz_xz, g_0_z_xxz_yy, g_0_z_xxz_yz, g_0_z_xxz_zz, g_x_z_xxxz_xx, g_x_z_xxxz_xy, g_x_z_xxxz_xz, g_x_z_xxxz_yy, g_x_z_xxxz_yz, g_x_z_xxxz_zz, g_x_z_xxz_xx, g_x_z_xxz_xxx, g_x_z_xxz_xxy, g_x_z_xxz_xxz, g_x_z_xxz_xy, g_x_z_xxz_xyy, g_x_z_xxz_xyz, g_x_z_xxz_xz, g_x_z_xxz_xzz, g_x_z_xxz_yy, g_x_z_xxz_yz, g_x_z_xxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxxz_xx[k] = -g_0_z_xxz_xx[k] - g_x_z_xxz_xx[k] * ab_x + g_x_z_xxz_xxx[k];

                g_x_z_xxxz_xy[k] = -g_0_z_xxz_xy[k] - g_x_z_xxz_xy[k] * ab_x + g_x_z_xxz_xxy[k];

                g_x_z_xxxz_xz[k] = -g_0_z_xxz_xz[k] - g_x_z_xxz_xz[k] * ab_x + g_x_z_xxz_xxz[k];

                g_x_z_xxxz_yy[k] = -g_0_z_xxz_yy[k] - g_x_z_xxz_yy[k] * ab_x + g_x_z_xxz_xyy[k];

                g_x_z_xxxz_yz[k] = -g_0_z_xxz_yz[k] - g_x_z_xxz_yz[k] * ab_x + g_x_z_xxz_xyz[k];

                g_x_z_xxxz_zz[k] = -g_0_z_xxz_zz[k] - g_x_z_xxz_zz[k] * ab_x + g_x_z_xxz_xzz[k];
            }

            /// Set up 198-204 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxyy_xx = cbuffer.data(gd_geom_11_off + 198 * ccomps * dcomps);

            auto g_x_z_xxyy_xy = cbuffer.data(gd_geom_11_off + 199 * ccomps * dcomps);

            auto g_x_z_xxyy_xz = cbuffer.data(gd_geom_11_off + 200 * ccomps * dcomps);

            auto g_x_z_xxyy_yy = cbuffer.data(gd_geom_11_off + 201 * ccomps * dcomps);

            auto g_x_z_xxyy_yz = cbuffer.data(gd_geom_11_off + 202 * ccomps * dcomps);

            auto g_x_z_xxyy_zz = cbuffer.data(gd_geom_11_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xxy_xx, g_x_z_xxy_xxy, g_x_z_xxy_xy, g_x_z_xxy_xyy, g_x_z_xxy_xyz, g_x_z_xxy_xz, g_x_z_xxy_yy, g_x_z_xxy_yyy, g_x_z_xxy_yyz, g_x_z_xxy_yz, g_x_z_xxy_yzz, g_x_z_xxy_zz, g_x_z_xxyy_xx, g_x_z_xxyy_xy, g_x_z_xxyy_xz, g_x_z_xxyy_yy, g_x_z_xxyy_yz, g_x_z_xxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxyy_xx[k] = -g_x_z_xxy_xx[k] * ab_y + g_x_z_xxy_xxy[k];

                g_x_z_xxyy_xy[k] = -g_x_z_xxy_xy[k] * ab_y + g_x_z_xxy_xyy[k];

                g_x_z_xxyy_xz[k] = -g_x_z_xxy_xz[k] * ab_y + g_x_z_xxy_xyz[k];

                g_x_z_xxyy_yy[k] = -g_x_z_xxy_yy[k] * ab_y + g_x_z_xxy_yyy[k];

                g_x_z_xxyy_yz[k] = -g_x_z_xxy_yz[k] * ab_y + g_x_z_xxy_yyz[k];

                g_x_z_xxyy_zz[k] = -g_x_z_xxy_zz[k] * ab_y + g_x_z_xxy_yzz[k];
            }

            /// Set up 204-210 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxyz_xx = cbuffer.data(gd_geom_11_off + 204 * ccomps * dcomps);

            auto g_x_z_xxyz_xy = cbuffer.data(gd_geom_11_off + 205 * ccomps * dcomps);

            auto g_x_z_xxyz_xz = cbuffer.data(gd_geom_11_off + 206 * ccomps * dcomps);

            auto g_x_z_xxyz_yy = cbuffer.data(gd_geom_11_off + 207 * ccomps * dcomps);

            auto g_x_z_xxyz_yz = cbuffer.data(gd_geom_11_off + 208 * ccomps * dcomps);

            auto g_x_z_xxyz_zz = cbuffer.data(gd_geom_11_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xxyz_xx, g_x_z_xxyz_xy, g_x_z_xxyz_xz, g_x_z_xxyz_yy, g_x_z_xxyz_yz, g_x_z_xxyz_zz, g_x_z_xxz_xx, g_x_z_xxz_xxy, g_x_z_xxz_xy, g_x_z_xxz_xyy, g_x_z_xxz_xyz, g_x_z_xxz_xz, g_x_z_xxz_yy, g_x_z_xxz_yyy, g_x_z_xxz_yyz, g_x_z_xxz_yz, g_x_z_xxz_yzz, g_x_z_xxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxyz_xx[k] = -g_x_z_xxz_xx[k] * ab_y + g_x_z_xxz_xxy[k];

                g_x_z_xxyz_xy[k] = -g_x_z_xxz_xy[k] * ab_y + g_x_z_xxz_xyy[k];

                g_x_z_xxyz_xz[k] = -g_x_z_xxz_xz[k] * ab_y + g_x_z_xxz_xyz[k];

                g_x_z_xxyz_yy[k] = -g_x_z_xxz_yy[k] * ab_y + g_x_z_xxz_yyy[k];

                g_x_z_xxyz_yz[k] = -g_x_z_xxz_yz[k] * ab_y + g_x_z_xxz_yyz[k];

                g_x_z_xxyz_zz[k] = -g_x_z_xxz_zz[k] * ab_y + g_x_z_xxz_yzz[k];
            }

            /// Set up 210-216 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxzz_xx = cbuffer.data(gd_geom_11_off + 210 * ccomps * dcomps);

            auto g_x_z_xxzz_xy = cbuffer.data(gd_geom_11_off + 211 * ccomps * dcomps);

            auto g_x_z_xxzz_xz = cbuffer.data(gd_geom_11_off + 212 * ccomps * dcomps);

            auto g_x_z_xxzz_yy = cbuffer.data(gd_geom_11_off + 213 * ccomps * dcomps);

            auto g_x_z_xxzz_yz = cbuffer.data(gd_geom_11_off + 214 * ccomps * dcomps);

            auto g_x_z_xxzz_zz = cbuffer.data(gd_geom_11_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzz_xx, g_0_z_xzz_xy, g_0_z_xzz_xz, g_0_z_xzz_yy, g_0_z_xzz_yz, g_0_z_xzz_zz, g_x_z_xxzz_xx, g_x_z_xxzz_xy, g_x_z_xxzz_xz, g_x_z_xxzz_yy, g_x_z_xxzz_yz, g_x_z_xxzz_zz, g_x_z_xzz_xx, g_x_z_xzz_xxx, g_x_z_xzz_xxy, g_x_z_xzz_xxz, g_x_z_xzz_xy, g_x_z_xzz_xyy, g_x_z_xzz_xyz, g_x_z_xzz_xz, g_x_z_xzz_xzz, g_x_z_xzz_yy, g_x_z_xzz_yz, g_x_z_xzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxzz_xx[k] = -g_0_z_xzz_xx[k] - g_x_z_xzz_xx[k] * ab_x + g_x_z_xzz_xxx[k];

                g_x_z_xxzz_xy[k] = -g_0_z_xzz_xy[k] - g_x_z_xzz_xy[k] * ab_x + g_x_z_xzz_xxy[k];

                g_x_z_xxzz_xz[k] = -g_0_z_xzz_xz[k] - g_x_z_xzz_xz[k] * ab_x + g_x_z_xzz_xxz[k];

                g_x_z_xxzz_yy[k] = -g_0_z_xzz_yy[k] - g_x_z_xzz_yy[k] * ab_x + g_x_z_xzz_xyy[k];

                g_x_z_xxzz_yz[k] = -g_0_z_xzz_yz[k] - g_x_z_xzz_yz[k] * ab_x + g_x_z_xzz_xyz[k];

                g_x_z_xxzz_zz[k] = -g_0_z_xzz_zz[k] - g_x_z_xzz_zz[k] * ab_x + g_x_z_xzz_xzz[k];
            }

            /// Set up 216-222 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyyy_xx = cbuffer.data(gd_geom_11_off + 216 * ccomps * dcomps);

            auto g_x_z_xyyy_xy = cbuffer.data(gd_geom_11_off + 217 * ccomps * dcomps);

            auto g_x_z_xyyy_xz = cbuffer.data(gd_geom_11_off + 218 * ccomps * dcomps);

            auto g_x_z_xyyy_yy = cbuffer.data(gd_geom_11_off + 219 * ccomps * dcomps);

            auto g_x_z_xyyy_yz = cbuffer.data(gd_geom_11_off + 220 * ccomps * dcomps);

            auto g_x_z_xyyy_zz = cbuffer.data(gd_geom_11_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xyy_xx, g_x_z_xyy_xxy, g_x_z_xyy_xy, g_x_z_xyy_xyy, g_x_z_xyy_xyz, g_x_z_xyy_xz, g_x_z_xyy_yy, g_x_z_xyy_yyy, g_x_z_xyy_yyz, g_x_z_xyy_yz, g_x_z_xyy_yzz, g_x_z_xyy_zz, g_x_z_xyyy_xx, g_x_z_xyyy_xy, g_x_z_xyyy_xz, g_x_z_xyyy_yy, g_x_z_xyyy_yz, g_x_z_xyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyyy_xx[k] = -g_x_z_xyy_xx[k] * ab_y + g_x_z_xyy_xxy[k];

                g_x_z_xyyy_xy[k] = -g_x_z_xyy_xy[k] * ab_y + g_x_z_xyy_xyy[k];

                g_x_z_xyyy_xz[k] = -g_x_z_xyy_xz[k] * ab_y + g_x_z_xyy_xyz[k];

                g_x_z_xyyy_yy[k] = -g_x_z_xyy_yy[k] * ab_y + g_x_z_xyy_yyy[k];

                g_x_z_xyyy_yz[k] = -g_x_z_xyy_yz[k] * ab_y + g_x_z_xyy_yyz[k];

                g_x_z_xyyy_zz[k] = -g_x_z_xyy_zz[k] * ab_y + g_x_z_xyy_yzz[k];
            }

            /// Set up 222-228 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyyz_xx = cbuffer.data(gd_geom_11_off + 222 * ccomps * dcomps);

            auto g_x_z_xyyz_xy = cbuffer.data(gd_geom_11_off + 223 * ccomps * dcomps);

            auto g_x_z_xyyz_xz = cbuffer.data(gd_geom_11_off + 224 * ccomps * dcomps);

            auto g_x_z_xyyz_yy = cbuffer.data(gd_geom_11_off + 225 * ccomps * dcomps);

            auto g_x_z_xyyz_yz = cbuffer.data(gd_geom_11_off + 226 * ccomps * dcomps);

            auto g_x_z_xyyz_zz = cbuffer.data(gd_geom_11_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xyyz_xx, g_x_z_xyyz_xy, g_x_z_xyyz_xz, g_x_z_xyyz_yy, g_x_z_xyyz_yz, g_x_z_xyyz_zz, g_x_z_xyz_xx, g_x_z_xyz_xxy, g_x_z_xyz_xy, g_x_z_xyz_xyy, g_x_z_xyz_xyz, g_x_z_xyz_xz, g_x_z_xyz_yy, g_x_z_xyz_yyy, g_x_z_xyz_yyz, g_x_z_xyz_yz, g_x_z_xyz_yzz, g_x_z_xyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyyz_xx[k] = -g_x_z_xyz_xx[k] * ab_y + g_x_z_xyz_xxy[k];

                g_x_z_xyyz_xy[k] = -g_x_z_xyz_xy[k] * ab_y + g_x_z_xyz_xyy[k];

                g_x_z_xyyz_xz[k] = -g_x_z_xyz_xz[k] * ab_y + g_x_z_xyz_xyz[k];

                g_x_z_xyyz_yy[k] = -g_x_z_xyz_yy[k] * ab_y + g_x_z_xyz_yyy[k];

                g_x_z_xyyz_yz[k] = -g_x_z_xyz_yz[k] * ab_y + g_x_z_xyz_yyz[k];

                g_x_z_xyyz_zz[k] = -g_x_z_xyz_zz[k] * ab_y + g_x_z_xyz_yzz[k];
            }

            /// Set up 228-234 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyzz_xx = cbuffer.data(gd_geom_11_off + 228 * ccomps * dcomps);

            auto g_x_z_xyzz_xy = cbuffer.data(gd_geom_11_off + 229 * ccomps * dcomps);

            auto g_x_z_xyzz_xz = cbuffer.data(gd_geom_11_off + 230 * ccomps * dcomps);

            auto g_x_z_xyzz_yy = cbuffer.data(gd_geom_11_off + 231 * ccomps * dcomps);

            auto g_x_z_xyzz_yz = cbuffer.data(gd_geom_11_off + 232 * ccomps * dcomps);

            auto g_x_z_xyzz_zz = cbuffer.data(gd_geom_11_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xyzz_xx, g_x_z_xyzz_xy, g_x_z_xyzz_xz, g_x_z_xyzz_yy, g_x_z_xyzz_yz, g_x_z_xyzz_zz, g_x_z_xzz_xx, g_x_z_xzz_xxy, g_x_z_xzz_xy, g_x_z_xzz_xyy, g_x_z_xzz_xyz, g_x_z_xzz_xz, g_x_z_xzz_yy, g_x_z_xzz_yyy, g_x_z_xzz_yyz, g_x_z_xzz_yz, g_x_z_xzz_yzz, g_x_z_xzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyzz_xx[k] = -g_x_z_xzz_xx[k] * ab_y + g_x_z_xzz_xxy[k];

                g_x_z_xyzz_xy[k] = -g_x_z_xzz_xy[k] * ab_y + g_x_z_xzz_xyy[k];

                g_x_z_xyzz_xz[k] = -g_x_z_xzz_xz[k] * ab_y + g_x_z_xzz_xyz[k];

                g_x_z_xyzz_yy[k] = -g_x_z_xzz_yy[k] * ab_y + g_x_z_xzz_yyy[k];

                g_x_z_xyzz_yz[k] = -g_x_z_xzz_yz[k] * ab_y + g_x_z_xzz_yyz[k];

                g_x_z_xyzz_zz[k] = -g_x_z_xzz_zz[k] * ab_y + g_x_z_xzz_yzz[k];
            }

            /// Set up 234-240 components of targeted buffer : cbuffer.data(

            auto g_x_z_xzzz_xx = cbuffer.data(gd_geom_11_off + 234 * ccomps * dcomps);

            auto g_x_z_xzzz_xy = cbuffer.data(gd_geom_11_off + 235 * ccomps * dcomps);

            auto g_x_z_xzzz_xz = cbuffer.data(gd_geom_11_off + 236 * ccomps * dcomps);

            auto g_x_z_xzzz_yy = cbuffer.data(gd_geom_11_off + 237 * ccomps * dcomps);

            auto g_x_z_xzzz_yz = cbuffer.data(gd_geom_11_off + 238 * ccomps * dcomps);

            auto g_x_z_xzzz_zz = cbuffer.data(gd_geom_11_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzz_xx, g_0_z_zzz_xy, g_0_z_zzz_xz, g_0_z_zzz_yy, g_0_z_zzz_yz, g_0_z_zzz_zz, g_x_z_xzzz_xx, g_x_z_xzzz_xy, g_x_z_xzzz_xz, g_x_z_xzzz_yy, g_x_z_xzzz_yz, g_x_z_xzzz_zz, g_x_z_zzz_xx, g_x_z_zzz_xxx, g_x_z_zzz_xxy, g_x_z_zzz_xxz, g_x_z_zzz_xy, g_x_z_zzz_xyy, g_x_z_zzz_xyz, g_x_z_zzz_xz, g_x_z_zzz_xzz, g_x_z_zzz_yy, g_x_z_zzz_yz, g_x_z_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xzzz_xx[k] = -g_0_z_zzz_xx[k] - g_x_z_zzz_xx[k] * ab_x + g_x_z_zzz_xxx[k];

                g_x_z_xzzz_xy[k] = -g_0_z_zzz_xy[k] - g_x_z_zzz_xy[k] * ab_x + g_x_z_zzz_xxy[k];

                g_x_z_xzzz_xz[k] = -g_0_z_zzz_xz[k] - g_x_z_zzz_xz[k] * ab_x + g_x_z_zzz_xxz[k];

                g_x_z_xzzz_yy[k] = -g_0_z_zzz_yy[k] - g_x_z_zzz_yy[k] * ab_x + g_x_z_zzz_xyy[k];

                g_x_z_xzzz_yz[k] = -g_0_z_zzz_yz[k] - g_x_z_zzz_yz[k] * ab_x + g_x_z_zzz_xyz[k];

                g_x_z_xzzz_zz[k] = -g_0_z_zzz_zz[k] - g_x_z_zzz_zz[k] * ab_x + g_x_z_zzz_xzz[k];
            }

            /// Set up 240-246 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyyy_xx = cbuffer.data(gd_geom_11_off + 240 * ccomps * dcomps);

            auto g_x_z_yyyy_xy = cbuffer.data(gd_geom_11_off + 241 * ccomps * dcomps);

            auto g_x_z_yyyy_xz = cbuffer.data(gd_geom_11_off + 242 * ccomps * dcomps);

            auto g_x_z_yyyy_yy = cbuffer.data(gd_geom_11_off + 243 * ccomps * dcomps);

            auto g_x_z_yyyy_yz = cbuffer.data(gd_geom_11_off + 244 * ccomps * dcomps);

            auto g_x_z_yyyy_zz = cbuffer.data(gd_geom_11_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yyy_xx, g_x_z_yyy_xxy, g_x_z_yyy_xy, g_x_z_yyy_xyy, g_x_z_yyy_xyz, g_x_z_yyy_xz, g_x_z_yyy_yy, g_x_z_yyy_yyy, g_x_z_yyy_yyz, g_x_z_yyy_yz, g_x_z_yyy_yzz, g_x_z_yyy_zz, g_x_z_yyyy_xx, g_x_z_yyyy_xy, g_x_z_yyyy_xz, g_x_z_yyyy_yy, g_x_z_yyyy_yz, g_x_z_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyyy_xx[k] = -g_x_z_yyy_xx[k] * ab_y + g_x_z_yyy_xxy[k];

                g_x_z_yyyy_xy[k] = -g_x_z_yyy_xy[k] * ab_y + g_x_z_yyy_xyy[k];

                g_x_z_yyyy_xz[k] = -g_x_z_yyy_xz[k] * ab_y + g_x_z_yyy_xyz[k];

                g_x_z_yyyy_yy[k] = -g_x_z_yyy_yy[k] * ab_y + g_x_z_yyy_yyy[k];

                g_x_z_yyyy_yz[k] = -g_x_z_yyy_yz[k] * ab_y + g_x_z_yyy_yyz[k];

                g_x_z_yyyy_zz[k] = -g_x_z_yyy_zz[k] * ab_y + g_x_z_yyy_yzz[k];
            }

            /// Set up 246-252 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyyz_xx = cbuffer.data(gd_geom_11_off + 246 * ccomps * dcomps);

            auto g_x_z_yyyz_xy = cbuffer.data(gd_geom_11_off + 247 * ccomps * dcomps);

            auto g_x_z_yyyz_xz = cbuffer.data(gd_geom_11_off + 248 * ccomps * dcomps);

            auto g_x_z_yyyz_yy = cbuffer.data(gd_geom_11_off + 249 * ccomps * dcomps);

            auto g_x_z_yyyz_yz = cbuffer.data(gd_geom_11_off + 250 * ccomps * dcomps);

            auto g_x_z_yyyz_zz = cbuffer.data(gd_geom_11_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yyyz_xx, g_x_z_yyyz_xy, g_x_z_yyyz_xz, g_x_z_yyyz_yy, g_x_z_yyyz_yz, g_x_z_yyyz_zz, g_x_z_yyz_xx, g_x_z_yyz_xxy, g_x_z_yyz_xy, g_x_z_yyz_xyy, g_x_z_yyz_xyz, g_x_z_yyz_xz, g_x_z_yyz_yy, g_x_z_yyz_yyy, g_x_z_yyz_yyz, g_x_z_yyz_yz, g_x_z_yyz_yzz, g_x_z_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyyz_xx[k] = -g_x_z_yyz_xx[k] * ab_y + g_x_z_yyz_xxy[k];

                g_x_z_yyyz_xy[k] = -g_x_z_yyz_xy[k] * ab_y + g_x_z_yyz_xyy[k];

                g_x_z_yyyz_xz[k] = -g_x_z_yyz_xz[k] * ab_y + g_x_z_yyz_xyz[k];

                g_x_z_yyyz_yy[k] = -g_x_z_yyz_yy[k] * ab_y + g_x_z_yyz_yyy[k];

                g_x_z_yyyz_yz[k] = -g_x_z_yyz_yz[k] * ab_y + g_x_z_yyz_yyz[k];

                g_x_z_yyyz_zz[k] = -g_x_z_yyz_zz[k] * ab_y + g_x_z_yyz_yzz[k];
            }

            /// Set up 252-258 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyzz_xx = cbuffer.data(gd_geom_11_off + 252 * ccomps * dcomps);

            auto g_x_z_yyzz_xy = cbuffer.data(gd_geom_11_off + 253 * ccomps * dcomps);

            auto g_x_z_yyzz_xz = cbuffer.data(gd_geom_11_off + 254 * ccomps * dcomps);

            auto g_x_z_yyzz_yy = cbuffer.data(gd_geom_11_off + 255 * ccomps * dcomps);

            auto g_x_z_yyzz_yz = cbuffer.data(gd_geom_11_off + 256 * ccomps * dcomps);

            auto g_x_z_yyzz_zz = cbuffer.data(gd_geom_11_off + 257 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yyzz_xx, g_x_z_yyzz_xy, g_x_z_yyzz_xz, g_x_z_yyzz_yy, g_x_z_yyzz_yz, g_x_z_yyzz_zz, g_x_z_yzz_xx, g_x_z_yzz_xxy, g_x_z_yzz_xy, g_x_z_yzz_xyy, g_x_z_yzz_xyz, g_x_z_yzz_xz, g_x_z_yzz_yy, g_x_z_yzz_yyy, g_x_z_yzz_yyz, g_x_z_yzz_yz, g_x_z_yzz_yzz, g_x_z_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyzz_xx[k] = -g_x_z_yzz_xx[k] * ab_y + g_x_z_yzz_xxy[k];

                g_x_z_yyzz_xy[k] = -g_x_z_yzz_xy[k] * ab_y + g_x_z_yzz_xyy[k];

                g_x_z_yyzz_xz[k] = -g_x_z_yzz_xz[k] * ab_y + g_x_z_yzz_xyz[k];

                g_x_z_yyzz_yy[k] = -g_x_z_yzz_yy[k] * ab_y + g_x_z_yzz_yyy[k];

                g_x_z_yyzz_yz[k] = -g_x_z_yzz_yz[k] * ab_y + g_x_z_yzz_yyz[k];

                g_x_z_yyzz_zz[k] = -g_x_z_yzz_zz[k] * ab_y + g_x_z_yzz_yzz[k];
            }

            /// Set up 258-264 components of targeted buffer : cbuffer.data(

            auto g_x_z_yzzz_xx = cbuffer.data(gd_geom_11_off + 258 * ccomps * dcomps);

            auto g_x_z_yzzz_xy = cbuffer.data(gd_geom_11_off + 259 * ccomps * dcomps);

            auto g_x_z_yzzz_xz = cbuffer.data(gd_geom_11_off + 260 * ccomps * dcomps);

            auto g_x_z_yzzz_yy = cbuffer.data(gd_geom_11_off + 261 * ccomps * dcomps);

            auto g_x_z_yzzz_yz = cbuffer.data(gd_geom_11_off + 262 * ccomps * dcomps);

            auto g_x_z_yzzz_zz = cbuffer.data(gd_geom_11_off + 263 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yzzz_xx, g_x_z_yzzz_xy, g_x_z_yzzz_xz, g_x_z_yzzz_yy, g_x_z_yzzz_yz, g_x_z_yzzz_zz, g_x_z_zzz_xx, g_x_z_zzz_xxy, g_x_z_zzz_xy, g_x_z_zzz_xyy, g_x_z_zzz_xyz, g_x_z_zzz_xz, g_x_z_zzz_yy, g_x_z_zzz_yyy, g_x_z_zzz_yyz, g_x_z_zzz_yz, g_x_z_zzz_yzz, g_x_z_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yzzz_xx[k] = -g_x_z_zzz_xx[k] * ab_y + g_x_z_zzz_xxy[k];

                g_x_z_yzzz_xy[k] = -g_x_z_zzz_xy[k] * ab_y + g_x_z_zzz_xyy[k];

                g_x_z_yzzz_xz[k] = -g_x_z_zzz_xz[k] * ab_y + g_x_z_zzz_xyz[k];

                g_x_z_yzzz_yy[k] = -g_x_z_zzz_yy[k] * ab_y + g_x_z_zzz_yyy[k];

                g_x_z_yzzz_yz[k] = -g_x_z_zzz_yz[k] * ab_y + g_x_z_zzz_yyz[k];

                g_x_z_yzzz_zz[k] = -g_x_z_zzz_zz[k] * ab_y + g_x_z_zzz_yzz[k];
            }

            /// Set up 264-270 components of targeted buffer : cbuffer.data(

            auto g_x_z_zzzz_xx = cbuffer.data(gd_geom_11_off + 264 * ccomps * dcomps);

            auto g_x_z_zzzz_xy = cbuffer.data(gd_geom_11_off + 265 * ccomps * dcomps);

            auto g_x_z_zzzz_xz = cbuffer.data(gd_geom_11_off + 266 * ccomps * dcomps);

            auto g_x_z_zzzz_yy = cbuffer.data(gd_geom_11_off + 267 * ccomps * dcomps);

            auto g_x_z_zzzz_yz = cbuffer.data(gd_geom_11_off + 268 * ccomps * dcomps);

            auto g_x_z_zzzz_zz = cbuffer.data(gd_geom_11_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzz_xx, g_x_0_zzz_xy, g_x_0_zzz_xz, g_x_0_zzz_yy, g_x_0_zzz_yz, g_x_0_zzz_zz, g_x_z_zzz_xx, g_x_z_zzz_xxz, g_x_z_zzz_xy, g_x_z_zzz_xyz, g_x_z_zzz_xz, g_x_z_zzz_xzz, g_x_z_zzz_yy, g_x_z_zzz_yyz, g_x_z_zzz_yz, g_x_z_zzz_yzz, g_x_z_zzz_zz, g_x_z_zzz_zzz, g_x_z_zzzz_xx, g_x_z_zzzz_xy, g_x_z_zzzz_xz, g_x_z_zzzz_yy, g_x_z_zzzz_yz, g_x_z_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_zzzz_xx[k] = g_x_0_zzz_xx[k] - g_x_z_zzz_xx[k] * ab_z + g_x_z_zzz_xxz[k];

                g_x_z_zzzz_xy[k] = g_x_0_zzz_xy[k] - g_x_z_zzz_xy[k] * ab_z + g_x_z_zzz_xyz[k];

                g_x_z_zzzz_xz[k] = g_x_0_zzz_xz[k] - g_x_z_zzz_xz[k] * ab_z + g_x_z_zzz_xzz[k];

                g_x_z_zzzz_yy[k] = g_x_0_zzz_yy[k] - g_x_z_zzz_yy[k] * ab_z + g_x_z_zzz_yyz[k];

                g_x_z_zzzz_yz[k] = g_x_0_zzz_yz[k] - g_x_z_zzz_yz[k] * ab_z + g_x_z_zzz_yzz[k];

                g_x_z_zzzz_zz[k] = g_x_0_zzz_zz[k] - g_x_z_zzz_zz[k] * ab_z + g_x_z_zzz_zzz[k];
            }

            /// Set up 270-276 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxxx_xx = cbuffer.data(gd_geom_11_off + 270 * ccomps * dcomps);

            auto g_y_x_xxxx_xy = cbuffer.data(gd_geom_11_off + 271 * ccomps * dcomps);

            auto g_y_x_xxxx_xz = cbuffer.data(gd_geom_11_off + 272 * ccomps * dcomps);

            auto g_y_x_xxxx_yy = cbuffer.data(gd_geom_11_off + 273 * ccomps * dcomps);

            auto g_y_x_xxxx_yz = cbuffer.data(gd_geom_11_off + 274 * ccomps * dcomps);

            auto g_y_x_xxxx_zz = cbuffer.data(gd_geom_11_off + 275 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxx_xx, g_y_0_xxx_xy, g_y_0_xxx_xz, g_y_0_xxx_yy, g_y_0_xxx_yz, g_y_0_xxx_zz, g_y_x_xxx_xx, g_y_x_xxx_xxx, g_y_x_xxx_xxy, g_y_x_xxx_xxz, g_y_x_xxx_xy, g_y_x_xxx_xyy, g_y_x_xxx_xyz, g_y_x_xxx_xz, g_y_x_xxx_xzz, g_y_x_xxx_yy, g_y_x_xxx_yz, g_y_x_xxx_zz, g_y_x_xxxx_xx, g_y_x_xxxx_xy, g_y_x_xxxx_xz, g_y_x_xxxx_yy, g_y_x_xxxx_yz, g_y_x_xxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxxx_xx[k] = g_y_0_xxx_xx[k] - g_y_x_xxx_xx[k] * ab_x + g_y_x_xxx_xxx[k];

                g_y_x_xxxx_xy[k] = g_y_0_xxx_xy[k] - g_y_x_xxx_xy[k] * ab_x + g_y_x_xxx_xxy[k];

                g_y_x_xxxx_xz[k] = g_y_0_xxx_xz[k] - g_y_x_xxx_xz[k] * ab_x + g_y_x_xxx_xxz[k];

                g_y_x_xxxx_yy[k] = g_y_0_xxx_yy[k] - g_y_x_xxx_yy[k] * ab_x + g_y_x_xxx_xyy[k];

                g_y_x_xxxx_yz[k] = g_y_0_xxx_yz[k] - g_y_x_xxx_yz[k] * ab_x + g_y_x_xxx_xyz[k];

                g_y_x_xxxx_zz[k] = g_y_0_xxx_zz[k] - g_y_x_xxx_zz[k] * ab_x + g_y_x_xxx_xzz[k];
            }

            /// Set up 276-282 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxxy_xx = cbuffer.data(gd_geom_11_off + 276 * ccomps * dcomps);

            auto g_y_x_xxxy_xy = cbuffer.data(gd_geom_11_off + 277 * ccomps * dcomps);

            auto g_y_x_xxxy_xz = cbuffer.data(gd_geom_11_off + 278 * ccomps * dcomps);

            auto g_y_x_xxxy_yy = cbuffer.data(gd_geom_11_off + 279 * ccomps * dcomps);

            auto g_y_x_xxxy_yz = cbuffer.data(gd_geom_11_off + 280 * ccomps * dcomps);

            auto g_y_x_xxxy_zz = cbuffer.data(gd_geom_11_off + 281 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxy_xx, g_y_0_xxy_xy, g_y_0_xxy_xz, g_y_0_xxy_yy, g_y_0_xxy_yz, g_y_0_xxy_zz, g_y_x_xxxy_xx, g_y_x_xxxy_xy, g_y_x_xxxy_xz, g_y_x_xxxy_yy, g_y_x_xxxy_yz, g_y_x_xxxy_zz, g_y_x_xxy_xx, g_y_x_xxy_xxx, g_y_x_xxy_xxy, g_y_x_xxy_xxz, g_y_x_xxy_xy, g_y_x_xxy_xyy, g_y_x_xxy_xyz, g_y_x_xxy_xz, g_y_x_xxy_xzz, g_y_x_xxy_yy, g_y_x_xxy_yz, g_y_x_xxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxxy_xx[k] = g_y_0_xxy_xx[k] - g_y_x_xxy_xx[k] * ab_x + g_y_x_xxy_xxx[k];

                g_y_x_xxxy_xy[k] = g_y_0_xxy_xy[k] - g_y_x_xxy_xy[k] * ab_x + g_y_x_xxy_xxy[k];

                g_y_x_xxxy_xz[k] = g_y_0_xxy_xz[k] - g_y_x_xxy_xz[k] * ab_x + g_y_x_xxy_xxz[k];

                g_y_x_xxxy_yy[k] = g_y_0_xxy_yy[k] - g_y_x_xxy_yy[k] * ab_x + g_y_x_xxy_xyy[k];

                g_y_x_xxxy_yz[k] = g_y_0_xxy_yz[k] - g_y_x_xxy_yz[k] * ab_x + g_y_x_xxy_xyz[k];

                g_y_x_xxxy_zz[k] = g_y_0_xxy_zz[k] - g_y_x_xxy_zz[k] * ab_x + g_y_x_xxy_xzz[k];
            }

            /// Set up 282-288 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxxz_xx = cbuffer.data(gd_geom_11_off + 282 * ccomps * dcomps);

            auto g_y_x_xxxz_xy = cbuffer.data(gd_geom_11_off + 283 * ccomps * dcomps);

            auto g_y_x_xxxz_xz = cbuffer.data(gd_geom_11_off + 284 * ccomps * dcomps);

            auto g_y_x_xxxz_yy = cbuffer.data(gd_geom_11_off + 285 * ccomps * dcomps);

            auto g_y_x_xxxz_yz = cbuffer.data(gd_geom_11_off + 286 * ccomps * dcomps);

            auto g_y_x_xxxz_zz = cbuffer.data(gd_geom_11_off + 287 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xxx_xx, g_y_x_xxx_xxz, g_y_x_xxx_xy, g_y_x_xxx_xyz, g_y_x_xxx_xz, g_y_x_xxx_xzz, g_y_x_xxx_yy, g_y_x_xxx_yyz, g_y_x_xxx_yz, g_y_x_xxx_yzz, g_y_x_xxx_zz, g_y_x_xxx_zzz, g_y_x_xxxz_xx, g_y_x_xxxz_xy, g_y_x_xxxz_xz, g_y_x_xxxz_yy, g_y_x_xxxz_yz, g_y_x_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxxz_xx[k] = -g_y_x_xxx_xx[k] * ab_z + g_y_x_xxx_xxz[k];

                g_y_x_xxxz_xy[k] = -g_y_x_xxx_xy[k] * ab_z + g_y_x_xxx_xyz[k];

                g_y_x_xxxz_xz[k] = -g_y_x_xxx_xz[k] * ab_z + g_y_x_xxx_xzz[k];

                g_y_x_xxxz_yy[k] = -g_y_x_xxx_yy[k] * ab_z + g_y_x_xxx_yyz[k];

                g_y_x_xxxz_yz[k] = -g_y_x_xxx_yz[k] * ab_z + g_y_x_xxx_yzz[k];

                g_y_x_xxxz_zz[k] = -g_y_x_xxx_zz[k] * ab_z + g_y_x_xxx_zzz[k];
            }

            /// Set up 288-294 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxyy_xx = cbuffer.data(gd_geom_11_off + 288 * ccomps * dcomps);

            auto g_y_x_xxyy_xy = cbuffer.data(gd_geom_11_off + 289 * ccomps * dcomps);

            auto g_y_x_xxyy_xz = cbuffer.data(gd_geom_11_off + 290 * ccomps * dcomps);

            auto g_y_x_xxyy_yy = cbuffer.data(gd_geom_11_off + 291 * ccomps * dcomps);

            auto g_y_x_xxyy_yz = cbuffer.data(gd_geom_11_off + 292 * ccomps * dcomps);

            auto g_y_x_xxyy_zz = cbuffer.data(gd_geom_11_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyy_xx, g_y_0_xyy_xy, g_y_0_xyy_xz, g_y_0_xyy_yy, g_y_0_xyy_yz, g_y_0_xyy_zz, g_y_x_xxyy_xx, g_y_x_xxyy_xy, g_y_x_xxyy_xz, g_y_x_xxyy_yy, g_y_x_xxyy_yz, g_y_x_xxyy_zz, g_y_x_xyy_xx, g_y_x_xyy_xxx, g_y_x_xyy_xxy, g_y_x_xyy_xxz, g_y_x_xyy_xy, g_y_x_xyy_xyy, g_y_x_xyy_xyz, g_y_x_xyy_xz, g_y_x_xyy_xzz, g_y_x_xyy_yy, g_y_x_xyy_yz, g_y_x_xyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxyy_xx[k] = g_y_0_xyy_xx[k] - g_y_x_xyy_xx[k] * ab_x + g_y_x_xyy_xxx[k];

                g_y_x_xxyy_xy[k] = g_y_0_xyy_xy[k] - g_y_x_xyy_xy[k] * ab_x + g_y_x_xyy_xxy[k];

                g_y_x_xxyy_xz[k] = g_y_0_xyy_xz[k] - g_y_x_xyy_xz[k] * ab_x + g_y_x_xyy_xxz[k];

                g_y_x_xxyy_yy[k] = g_y_0_xyy_yy[k] - g_y_x_xyy_yy[k] * ab_x + g_y_x_xyy_xyy[k];

                g_y_x_xxyy_yz[k] = g_y_0_xyy_yz[k] - g_y_x_xyy_yz[k] * ab_x + g_y_x_xyy_xyz[k];

                g_y_x_xxyy_zz[k] = g_y_0_xyy_zz[k] - g_y_x_xyy_zz[k] * ab_x + g_y_x_xyy_xzz[k];
            }

            /// Set up 294-300 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxyz_xx = cbuffer.data(gd_geom_11_off + 294 * ccomps * dcomps);

            auto g_y_x_xxyz_xy = cbuffer.data(gd_geom_11_off + 295 * ccomps * dcomps);

            auto g_y_x_xxyz_xz = cbuffer.data(gd_geom_11_off + 296 * ccomps * dcomps);

            auto g_y_x_xxyz_yy = cbuffer.data(gd_geom_11_off + 297 * ccomps * dcomps);

            auto g_y_x_xxyz_yz = cbuffer.data(gd_geom_11_off + 298 * ccomps * dcomps);

            auto g_y_x_xxyz_zz = cbuffer.data(gd_geom_11_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xxy_xx, g_y_x_xxy_xxz, g_y_x_xxy_xy, g_y_x_xxy_xyz, g_y_x_xxy_xz, g_y_x_xxy_xzz, g_y_x_xxy_yy, g_y_x_xxy_yyz, g_y_x_xxy_yz, g_y_x_xxy_yzz, g_y_x_xxy_zz, g_y_x_xxy_zzz, g_y_x_xxyz_xx, g_y_x_xxyz_xy, g_y_x_xxyz_xz, g_y_x_xxyz_yy, g_y_x_xxyz_yz, g_y_x_xxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxyz_xx[k] = -g_y_x_xxy_xx[k] * ab_z + g_y_x_xxy_xxz[k];

                g_y_x_xxyz_xy[k] = -g_y_x_xxy_xy[k] * ab_z + g_y_x_xxy_xyz[k];

                g_y_x_xxyz_xz[k] = -g_y_x_xxy_xz[k] * ab_z + g_y_x_xxy_xzz[k];

                g_y_x_xxyz_yy[k] = -g_y_x_xxy_yy[k] * ab_z + g_y_x_xxy_yyz[k];

                g_y_x_xxyz_yz[k] = -g_y_x_xxy_yz[k] * ab_z + g_y_x_xxy_yzz[k];

                g_y_x_xxyz_zz[k] = -g_y_x_xxy_zz[k] * ab_z + g_y_x_xxy_zzz[k];
            }

            /// Set up 300-306 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxzz_xx = cbuffer.data(gd_geom_11_off + 300 * ccomps * dcomps);

            auto g_y_x_xxzz_xy = cbuffer.data(gd_geom_11_off + 301 * ccomps * dcomps);

            auto g_y_x_xxzz_xz = cbuffer.data(gd_geom_11_off + 302 * ccomps * dcomps);

            auto g_y_x_xxzz_yy = cbuffer.data(gd_geom_11_off + 303 * ccomps * dcomps);

            auto g_y_x_xxzz_yz = cbuffer.data(gd_geom_11_off + 304 * ccomps * dcomps);

            auto g_y_x_xxzz_zz = cbuffer.data(gd_geom_11_off + 305 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xxz_xx, g_y_x_xxz_xxz, g_y_x_xxz_xy, g_y_x_xxz_xyz, g_y_x_xxz_xz, g_y_x_xxz_xzz, g_y_x_xxz_yy, g_y_x_xxz_yyz, g_y_x_xxz_yz, g_y_x_xxz_yzz, g_y_x_xxz_zz, g_y_x_xxz_zzz, g_y_x_xxzz_xx, g_y_x_xxzz_xy, g_y_x_xxzz_xz, g_y_x_xxzz_yy, g_y_x_xxzz_yz, g_y_x_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxzz_xx[k] = -g_y_x_xxz_xx[k] * ab_z + g_y_x_xxz_xxz[k];

                g_y_x_xxzz_xy[k] = -g_y_x_xxz_xy[k] * ab_z + g_y_x_xxz_xyz[k];

                g_y_x_xxzz_xz[k] = -g_y_x_xxz_xz[k] * ab_z + g_y_x_xxz_xzz[k];

                g_y_x_xxzz_yy[k] = -g_y_x_xxz_yy[k] * ab_z + g_y_x_xxz_yyz[k];

                g_y_x_xxzz_yz[k] = -g_y_x_xxz_yz[k] * ab_z + g_y_x_xxz_yzz[k];

                g_y_x_xxzz_zz[k] = -g_y_x_xxz_zz[k] * ab_z + g_y_x_xxz_zzz[k];
            }

            /// Set up 306-312 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyyy_xx = cbuffer.data(gd_geom_11_off + 306 * ccomps * dcomps);

            auto g_y_x_xyyy_xy = cbuffer.data(gd_geom_11_off + 307 * ccomps * dcomps);

            auto g_y_x_xyyy_xz = cbuffer.data(gd_geom_11_off + 308 * ccomps * dcomps);

            auto g_y_x_xyyy_yy = cbuffer.data(gd_geom_11_off + 309 * ccomps * dcomps);

            auto g_y_x_xyyy_yz = cbuffer.data(gd_geom_11_off + 310 * ccomps * dcomps);

            auto g_y_x_xyyy_zz = cbuffer.data(gd_geom_11_off + 311 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyy_xx, g_y_0_yyy_xy, g_y_0_yyy_xz, g_y_0_yyy_yy, g_y_0_yyy_yz, g_y_0_yyy_zz, g_y_x_xyyy_xx, g_y_x_xyyy_xy, g_y_x_xyyy_xz, g_y_x_xyyy_yy, g_y_x_xyyy_yz, g_y_x_xyyy_zz, g_y_x_yyy_xx, g_y_x_yyy_xxx, g_y_x_yyy_xxy, g_y_x_yyy_xxz, g_y_x_yyy_xy, g_y_x_yyy_xyy, g_y_x_yyy_xyz, g_y_x_yyy_xz, g_y_x_yyy_xzz, g_y_x_yyy_yy, g_y_x_yyy_yz, g_y_x_yyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyyy_xx[k] = g_y_0_yyy_xx[k] - g_y_x_yyy_xx[k] * ab_x + g_y_x_yyy_xxx[k];

                g_y_x_xyyy_xy[k] = g_y_0_yyy_xy[k] - g_y_x_yyy_xy[k] * ab_x + g_y_x_yyy_xxy[k];

                g_y_x_xyyy_xz[k] = g_y_0_yyy_xz[k] - g_y_x_yyy_xz[k] * ab_x + g_y_x_yyy_xxz[k];

                g_y_x_xyyy_yy[k] = g_y_0_yyy_yy[k] - g_y_x_yyy_yy[k] * ab_x + g_y_x_yyy_xyy[k];

                g_y_x_xyyy_yz[k] = g_y_0_yyy_yz[k] - g_y_x_yyy_yz[k] * ab_x + g_y_x_yyy_xyz[k];

                g_y_x_xyyy_zz[k] = g_y_0_yyy_zz[k] - g_y_x_yyy_zz[k] * ab_x + g_y_x_yyy_xzz[k];
            }

            /// Set up 312-318 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyyz_xx = cbuffer.data(gd_geom_11_off + 312 * ccomps * dcomps);

            auto g_y_x_xyyz_xy = cbuffer.data(gd_geom_11_off + 313 * ccomps * dcomps);

            auto g_y_x_xyyz_xz = cbuffer.data(gd_geom_11_off + 314 * ccomps * dcomps);

            auto g_y_x_xyyz_yy = cbuffer.data(gd_geom_11_off + 315 * ccomps * dcomps);

            auto g_y_x_xyyz_yz = cbuffer.data(gd_geom_11_off + 316 * ccomps * dcomps);

            auto g_y_x_xyyz_zz = cbuffer.data(gd_geom_11_off + 317 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xyy_xx, g_y_x_xyy_xxz, g_y_x_xyy_xy, g_y_x_xyy_xyz, g_y_x_xyy_xz, g_y_x_xyy_xzz, g_y_x_xyy_yy, g_y_x_xyy_yyz, g_y_x_xyy_yz, g_y_x_xyy_yzz, g_y_x_xyy_zz, g_y_x_xyy_zzz, g_y_x_xyyz_xx, g_y_x_xyyz_xy, g_y_x_xyyz_xz, g_y_x_xyyz_yy, g_y_x_xyyz_yz, g_y_x_xyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyyz_xx[k] = -g_y_x_xyy_xx[k] * ab_z + g_y_x_xyy_xxz[k];

                g_y_x_xyyz_xy[k] = -g_y_x_xyy_xy[k] * ab_z + g_y_x_xyy_xyz[k];

                g_y_x_xyyz_xz[k] = -g_y_x_xyy_xz[k] * ab_z + g_y_x_xyy_xzz[k];

                g_y_x_xyyz_yy[k] = -g_y_x_xyy_yy[k] * ab_z + g_y_x_xyy_yyz[k];

                g_y_x_xyyz_yz[k] = -g_y_x_xyy_yz[k] * ab_z + g_y_x_xyy_yzz[k];

                g_y_x_xyyz_zz[k] = -g_y_x_xyy_zz[k] * ab_z + g_y_x_xyy_zzz[k];
            }

            /// Set up 318-324 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyzz_xx = cbuffer.data(gd_geom_11_off + 318 * ccomps * dcomps);

            auto g_y_x_xyzz_xy = cbuffer.data(gd_geom_11_off + 319 * ccomps * dcomps);

            auto g_y_x_xyzz_xz = cbuffer.data(gd_geom_11_off + 320 * ccomps * dcomps);

            auto g_y_x_xyzz_yy = cbuffer.data(gd_geom_11_off + 321 * ccomps * dcomps);

            auto g_y_x_xyzz_yz = cbuffer.data(gd_geom_11_off + 322 * ccomps * dcomps);

            auto g_y_x_xyzz_zz = cbuffer.data(gd_geom_11_off + 323 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xyz_xx, g_y_x_xyz_xxz, g_y_x_xyz_xy, g_y_x_xyz_xyz, g_y_x_xyz_xz, g_y_x_xyz_xzz, g_y_x_xyz_yy, g_y_x_xyz_yyz, g_y_x_xyz_yz, g_y_x_xyz_yzz, g_y_x_xyz_zz, g_y_x_xyz_zzz, g_y_x_xyzz_xx, g_y_x_xyzz_xy, g_y_x_xyzz_xz, g_y_x_xyzz_yy, g_y_x_xyzz_yz, g_y_x_xyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyzz_xx[k] = -g_y_x_xyz_xx[k] * ab_z + g_y_x_xyz_xxz[k];

                g_y_x_xyzz_xy[k] = -g_y_x_xyz_xy[k] * ab_z + g_y_x_xyz_xyz[k];

                g_y_x_xyzz_xz[k] = -g_y_x_xyz_xz[k] * ab_z + g_y_x_xyz_xzz[k];

                g_y_x_xyzz_yy[k] = -g_y_x_xyz_yy[k] * ab_z + g_y_x_xyz_yyz[k];

                g_y_x_xyzz_yz[k] = -g_y_x_xyz_yz[k] * ab_z + g_y_x_xyz_yzz[k];

                g_y_x_xyzz_zz[k] = -g_y_x_xyz_zz[k] * ab_z + g_y_x_xyz_zzz[k];
            }

            /// Set up 324-330 components of targeted buffer : cbuffer.data(

            auto g_y_x_xzzz_xx = cbuffer.data(gd_geom_11_off + 324 * ccomps * dcomps);

            auto g_y_x_xzzz_xy = cbuffer.data(gd_geom_11_off + 325 * ccomps * dcomps);

            auto g_y_x_xzzz_xz = cbuffer.data(gd_geom_11_off + 326 * ccomps * dcomps);

            auto g_y_x_xzzz_yy = cbuffer.data(gd_geom_11_off + 327 * ccomps * dcomps);

            auto g_y_x_xzzz_yz = cbuffer.data(gd_geom_11_off + 328 * ccomps * dcomps);

            auto g_y_x_xzzz_zz = cbuffer.data(gd_geom_11_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xzz_xx, g_y_x_xzz_xxz, g_y_x_xzz_xy, g_y_x_xzz_xyz, g_y_x_xzz_xz, g_y_x_xzz_xzz, g_y_x_xzz_yy, g_y_x_xzz_yyz, g_y_x_xzz_yz, g_y_x_xzz_yzz, g_y_x_xzz_zz, g_y_x_xzz_zzz, g_y_x_xzzz_xx, g_y_x_xzzz_xy, g_y_x_xzzz_xz, g_y_x_xzzz_yy, g_y_x_xzzz_yz, g_y_x_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xzzz_xx[k] = -g_y_x_xzz_xx[k] * ab_z + g_y_x_xzz_xxz[k];

                g_y_x_xzzz_xy[k] = -g_y_x_xzz_xy[k] * ab_z + g_y_x_xzz_xyz[k];

                g_y_x_xzzz_xz[k] = -g_y_x_xzz_xz[k] * ab_z + g_y_x_xzz_xzz[k];

                g_y_x_xzzz_yy[k] = -g_y_x_xzz_yy[k] * ab_z + g_y_x_xzz_yyz[k];

                g_y_x_xzzz_yz[k] = -g_y_x_xzz_yz[k] * ab_z + g_y_x_xzz_yzz[k];

                g_y_x_xzzz_zz[k] = -g_y_x_xzz_zz[k] * ab_z + g_y_x_xzz_zzz[k];
            }

            /// Set up 330-336 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyyy_xx = cbuffer.data(gd_geom_11_off + 330 * ccomps * dcomps);

            auto g_y_x_yyyy_xy = cbuffer.data(gd_geom_11_off + 331 * ccomps * dcomps);

            auto g_y_x_yyyy_xz = cbuffer.data(gd_geom_11_off + 332 * ccomps * dcomps);

            auto g_y_x_yyyy_yy = cbuffer.data(gd_geom_11_off + 333 * ccomps * dcomps);

            auto g_y_x_yyyy_yz = cbuffer.data(gd_geom_11_off + 334 * ccomps * dcomps);

            auto g_y_x_yyyy_zz = cbuffer.data(gd_geom_11_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyy_xx, g_0_x_yyy_xy, g_0_x_yyy_xz, g_0_x_yyy_yy, g_0_x_yyy_yz, g_0_x_yyy_zz, g_y_x_yyy_xx, g_y_x_yyy_xxy, g_y_x_yyy_xy, g_y_x_yyy_xyy, g_y_x_yyy_xyz, g_y_x_yyy_xz, g_y_x_yyy_yy, g_y_x_yyy_yyy, g_y_x_yyy_yyz, g_y_x_yyy_yz, g_y_x_yyy_yzz, g_y_x_yyy_zz, g_y_x_yyyy_xx, g_y_x_yyyy_xy, g_y_x_yyyy_xz, g_y_x_yyyy_yy, g_y_x_yyyy_yz, g_y_x_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyyy_xx[k] = -g_0_x_yyy_xx[k] - g_y_x_yyy_xx[k] * ab_y + g_y_x_yyy_xxy[k];

                g_y_x_yyyy_xy[k] = -g_0_x_yyy_xy[k] - g_y_x_yyy_xy[k] * ab_y + g_y_x_yyy_xyy[k];

                g_y_x_yyyy_xz[k] = -g_0_x_yyy_xz[k] - g_y_x_yyy_xz[k] * ab_y + g_y_x_yyy_xyz[k];

                g_y_x_yyyy_yy[k] = -g_0_x_yyy_yy[k] - g_y_x_yyy_yy[k] * ab_y + g_y_x_yyy_yyy[k];

                g_y_x_yyyy_yz[k] = -g_0_x_yyy_yz[k] - g_y_x_yyy_yz[k] * ab_y + g_y_x_yyy_yyz[k];

                g_y_x_yyyy_zz[k] = -g_0_x_yyy_zz[k] - g_y_x_yyy_zz[k] * ab_y + g_y_x_yyy_yzz[k];
            }

            /// Set up 336-342 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyyz_xx = cbuffer.data(gd_geom_11_off + 336 * ccomps * dcomps);

            auto g_y_x_yyyz_xy = cbuffer.data(gd_geom_11_off + 337 * ccomps * dcomps);

            auto g_y_x_yyyz_xz = cbuffer.data(gd_geom_11_off + 338 * ccomps * dcomps);

            auto g_y_x_yyyz_yy = cbuffer.data(gd_geom_11_off + 339 * ccomps * dcomps);

            auto g_y_x_yyyz_yz = cbuffer.data(gd_geom_11_off + 340 * ccomps * dcomps);

            auto g_y_x_yyyz_zz = cbuffer.data(gd_geom_11_off + 341 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yyy_xx, g_y_x_yyy_xxz, g_y_x_yyy_xy, g_y_x_yyy_xyz, g_y_x_yyy_xz, g_y_x_yyy_xzz, g_y_x_yyy_yy, g_y_x_yyy_yyz, g_y_x_yyy_yz, g_y_x_yyy_yzz, g_y_x_yyy_zz, g_y_x_yyy_zzz, g_y_x_yyyz_xx, g_y_x_yyyz_xy, g_y_x_yyyz_xz, g_y_x_yyyz_yy, g_y_x_yyyz_yz, g_y_x_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyyz_xx[k] = -g_y_x_yyy_xx[k] * ab_z + g_y_x_yyy_xxz[k];

                g_y_x_yyyz_xy[k] = -g_y_x_yyy_xy[k] * ab_z + g_y_x_yyy_xyz[k];

                g_y_x_yyyz_xz[k] = -g_y_x_yyy_xz[k] * ab_z + g_y_x_yyy_xzz[k];

                g_y_x_yyyz_yy[k] = -g_y_x_yyy_yy[k] * ab_z + g_y_x_yyy_yyz[k];

                g_y_x_yyyz_yz[k] = -g_y_x_yyy_yz[k] * ab_z + g_y_x_yyy_yzz[k];

                g_y_x_yyyz_zz[k] = -g_y_x_yyy_zz[k] * ab_z + g_y_x_yyy_zzz[k];
            }

            /// Set up 342-348 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyzz_xx = cbuffer.data(gd_geom_11_off + 342 * ccomps * dcomps);

            auto g_y_x_yyzz_xy = cbuffer.data(gd_geom_11_off + 343 * ccomps * dcomps);

            auto g_y_x_yyzz_xz = cbuffer.data(gd_geom_11_off + 344 * ccomps * dcomps);

            auto g_y_x_yyzz_yy = cbuffer.data(gd_geom_11_off + 345 * ccomps * dcomps);

            auto g_y_x_yyzz_yz = cbuffer.data(gd_geom_11_off + 346 * ccomps * dcomps);

            auto g_y_x_yyzz_zz = cbuffer.data(gd_geom_11_off + 347 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yyz_xx, g_y_x_yyz_xxz, g_y_x_yyz_xy, g_y_x_yyz_xyz, g_y_x_yyz_xz, g_y_x_yyz_xzz, g_y_x_yyz_yy, g_y_x_yyz_yyz, g_y_x_yyz_yz, g_y_x_yyz_yzz, g_y_x_yyz_zz, g_y_x_yyz_zzz, g_y_x_yyzz_xx, g_y_x_yyzz_xy, g_y_x_yyzz_xz, g_y_x_yyzz_yy, g_y_x_yyzz_yz, g_y_x_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyzz_xx[k] = -g_y_x_yyz_xx[k] * ab_z + g_y_x_yyz_xxz[k];

                g_y_x_yyzz_xy[k] = -g_y_x_yyz_xy[k] * ab_z + g_y_x_yyz_xyz[k];

                g_y_x_yyzz_xz[k] = -g_y_x_yyz_xz[k] * ab_z + g_y_x_yyz_xzz[k];

                g_y_x_yyzz_yy[k] = -g_y_x_yyz_yy[k] * ab_z + g_y_x_yyz_yyz[k];

                g_y_x_yyzz_yz[k] = -g_y_x_yyz_yz[k] * ab_z + g_y_x_yyz_yzz[k];

                g_y_x_yyzz_zz[k] = -g_y_x_yyz_zz[k] * ab_z + g_y_x_yyz_zzz[k];
            }

            /// Set up 348-354 components of targeted buffer : cbuffer.data(

            auto g_y_x_yzzz_xx = cbuffer.data(gd_geom_11_off + 348 * ccomps * dcomps);

            auto g_y_x_yzzz_xy = cbuffer.data(gd_geom_11_off + 349 * ccomps * dcomps);

            auto g_y_x_yzzz_xz = cbuffer.data(gd_geom_11_off + 350 * ccomps * dcomps);

            auto g_y_x_yzzz_yy = cbuffer.data(gd_geom_11_off + 351 * ccomps * dcomps);

            auto g_y_x_yzzz_yz = cbuffer.data(gd_geom_11_off + 352 * ccomps * dcomps);

            auto g_y_x_yzzz_zz = cbuffer.data(gd_geom_11_off + 353 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yzz_xx, g_y_x_yzz_xxz, g_y_x_yzz_xy, g_y_x_yzz_xyz, g_y_x_yzz_xz, g_y_x_yzz_xzz, g_y_x_yzz_yy, g_y_x_yzz_yyz, g_y_x_yzz_yz, g_y_x_yzz_yzz, g_y_x_yzz_zz, g_y_x_yzz_zzz, g_y_x_yzzz_xx, g_y_x_yzzz_xy, g_y_x_yzzz_xz, g_y_x_yzzz_yy, g_y_x_yzzz_yz, g_y_x_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yzzz_xx[k] = -g_y_x_yzz_xx[k] * ab_z + g_y_x_yzz_xxz[k];

                g_y_x_yzzz_xy[k] = -g_y_x_yzz_xy[k] * ab_z + g_y_x_yzz_xyz[k];

                g_y_x_yzzz_xz[k] = -g_y_x_yzz_xz[k] * ab_z + g_y_x_yzz_xzz[k];

                g_y_x_yzzz_yy[k] = -g_y_x_yzz_yy[k] * ab_z + g_y_x_yzz_yyz[k];

                g_y_x_yzzz_yz[k] = -g_y_x_yzz_yz[k] * ab_z + g_y_x_yzz_yzz[k];

                g_y_x_yzzz_zz[k] = -g_y_x_yzz_zz[k] * ab_z + g_y_x_yzz_zzz[k];
            }

            /// Set up 354-360 components of targeted buffer : cbuffer.data(

            auto g_y_x_zzzz_xx = cbuffer.data(gd_geom_11_off + 354 * ccomps * dcomps);

            auto g_y_x_zzzz_xy = cbuffer.data(gd_geom_11_off + 355 * ccomps * dcomps);

            auto g_y_x_zzzz_xz = cbuffer.data(gd_geom_11_off + 356 * ccomps * dcomps);

            auto g_y_x_zzzz_yy = cbuffer.data(gd_geom_11_off + 357 * ccomps * dcomps);

            auto g_y_x_zzzz_yz = cbuffer.data(gd_geom_11_off + 358 * ccomps * dcomps);

            auto g_y_x_zzzz_zz = cbuffer.data(gd_geom_11_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_zzz_xx, g_y_x_zzz_xxz, g_y_x_zzz_xy, g_y_x_zzz_xyz, g_y_x_zzz_xz, g_y_x_zzz_xzz, g_y_x_zzz_yy, g_y_x_zzz_yyz, g_y_x_zzz_yz, g_y_x_zzz_yzz, g_y_x_zzz_zz, g_y_x_zzz_zzz, g_y_x_zzzz_xx, g_y_x_zzzz_xy, g_y_x_zzzz_xz, g_y_x_zzzz_yy, g_y_x_zzzz_yz, g_y_x_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_zzzz_xx[k] = -g_y_x_zzz_xx[k] * ab_z + g_y_x_zzz_xxz[k];

                g_y_x_zzzz_xy[k] = -g_y_x_zzz_xy[k] * ab_z + g_y_x_zzz_xyz[k];

                g_y_x_zzzz_xz[k] = -g_y_x_zzz_xz[k] * ab_z + g_y_x_zzz_xzz[k];

                g_y_x_zzzz_yy[k] = -g_y_x_zzz_yy[k] * ab_z + g_y_x_zzz_yyz[k];

                g_y_x_zzzz_yz[k] = -g_y_x_zzz_yz[k] * ab_z + g_y_x_zzz_yzz[k];

                g_y_x_zzzz_zz[k] = -g_y_x_zzz_zz[k] * ab_z + g_y_x_zzz_zzz[k];
            }

            /// Set up 360-366 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxxx_xx = cbuffer.data(gd_geom_11_off + 360 * ccomps * dcomps);

            auto g_y_y_xxxx_xy = cbuffer.data(gd_geom_11_off + 361 * ccomps * dcomps);

            auto g_y_y_xxxx_xz = cbuffer.data(gd_geom_11_off + 362 * ccomps * dcomps);

            auto g_y_y_xxxx_yy = cbuffer.data(gd_geom_11_off + 363 * ccomps * dcomps);

            auto g_y_y_xxxx_yz = cbuffer.data(gd_geom_11_off + 364 * ccomps * dcomps);

            auto g_y_y_xxxx_zz = cbuffer.data(gd_geom_11_off + 365 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxx_xx, g_y_y_xxx_xxx, g_y_y_xxx_xxy, g_y_y_xxx_xxz, g_y_y_xxx_xy, g_y_y_xxx_xyy, g_y_y_xxx_xyz, g_y_y_xxx_xz, g_y_y_xxx_xzz, g_y_y_xxx_yy, g_y_y_xxx_yz, g_y_y_xxx_zz, g_y_y_xxxx_xx, g_y_y_xxxx_xy, g_y_y_xxxx_xz, g_y_y_xxxx_yy, g_y_y_xxxx_yz, g_y_y_xxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxxx_xx[k] = -g_y_y_xxx_xx[k] * ab_x + g_y_y_xxx_xxx[k];

                g_y_y_xxxx_xy[k] = -g_y_y_xxx_xy[k] * ab_x + g_y_y_xxx_xxy[k];

                g_y_y_xxxx_xz[k] = -g_y_y_xxx_xz[k] * ab_x + g_y_y_xxx_xxz[k];

                g_y_y_xxxx_yy[k] = -g_y_y_xxx_yy[k] * ab_x + g_y_y_xxx_xyy[k];

                g_y_y_xxxx_yz[k] = -g_y_y_xxx_yz[k] * ab_x + g_y_y_xxx_xyz[k];

                g_y_y_xxxx_zz[k] = -g_y_y_xxx_zz[k] * ab_x + g_y_y_xxx_xzz[k];
            }

            /// Set up 366-372 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxxy_xx = cbuffer.data(gd_geom_11_off + 366 * ccomps * dcomps);

            auto g_y_y_xxxy_xy = cbuffer.data(gd_geom_11_off + 367 * ccomps * dcomps);

            auto g_y_y_xxxy_xz = cbuffer.data(gd_geom_11_off + 368 * ccomps * dcomps);

            auto g_y_y_xxxy_yy = cbuffer.data(gd_geom_11_off + 369 * ccomps * dcomps);

            auto g_y_y_xxxy_yz = cbuffer.data(gd_geom_11_off + 370 * ccomps * dcomps);

            auto g_y_y_xxxy_zz = cbuffer.data(gd_geom_11_off + 371 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxxy_xx, g_y_y_xxxy_xy, g_y_y_xxxy_xz, g_y_y_xxxy_yy, g_y_y_xxxy_yz, g_y_y_xxxy_zz, g_y_y_xxy_xx, g_y_y_xxy_xxx, g_y_y_xxy_xxy, g_y_y_xxy_xxz, g_y_y_xxy_xy, g_y_y_xxy_xyy, g_y_y_xxy_xyz, g_y_y_xxy_xz, g_y_y_xxy_xzz, g_y_y_xxy_yy, g_y_y_xxy_yz, g_y_y_xxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxxy_xx[k] = -g_y_y_xxy_xx[k] * ab_x + g_y_y_xxy_xxx[k];

                g_y_y_xxxy_xy[k] = -g_y_y_xxy_xy[k] * ab_x + g_y_y_xxy_xxy[k];

                g_y_y_xxxy_xz[k] = -g_y_y_xxy_xz[k] * ab_x + g_y_y_xxy_xxz[k];

                g_y_y_xxxy_yy[k] = -g_y_y_xxy_yy[k] * ab_x + g_y_y_xxy_xyy[k];

                g_y_y_xxxy_yz[k] = -g_y_y_xxy_yz[k] * ab_x + g_y_y_xxy_xyz[k];

                g_y_y_xxxy_zz[k] = -g_y_y_xxy_zz[k] * ab_x + g_y_y_xxy_xzz[k];
            }

            /// Set up 372-378 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxxz_xx = cbuffer.data(gd_geom_11_off + 372 * ccomps * dcomps);

            auto g_y_y_xxxz_xy = cbuffer.data(gd_geom_11_off + 373 * ccomps * dcomps);

            auto g_y_y_xxxz_xz = cbuffer.data(gd_geom_11_off + 374 * ccomps * dcomps);

            auto g_y_y_xxxz_yy = cbuffer.data(gd_geom_11_off + 375 * ccomps * dcomps);

            auto g_y_y_xxxz_yz = cbuffer.data(gd_geom_11_off + 376 * ccomps * dcomps);

            auto g_y_y_xxxz_zz = cbuffer.data(gd_geom_11_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxxz_xx, g_y_y_xxxz_xy, g_y_y_xxxz_xz, g_y_y_xxxz_yy, g_y_y_xxxz_yz, g_y_y_xxxz_zz, g_y_y_xxz_xx, g_y_y_xxz_xxx, g_y_y_xxz_xxy, g_y_y_xxz_xxz, g_y_y_xxz_xy, g_y_y_xxz_xyy, g_y_y_xxz_xyz, g_y_y_xxz_xz, g_y_y_xxz_xzz, g_y_y_xxz_yy, g_y_y_xxz_yz, g_y_y_xxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxxz_xx[k] = -g_y_y_xxz_xx[k] * ab_x + g_y_y_xxz_xxx[k];

                g_y_y_xxxz_xy[k] = -g_y_y_xxz_xy[k] * ab_x + g_y_y_xxz_xxy[k];

                g_y_y_xxxz_xz[k] = -g_y_y_xxz_xz[k] * ab_x + g_y_y_xxz_xxz[k];

                g_y_y_xxxz_yy[k] = -g_y_y_xxz_yy[k] * ab_x + g_y_y_xxz_xyy[k];

                g_y_y_xxxz_yz[k] = -g_y_y_xxz_yz[k] * ab_x + g_y_y_xxz_xyz[k];

                g_y_y_xxxz_zz[k] = -g_y_y_xxz_zz[k] * ab_x + g_y_y_xxz_xzz[k];
            }

            /// Set up 378-384 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxyy_xx = cbuffer.data(gd_geom_11_off + 378 * ccomps * dcomps);

            auto g_y_y_xxyy_xy = cbuffer.data(gd_geom_11_off + 379 * ccomps * dcomps);

            auto g_y_y_xxyy_xz = cbuffer.data(gd_geom_11_off + 380 * ccomps * dcomps);

            auto g_y_y_xxyy_yy = cbuffer.data(gd_geom_11_off + 381 * ccomps * dcomps);

            auto g_y_y_xxyy_yz = cbuffer.data(gd_geom_11_off + 382 * ccomps * dcomps);

            auto g_y_y_xxyy_zz = cbuffer.data(gd_geom_11_off + 383 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxyy_xx, g_y_y_xxyy_xy, g_y_y_xxyy_xz, g_y_y_xxyy_yy, g_y_y_xxyy_yz, g_y_y_xxyy_zz, g_y_y_xyy_xx, g_y_y_xyy_xxx, g_y_y_xyy_xxy, g_y_y_xyy_xxz, g_y_y_xyy_xy, g_y_y_xyy_xyy, g_y_y_xyy_xyz, g_y_y_xyy_xz, g_y_y_xyy_xzz, g_y_y_xyy_yy, g_y_y_xyy_yz, g_y_y_xyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxyy_xx[k] = -g_y_y_xyy_xx[k] * ab_x + g_y_y_xyy_xxx[k];

                g_y_y_xxyy_xy[k] = -g_y_y_xyy_xy[k] * ab_x + g_y_y_xyy_xxy[k];

                g_y_y_xxyy_xz[k] = -g_y_y_xyy_xz[k] * ab_x + g_y_y_xyy_xxz[k];

                g_y_y_xxyy_yy[k] = -g_y_y_xyy_yy[k] * ab_x + g_y_y_xyy_xyy[k];

                g_y_y_xxyy_yz[k] = -g_y_y_xyy_yz[k] * ab_x + g_y_y_xyy_xyz[k];

                g_y_y_xxyy_zz[k] = -g_y_y_xyy_zz[k] * ab_x + g_y_y_xyy_xzz[k];
            }

            /// Set up 384-390 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxyz_xx = cbuffer.data(gd_geom_11_off + 384 * ccomps * dcomps);

            auto g_y_y_xxyz_xy = cbuffer.data(gd_geom_11_off + 385 * ccomps * dcomps);

            auto g_y_y_xxyz_xz = cbuffer.data(gd_geom_11_off + 386 * ccomps * dcomps);

            auto g_y_y_xxyz_yy = cbuffer.data(gd_geom_11_off + 387 * ccomps * dcomps);

            auto g_y_y_xxyz_yz = cbuffer.data(gd_geom_11_off + 388 * ccomps * dcomps);

            auto g_y_y_xxyz_zz = cbuffer.data(gd_geom_11_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxyz_xx, g_y_y_xxyz_xy, g_y_y_xxyz_xz, g_y_y_xxyz_yy, g_y_y_xxyz_yz, g_y_y_xxyz_zz, g_y_y_xyz_xx, g_y_y_xyz_xxx, g_y_y_xyz_xxy, g_y_y_xyz_xxz, g_y_y_xyz_xy, g_y_y_xyz_xyy, g_y_y_xyz_xyz, g_y_y_xyz_xz, g_y_y_xyz_xzz, g_y_y_xyz_yy, g_y_y_xyz_yz, g_y_y_xyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxyz_xx[k] = -g_y_y_xyz_xx[k] * ab_x + g_y_y_xyz_xxx[k];

                g_y_y_xxyz_xy[k] = -g_y_y_xyz_xy[k] * ab_x + g_y_y_xyz_xxy[k];

                g_y_y_xxyz_xz[k] = -g_y_y_xyz_xz[k] * ab_x + g_y_y_xyz_xxz[k];

                g_y_y_xxyz_yy[k] = -g_y_y_xyz_yy[k] * ab_x + g_y_y_xyz_xyy[k];

                g_y_y_xxyz_yz[k] = -g_y_y_xyz_yz[k] * ab_x + g_y_y_xyz_xyz[k];

                g_y_y_xxyz_zz[k] = -g_y_y_xyz_zz[k] * ab_x + g_y_y_xyz_xzz[k];
            }

            /// Set up 390-396 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxzz_xx = cbuffer.data(gd_geom_11_off + 390 * ccomps * dcomps);

            auto g_y_y_xxzz_xy = cbuffer.data(gd_geom_11_off + 391 * ccomps * dcomps);

            auto g_y_y_xxzz_xz = cbuffer.data(gd_geom_11_off + 392 * ccomps * dcomps);

            auto g_y_y_xxzz_yy = cbuffer.data(gd_geom_11_off + 393 * ccomps * dcomps);

            auto g_y_y_xxzz_yz = cbuffer.data(gd_geom_11_off + 394 * ccomps * dcomps);

            auto g_y_y_xxzz_zz = cbuffer.data(gd_geom_11_off + 395 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxzz_xx, g_y_y_xxzz_xy, g_y_y_xxzz_xz, g_y_y_xxzz_yy, g_y_y_xxzz_yz, g_y_y_xxzz_zz, g_y_y_xzz_xx, g_y_y_xzz_xxx, g_y_y_xzz_xxy, g_y_y_xzz_xxz, g_y_y_xzz_xy, g_y_y_xzz_xyy, g_y_y_xzz_xyz, g_y_y_xzz_xz, g_y_y_xzz_xzz, g_y_y_xzz_yy, g_y_y_xzz_yz, g_y_y_xzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxzz_xx[k] = -g_y_y_xzz_xx[k] * ab_x + g_y_y_xzz_xxx[k];

                g_y_y_xxzz_xy[k] = -g_y_y_xzz_xy[k] * ab_x + g_y_y_xzz_xxy[k];

                g_y_y_xxzz_xz[k] = -g_y_y_xzz_xz[k] * ab_x + g_y_y_xzz_xxz[k];

                g_y_y_xxzz_yy[k] = -g_y_y_xzz_yy[k] * ab_x + g_y_y_xzz_xyy[k];

                g_y_y_xxzz_yz[k] = -g_y_y_xzz_yz[k] * ab_x + g_y_y_xzz_xyz[k];

                g_y_y_xxzz_zz[k] = -g_y_y_xzz_zz[k] * ab_x + g_y_y_xzz_xzz[k];
            }

            /// Set up 396-402 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyyy_xx = cbuffer.data(gd_geom_11_off + 396 * ccomps * dcomps);

            auto g_y_y_xyyy_xy = cbuffer.data(gd_geom_11_off + 397 * ccomps * dcomps);

            auto g_y_y_xyyy_xz = cbuffer.data(gd_geom_11_off + 398 * ccomps * dcomps);

            auto g_y_y_xyyy_yy = cbuffer.data(gd_geom_11_off + 399 * ccomps * dcomps);

            auto g_y_y_xyyy_yz = cbuffer.data(gd_geom_11_off + 400 * ccomps * dcomps);

            auto g_y_y_xyyy_zz = cbuffer.data(gd_geom_11_off + 401 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyyy_xx, g_y_y_xyyy_xy, g_y_y_xyyy_xz, g_y_y_xyyy_yy, g_y_y_xyyy_yz, g_y_y_xyyy_zz, g_y_y_yyy_xx, g_y_y_yyy_xxx, g_y_y_yyy_xxy, g_y_y_yyy_xxz, g_y_y_yyy_xy, g_y_y_yyy_xyy, g_y_y_yyy_xyz, g_y_y_yyy_xz, g_y_y_yyy_xzz, g_y_y_yyy_yy, g_y_y_yyy_yz, g_y_y_yyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyyy_xx[k] = -g_y_y_yyy_xx[k] * ab_x + g_y_y_yyy_xxx[k];

                g_y_y_xyyy_xy[k] = -g_y_y_yyy_xy[k] * ab_x + g_y_y_yyy_xxy[k];

                g_y_y_xyyy_xz[k] = -g_y_y_yyy_xz[k] * ab_x + g_y_y_yyy_xxz[k];

                g_y_y_xyyy_yy[k] = -g_y_y_yyy_yy[k] * ab_x + g_y_y_yyy_xyy[k];

                g_y_y_xyyy_yz[k] = -g_y_y_yyy_yz[k] * ab_x + g_y_y_yyy_xyz[k];

                g_y_y_xyyy_zz[k] = -g_y_y_yyy_zz[k] * ab_x + g_y_y_yyy_xzz[k];
            }

            /// Set up 402-408 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyyz_xx = cbuffer.data(gd_geom_11_off + 402 * ccomps * dcomps);

            auto g_y_y_xyyz_xy = cbuffer.data(gd_geom_11_off + 403 * ccomps * dcomps);

            auto g_y_y_xyyz_xz = cbuffer.data(gd_geom_11_off + 404 * ccomps * dcomps);

            auto g_y_y_xyyz_yy = cbuffer.data(gd_geom_11_off + 405 * ccomps * dcomps);

            auto g_y_y_xyyz_yz = cbuffer.data(gd_geom_11_off + 406 * ccomps * dcomps);

            auto g_y_y_xyyz_zz = cbuffer.data(gd_geom_11_off + 407 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyyz_xx, g_y_y_xyyz_xy, g_y_y_xyyz_xz, g_y_y_xyyz_yy, g_y_y_xyyz_yz, g_y_y_xyyz_zz, g_y_y_yyz_xx, g_y_y_yyz_xxx, g_y_y_yyz_xxy, g_y_y_yyz_xxz, g_y_y_yyz_xy, g_y_y_yyz_xyy, g_y_y_yyz_xyz, g_y_y_yyz_xz, g_y_y_yyz_xzz, g_y_y_yyz_yy, g_y_y_yyz_yz, g_y_y_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyyz_xx[k] = -g_y_y_yyz_xx[k] * ab_x + g_y_y_yyz_xxx[k];

                g_y_y_xyyz_xy[k] = -g_y_y_yyz_xy[k] * ab_x + g_y_y_yyz_xxy[k];

                g_y_y_xyyz_xz[k] = -g_y_y_yyz_xz[k] * ab_x + g_y_y_yyz_xxz[k];

                g_y_y_xyyz_yy[k] = -g_y_y_yyz_yy[k] * ab_x + g_y_y_yyz_xyy[k];

                g_y_y_xyyz_yz[k] = -g_y_y_yyz_yz[k] * ab_x + g_y_y_yyz_xyz[k];

                g_y_y_xyyz_zz[k] = -g_y_y_yyz_zz[k] * ab_x + g_y_y_yyz_xzz[k];
            }

            /// Set up 408-414 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyzz_xx = cbuffer.data(gd_geom_11_off + 408 * ccomps * dcomps);

            auto g_y_y_xyzz_xy = cbuffer.data(gd_geom_11_off + 409 * ccomps * dcomps);

            auto g_y_y_xyzz_xz = cbuffer.data(gd_geom_11_off + 410 * ccomps * dcomps);

            auto g_y_y_xyzz_yy = cbuffer.data(gd_geom_11_off + 411 * ccomps * dcomps);

            auto g_y_y_xyzz_yz = cbuffer.data(gd_geom_11_off + 412 * ccomps * dcomps);

            auto g_y_y_xyzz_zz = cbuffer.data(gd_geom_11_off + 413 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyzz_xx, g_y_y_xyzz_xy, g_y_y_xyzz_xz, g_y_y_xyzz_yy, g_y_y_xyzz_yz, g_y_y_xyzz_zz, g_y_y_yzz_xx, g_y_y_yzz_xxx, g_y_y_yzz_xxy, g_y_y_yzz_xxz, g_y_y_yzz_xy, g_y_y_yzz_xyy, g_y_y_yzz_xyz, g_y_y_yzz_xz, g_y_y_yzz_xzz, g_y_y_yzz_yy, g_y_y_yzz_yz, g_y_y_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyzz_xx[k] = -g_y_y_yzz_xx[k] * ab_x + g_y_y_yzz_xxx[k];

                g_y_y_xyzz_xy[k] = -g_y_y_yzz_xy[k] * ab_x + g_y_y_yzz_xxy[k];

                g_y_y_xyzz_xz[k] = -g_y_y_yzz_xz[k] * ab_x + g_y_y_yzz_xxz[k];

                g_y_y_xyzz_yy[k] = -g_y_y_yzz_yy[k] * ab_x + g_y_y_yzz_xyy[k];

                g_y_y_xyzz_yz[k] = -g_y_y_yzz_yz[k] * ab_x + g_y_y_yzz_xyz[k];

                g_y_y_xyzz_zz[k] = -g_y_y_yzz_zz[k] * ab_x + g_y_y_yzz_xzz[k];
            }

            /// Set up 414-420 components of targeted buffer : cbuffer.data(

            auto g_y_y_xzzz_xx = cbuffer.data(gd_geom_11_off + 414 * ccomps * dcomps);

            auto g_y_y_xzzz_xy = cbuffer.data(gd_geom_11_off + 415 * ccomps * dcomps);

            auto g_y_y_xzzz_xz = cbuffer.data(gd_geom_11_off + 416 * ccomps * dcomps);

            auto g_y_y_xzzz_yy = cbuffer.data(gd_geom_11_off + 417 * ccomps * dcomps);

            auto g_y_y_xzzz_yz = cbuffer.data(gd_geom_11_off + 418 * ccomps * dcomps);

            auto g_y_y_xzzz_zz = cbuffer.data(gd_geom_11_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xzzz_xx, g_y_y_xzzz_xy, g_y_y_xzzz_xz, g_y_y_xzzz_yy, g_y_y_xzzz_yz, g_y_y_xzzz_zz, g_y_y_zzz_xx, g_y_y_zzz_xxx, g_y_y_zzz_xxy, g_y_y_zzz_xxz, g_y_y_zzz_xy, g_y_y_zzz_xyy, g_y_y_zzz_xyz, g_y_y_zzz_xz, g_y_y_zzz_xzz, g_y_y_zzz_yy, g_y_y_zzz_yz, g_y_y_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xzzz_xx[k] = -g_y_y_zzz_xx[k] * ab_x + g_y_y_zzz_xxx[k];

                g_y_y_xzzz_xy[k] = -g_y_y_zzz_xy[k] * ab_x + g_y_y_zzz_xxy[k];

                g_y_y_xzzz_xz[k] = -g_y_y_zzz_xz[k] * ab_x + g_y_y_zzz_xxz[k];

                g_y_y_xzzz_yy[k] = -g_y_y_zzz_yy[k] * ab_x + g_y_y_zzz_xyy[k];

                g_y_y_xzzz_yz[k] = -g_y_y_zzz_yz[k] * ab_x + g_y_y_zzz_xyz[k];

                g_y_y_xzzz_zz[k] = -g_y_y_zzz_zz[k] * ab_x + g_y_y_zzz_xzz[k];
            }

            /// Set up 420-426 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyyy_xx = cbuffer.data(gd_geom_11_off + 420 * ccomps * dcomps);

            auto g_y_y_yyyy_xy = cbuffer.data(gd_geom_11_off + 421 * ccomps * dcomps);

            auto g_y_y_yyyy_xz = cbuffer.data(gd_geom_11_off + 422 * ccomps * dcomps);

            auto g_y_y_yyyy_yy = cbuffer.data(gd_geom_11_off + 423 * ccomps * dcomps);

            auto g_y_y_yyyy_yz = cbuffer.data(gd_geom_11_off + 424 * ccomps * dcomps);

            auto g_y_y_yyyy_zz = cbuffer.data(gd_geom_11_off + 425 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyy_xx, g_0_y_yyy_xy, g_0_y_yyy_xz, g_0_y_yyy_yy, g_0_y_yyy_yz, g_0_y_yyy_zz, g_y_0_yyy_xx, g_y_0_yyy_xy, g_y_0_yyy_xz, g_y_0_yyy_yy, g_y_0_yyy_yz, g_y_0_yyy_zz, g_y_y_yyy_xx, g_y_y_yyy_xxy, g_y_y_yyy_xy, g_y_y_yyy_xyy, g_y_y_yyy_xyz, g_y_y_yyy_xz, g_y_y_yyy_yy, g_y_y_yyy_yyy, g_y_y_yyy_yyz, g_y_y_yyy_yz, g_y_y_yyy_yzz, g_y_y_yyy_zz, g_y_y_yyyy_xx, g_y_y_yyyy_xy, g_y_y_yyyy_xz, g_y_y_yyyy_yy, g_y_y_yyyy_yz, g_y_y_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyyy_xx[k] = -g_0_y_yyy_xx[k] + g_y_0_yyy_xx[k] - g_y_y_yyy_xx[k] * ab_y + g_y_y_yyy_xxy[k];

                g_y_y_yyyy_xy[k] = -g_0_y_yyy_xy[k] + g_y_0_yyy_xy[k] - g_y_y_yyy_xy[k] * ab_y + g_y_y_yyy_xyy[k];

                g_y_y_yyyy_xz[k] = -g_0_y_yyy_xz[k] + g_y_0_yyy_xz[k] - g_y_y_yyy_xz[k] * ab_y + g_y_y_yyy_xyz[k];

                g_y_y_yyyy_yy[k] = -g_0_y_yyy_yy[k] + g_y_0_yyy_yy[k] - g_y_y_yyy_yy[k] * ab_y + g_y_y_yyy_yyy[k];

                g_y_y_yyyy_yz[k] = -g_0_y_yyy_yz[k] + g_y_0_yyy_yz[k] - g_y_y_yyy_yz[k] * ab_y + g_y_y_yyy_yyz[k];

                g_y_y_yyyy_zz[k] = -g_0_y_yyy_zz[k] + g_y_0_yyy_zz[k] - g_y_y_yyy_zz[k] * ab_y + g_y_y_yyy_yzz[k];
            }

            /// Set up 426-432 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyyz_xx = cbuffer.data(gd_geom_11_off + 426 * ccomps * dcomps);

            auto g_y_y_yyyz_xy = cbuffer.data(gd_geom_11_off + 427 * ccomps * dcomps);

            auto g_y_y_yyyz_xz = cbuffer.data(gd_geom_11_off + 428 * ccomps * dcomps);

            auto g_y_y_yyyz_yy = cbuffer.data(gd_geom_11_off + 429 * ccomps * dcomps);

            auto g_y_y_yyyz_yz = cbuffer.data(gd_geom_11_off + 430 * ccomps * dcomps);

            auto g_y_y_yyyz_zz = cbuffer.data(gd_geom_11_off + 431 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yyy_xx, g_y_y_yyy_xxz, g_y_y_yyy_xy, g_y_y_yyy_xyz, g_y_y_yyy_xz, g_y_y_yyy_xzz, g_y_y_yyy_yy, g_y_y_yyy_yyz, g_y_y_yyy_yz, g_y_y_yyy_yzz, g_y_y_yyy_zz, g_y_y_yyy_zzz, g_y_y_yyyz_xx, g_y_y_yyyz_xy, g_y_y_yyyz_xz, g_y_y_yyyz_yy, g_y_y_yyyz_yz, g_y_y_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyyz_xx[k] = -g_y_y_yyy_xx[k] * ab_z + g_y_y_yyy_xxz[k];

                g_y_y_yyyz_xy[k] = -g_y_y_yyy_xy[k] * ab_z + g_y_y_yyy_xyz[k];

                g_y_y_yyyz_xz[k] = -g_y_y_yyy_xz[k] * ab_z + g_y_y_yyy_xzz[k];

                g_y_y_yyyz_yy[k] = -g_y_y_yyy_yy[k] * ab_z + g_y_y_yyy_yyz[k];

                g_y_y_yyyz_yz[k] = -g_y_y_yyy_yz[k] * ab_z + g_y_y_yyy_yzz[k];

                g_y_y_yyyz_zz[k] = -g_y_y_yyy_zz[k] * ab_z + g_y_y_yyy_zzz[k];
            }

            /// Set up 432-438 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyzz_xx = cbuffer.data(gd_geom_11_off + 432 * ccomps * dcomps);

            auto g_y_y_yyzz_xy = cbuffer.data(gd_geom_11_off + 433 * ccomps * dcomps);

            auto g_y_y_yyzz_xz = cbuffer.data(gd_geom_11_off + 434 * ccomps * dcomps);

            auto g_y_y_yyzz_yy = cbuffer.data(gd_geom_11_off + 435 * ccomps * dcomps);

            auto g_y_y_yyzz_yz = cbuffer.data(gd_geom_11_off + 436 * ccomps * dcomps);

            auto g_y_y_yyzz_zz = cbuffer.data(gd_geom_11_off + 437 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yyz_xx, g_y_y_yyz_xxz, g_y_y_yyz_xy, g_y_y_yyz_xyz, g_y_y_yyz_xz, g_y_y_yyz_xzz, g_y_y_yyz_yy, g_y_y_yyz_yyz, g_y_y_yyz_yz, g_y_y_yyz_yzz, g_y_y_yyz_zz, g_y_y_yyz_zzz, g_y_y_yyzz_xx, g_y_y_yyzz_xy, g_y_y_yyzz_xz, g_y_y_yyzz_yy, g_y_y_yyzz_yz, g_y_y_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyzz_xx[k] = -g_y_y_yyz_xx[k] * ab_z + g_y_y_yyz_xxz[k];

                g_y_y_yyzz_xy[k] = -g_y_y_yyz_xy[k] * ab_z + g_y_y_yyz_xyz[k];

                g_y_y_yyzz_xz[k] = -g_y_y_yyz_xz[k] * ab_z + g_y_y_yyz_xzz[k];

                g_y_y_yyzz_yy[k] = -g_y_y_yyz_yy[k] * ab_z + g_y_y_yyz_yyz[k];

                g_y_y_yyzz_yz[k] = -g_y_y_yyz_yz[k] * ab_z + g_y_y_yyz_yzz[k];

                g_y_y_yyzz_zz[k] = -g_y_y_yyz_zz[k] * ab_z + g_y_y_yyz_zzz[k];
            }

            /// Set up 438-444 components of targeted buffer : cbuffer.data(

            auto g_y_y_yzzz_xx = cbuffer.data(gd_geom_11_off + 438 * ccomps * dcomps);

            auto g_y_y_yzzz_xy = cbuffer.data(gd_geom_11_off + 439 * ccomps * dcomps);

            auto g_y_y_yzzz_xz = cbuffer.data(gd_geom_11_off + 440 * ccomps * dcomps);

            auto g_y_y_yzzz_yy = cbuffer.data(gd_geom_11_off + 441 * ccomps * dcomps);

            auto g_y_y_yzzz_yz = cbuffer.data(gd_geom_11_off + 442 * ccomps * dcomps);

            auto g_y_y_yzzz_zz = cbuffer.data(gd_geom_11_off + 443 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yzz_xx, g_y_y_yzz_xxz, g_y_y_yzz_xy, g_y_y_yzz_xyz, g_y_y_yzz_xz, g_y_y_yzz_xzz, g_y_y_yzz_yy, g_y_y_yzz_yyz, g_y_y_yzz_yz, g_y_y_yzz_yzz, g_y_y_yzz_zz, g_y_y_yzz_zzz, g_y_y_yzzz_xx, g_y_y_yzzz_xy, g_y_y_yzzz_xz, g_y_y_yzzz_yy, g_y_y_yzzz_yz, g_y_y_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yzzz_xx[k] = -g_y_y_yzz_xx[k] * ab_z + g_y_y_yzz_xxz[k];

                g_y_y_yzzz_xy[k] = -g_y_y_yzz_xy[k] * ab_z + g_y_y_yzz_xyz[k];

                g_y_y_yzzz_xz[k] = -g_y_y_yzz_xz[k] * ab_z + g_y_y_yzz_xzz[k];

                g_y_y_yzzz_yy[k] = -g_y_y_yzz_yy[k] * ab_z + g_y_y_yzz_yyz[k];

                g_y_y_yzzz_yz[k] = -g_y_y_yzz_yz[k] * ab_z + g_y_y_yzz_yzz[k];

                g_y_y_yzzz_zz[k] = -g_y_y_yzz_zz[k] * ab_z + g_y_y_yzz_zzz[k];
            }

            /// Set up 444-450 components of targeted buffer : cbuffer.data(

            auto g_y_y_zzzz_xx = cbuffer.data(gd_geom_11_off + 444 * ccomps * dcomps);

            auto g_y_y_zzzz_xy = cbuffer.data(gd_geom_11_off + 445 * ccomps * dcomps);

            auto g_y_y_zzzz_xz = cbuffer.data(gd_geom_11_off + 446 * ccomps * dcomps);

            auto g_y_y_zzzz_yy = cbuffer.data(gd_geom_11_off + 447 * ccomps * dcomps);

            auto g_y_y_zzzz_yz = cbuffer.data(gd_geom_11_off + 448 * ccomps * dcomps);

            auto g_y_y_zzzz_zz = cbuffer.data(gd_geom_11_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_zzz_xx, g_y_y_zzz_xxz, g_y_y_zzz_xy, g_y_y_zzz_xyz, g_y_y_zzz_xz, g_y_y_zzz_xzz, g_y_y_zzz_yy, g_y_y_zzz_yyz, g_y_y_zzz_yz, g_y_y_zzz_yzz, g_y_y_zzz_zz, g_y_y_zzz_zzz, g_y_y_zzzz_xx, g_y_y_zzzz_xy, g_y_y_zzzz_xz, g_y_y_zzzz_yy, g_y_y_zzzz_yz, g_y_y_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_zzzz_xx[k] = -g_y_y_zzz_xx[k] * ab_z + g_y_y_zzz_xxz[k];

                g_y_y_zzzz_xy[k] = -g_y_y_zzz_xy[k] * ab_z + g_y_y_zzz_xyz[k];

                g_y_y_zzzz_xz[k] = -g_y_y_zzz_xz[k] * ab_z + g_y_y_zzz_xzz[k];

                g_y_y_zzzz_yy[k] = -g_y_y_zzz_yy[k] * ab_z + g_y_y_zzz_yyz[k];

                g_y_y_zzzz_yz[k] = -g_y_y_zzz_yz[k] * ab_z + g_y_y_zzz_yzz[k];

                g_y_y_zzzz_zz[k] = -g_y_y_zzz_zz[k] * ab_z + g_y_y_zzz_zzz[k];
            }

            /// Set up 450-456 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxxx_xx = cbuffer.data(gd_geom_11_off + 450 * ccomps * dcomps);

            auto g_y_z_xxxx_xy = cbuffer.data(gd_geom_11_off + 451 * ccomps * dcomps);

            auto g_y_z_xxxx_xz = cbuffer.data(gd_geom_11_off + 452 * ccomps * dcomps);

            auto g_y_z_xxxx_yy = cbuffer.data(gd_geom_11_off + 453 * ccomps * dcomps);

            auto g_y_z_xxxx_yz = cbuffer.data(gd_geom_11_off + 454 * ccomps * dcomps);

            auto g_y_z_xxxx_zz = cbuffer.data(gd_geom_11_off + 455 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxx_xx, g_y_z_xxx_xxx, g_y_z_xxx_xxy, g_y_z_xxx_xxz, g_y_z_xxx_xy, g_y_z_xxx_xyy, g_y_z_xxx_xyz, g_y_z_xxx_xz, g_y_z_xxx_xzz, g_y_z_xxx_yy, g_y_z_xxx_yz, g_y_z_xxx_zz, g_y_z_xxxx_xx, g_y_z_xxxx_xy, g_y_z_xxxx_xz, g_y_z_xxxx_yy, g_y_z_xxxx_yz, g_y_z_xxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxxx_xx[k] = -g_y_z_xxx_xx[k] * ab_x + g_y_z_xxx_xxx[k];

                g_y_z_xxxx_xy[k] = -g_y_z_xxx_xy[k] * ab_x + g_y_z_xxx_xxy[k];

                g_y_z_xxxx_xz[k] = -g_y_z_xxx_xz[k] * ab_x + g_y_z_xxx_xxz[k];

                g_y_z_xxxx_yy[k] = -g_y_z_xxx_yy[k] * ab_x + g_y_z_xxx_xyy[k];

                g_y_z_xxxx_yz[k] = -g_y_z_xxx_yz[k] * ab_x + g_y_z_xxx_xyz[k];

                g_y_z_xxxx_zz[k] = -g_y_z_xxx_zz[k] * ab_x + g_y_z_xxx_xzz[k];
            }

            /// Set up 456-462 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxxy_xx = cbuffer.data(gd_geom_11_off + 456 * ccomps * dcomps);

            auto g_y_z_xxxy_xy = cbuffer.data(gd_geom_11_off + 457 * ccomps * dcomps);

            auto g_y_z_xxxy_xz = cbuffer.data(gd_geom_11_off + 458 * ccomps * dcomps);

            auto g_y_z_xxxy_yy = cbuffer.data(gd_geom_11_off + 459 * ccomps * dcomps);

            auto g_y_z_xxxy_yz = cbuffer.data(gd_geom_11_off + 460 * ccomps * dcomps);

            auto g_y_z_xxxy_zz = cbuffer.data(gd_geom_11_off + 461 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxxy_xx, g_y_z_xxxy_xy, g_y_z_xxxy_xz, g_y_z_xxxy_yy, g_y_z_xxxy_yz, g_y_z_xxxy_zz, g_y_z_xxy_xx, g_y_z_xxy_xxx, g_y_z_xxy_xxy, g_y_z_xxy_xxz, g_y_z_xxy_xy, g_y_z_xxy_xyy, g_y_z_xxy_xyz, g_y_z_xxy_xz, g_y_z_xxy_xzz, g_y_z_xxy_yy, g_y_z_xxy_yz, g_y_z_xxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxxy_xx[k] = -g_y_z_xxy_xx[k] * ab_x + g_y_z_xxy_xxx[k];

                g_y_z_xxxy_xy[k] = -g_y_z_xxy_xy[k] * ab_x + g_y_z_xxy_xxy[k];

                g_y_z_xxxy_xz[k] = -g_y_z_xxy_xz[k] * ab_x + g_y_z_xxy_xxz[k];

                g_y_z_xxxy_yy[k] = -g_y_z_xxy_yy[k] * ab_x + g_y_z_xxy_xyy[k];

                g_y_z_xxxy_yz[k] = -g_y_z_xxy_yz[k] * ab_x + g_y_z_xxy_xyz[k];

                g_y_z_xxxy_zz[k] = -g_y_z_xxy_zz[k] * ab_x + g_y_z_xxy_xzz[k];
            }

            /// Set up 462-468 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxxz_xx = cbuffer.data(gd_geom_11_off + 462 * ccomps * dcomps);

            auto g_y_z_xxxz_xy = cbuffer.data(gd_geom_11_off + 463 * ccomps * dcomps);

            auto g_y_z_xxxz_xz = cbuffer.data(gd_geom_11_off + 464 * ccomps * dcomps);

            auto g_y_z_xxxz_yy = cbuffer.data(gd_geom_11_off + 465 * ccomps * dcomps);

            auto g_y_z_xxxz_yz = cbuffer.data(gd_geom_11_off + 466 * ccomps * dcomps);

            auto g_y_z_xxxz_zz = cbuffer.data(gd_geom_11_off + 467 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxxz_xx, g_y_z_xxxz_xy, g_y_z_xxxz_xz, g_y_z_xxxz_yy, g_y_z_xxxz_yz, g_y_z_xxxz_zz, g_y_z_xxz_xx, g_y_z_xxz_xxx, g_y_z_xxz_xxy, g_y_z_xxz_xxz, g_y_z_xxz_xy, g_y_z_xxz_xyy, g_y_z_xxz_xyz, g_y_z_xxz_xz, g_y_z_xxz_xzz, g_y_z_xxz_yy, g_y_z_xxz_yz, g_y_z_xxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxxz_xx[k] = -g_y_z_xxz_xx[k] * ab_x + g_y_z_xxz_xxx[k];

                g_y_z_xxxz_xy[k] = -g_y_z_xxz_xy[k] * ab_x + g_y_z_xxz_xxy[k];

                g_y_z_xxxz_xz[k] = -g_y_z_xxz_xz[k] * ab_x + g_y_z_xxz_xxz[k];

                g_y_z_xxxz_yy[k] = -g_y_z_xxz_yy[k] * ab_x + g_y_z_xxz_xyy[k];

                g_y_z_xxxz_yz[k] = -g_y_z_xxz_yz[k] * ab_x + g_y_z_xxz_xyz[k];

                g_y_z_xxxz_zz[k] = -g_y_z_xxz_zz[k] * ab_x + g_y_z_xxz_xzz[k];
            }

            /// Set up 468-474 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxyy_xx = cbuffer.data(gd_geom_11_off + 468 * ccomps * dcomps);

            auto g_y_z_xxyy_xy = cbuffer.data(gd_geom_11_off + 469 * ccomps * dcomps);

            auto g_y_z_xxyy_xz = cbuffer.data(gd_geom_11_off + 470 * ccomps * dcomps);

            auto g_y_z_xxyy_yy = cbuffer.data(gd_geom_11_off + 471 * ccomps * dcomps);

            auto g_y_z_xxyy_yz = cbuffer.data(gd_geom_11_off + 472 * ccomps * dcomps);

            auto g_y_z_xxyy_zz = cbuffer.data(gd_geom_11_off + 473 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxyy_xx, g_y_z_xxyy_xy, g_y_z_xxyy_xz, g_y_z_xxyy_yy, g_y_z_xxyy_yz, g_y_z_xxyy_zz, g_y_z_xyy_xx, g_y_z_xyy_xxx, g_y_z_xyy_xxy, g_y_z_xyy_xxz, g_y_z_xyy_xy, g_y_z_xyy_xyy, g_y_z_xyy_xyz, g_y_z_xyy_xz, g_y_z_xyy_xzz, g_y_z_xyy_yy, g_y_z_xyy_yz, g_y_z_xyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxyy_xx[k] = -g_y_z_xyy_xx[k] * ab_x + g_y_z_xyy_xxx[k];

                g_y_z_xxyy_xy[k] = -g_y_z_xyy_xy[k] * ab_x + g_y_z_xyy_xxy[k];

                g_y_z_xxyy_xz[k] = -g_y_z_xyy_xz[k] * ab_x + g_y_z_xyy_xxz[k];

                g_y_z_xxyy_yy[k] = -g_y_z_xyy_yy[k] * ab_x + g_y_z_xyy_xyy[k];

                g_y_z_xxyy_yz[k] = -g_y_z_xyy_yz[k] * ab_x + g_y_z_xyy_xyz[k];

                g_y_z_xxyy_zz[k] = -g_y_z_xyy_zz[k] * ab_x + g_y_z_xyy_xzz[k];
            }

            /// Set up 474-480 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxyz_xx = cbuffer.data(gd_geom_11_off + 474 * ccomps * dcomps);

            auto g_y_z_xxyz_xy = cbuffer.data(gd_geom_11_off + 475 * ccomps * dcomps);

            auto g_y_z_xxyz_xz = cbuffer.data(gd_geom_11_off + 476 * ccomps * dcomps);

            auto g_y_z_xxyz_yy = cbuffer.data(gd_geom_11_off + 477 * ccomps * dcomps);

            auto g_y_z_xxyz_yz = cbuffer.data(gd_geom_11_off + 478 * ccomps * dcomps);

            auto g_y_z_xxyz_zz = cbuffer.data(gd_geom_11_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxyz_xx, g_y_z_xxyz_xy, g_y_z_xxyz_xz, g_y_z_xxyz_yy, g_y_z_xxyz_yz, g_y_z_xxyz_zz, g_y_z_xyz_xx, g_y_z_xyz_xxx, g_y_z_xyz_xxy, g_y_z_xyz_xxz, g_y_z_xyz_xy, g_y_z_xyz_xyy, g_y_z_xyz_xyz, g_y_z_xyz_xz, g_y_z_xyz_xzz, g_y_z_xyz_yy, g_y_z_xyz_yz, g_y_z_xyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxyz_xx[k] = -g_y_z_xyz_xx[k] * ab_x + g_y_z_xyz_xxx[k];

                g_y_z_xxyz_xy[k] = -g_y_z_xyz_xy[k] * ab_x + g_y_z_xyz_xxy[k];

                g_y_z_xxyz_xz[k] = -g_y_z_xyz_xz[k] * ab_x + g_y_z_xyz_xxz[k];

                g_y_z_xxyz_yy[k] = -g_y_z_xyz_yy[k] * ab_x + g_y_z_xyz_xyy[k];

                g_y_z_xxyz_yz[k] = -g_y_z_xyz_yz[k] * ab_x + g_y_z_xyz_xyz[k];

                g_y_z_xxyz_zz[k] = -g_y_z_xyz_zz[k] * ab_x + g_y_z_xyz_xzz[k];
            }

            /// Set up 480-486 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxzz_xx = cbuffer.data(gd_geom_11_off + 480 * ccomps * dcomps);

            auto g_y_z_xxzz_xy = cbuffer.data(gd_geom_11_off + 481 * ccomps * dcomps);

            auto g_y_z_xxzz_xz = cbuffer.data(gd_geom_11_off + 482 * ccomps * dcomps);

            auto g_y_z_xxzz_yy = cbuffer.data(gd_geom_11_off + 483 * ccomps * dcomps);

            auto g_y_z_xxzz_yz = cbuffer.data(gd_geom_11_off + 484 * ccomps * dcomps);

            auto g_y_z_xxzz_zz = cbuffer.data(gd_geom_11_off + 485 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxzz_xx, g_y_z_xxzz_xy, g_y_z_xxzz_xz, g_y_z_xxzz_yy, g_y_z_xxzz_yz, g_y_z_xxzz_zz, g_y_z_xzz_xx, g_y_z_xzz_xxx, g_y_z_xzz_xxy, g_y_z_xzz_xxz, g_y_z_xzz_xy, g_y_z_xzz_xyy, g_y_z_xzz_xyz, g_y_z_xzz_xz, g_y_z_xzz_xzz, g_y_z_xzz_yy, g_y_z_xzz_yz, g_y_z_xzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxzz_xx[k] = -g_y_z_xzz_xx[k] * ab_x + g_y_z_xzz_xxx[k];

                g_y_z_xxzz_xy[k] = -g_y_z_xzz_xy[k] * ab_x + g_y_z_xzz_xxy[k];

                g_y_z_xxzz_xz[k] = -g_y_z_xzz_xz[k] * ab_x + g_y_z_xzz_xxz[k];

                g_y_z_xxzz_yy[k] = -g_y_z_xzz_yy[k] * ab_x + g_y_z_xzz_xyy[k];

                g_y_z_xxzz_yz[k] = -g_y_z_xzz_yz[k] * ab_x + g_y_z_xzz_xyz[k];

                g_y_z_xxzz_zz[k] = -g_y_z_xzz_zz[k] * ab_x + g_y_z_xzz_xzz[k];
            }

            /// Set up 486-492 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyyy_xx = cbuffer.data(gd_geom_11_off + 486 * ccomps * dcomps);

            auto g_y_z_xyyy_xy = cbuffer.data(gd_geom_11_off + 487 * ccomps * dcomps);

            auto g_y_z_xyyy_xz = cbuffer.data(gd_geom_11_off + 488 * ccomps * dcomps);

            auto g_y_z_xyyy_yy = cbuffer.data(gd_geom_11_off + 489 * ccomps * dcomps);

            auto g_y_z_xyyy_yz = cbuffer.data(gd_geom_11_off + 490 * ccomps * dcomps);

            auto g_y_z_xyyy_zz = cbuffer.data(gd_geom_11_off + 491 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyyy_xx, g_y_z_xyyy_xy, g_y_z_xyyy_xz, g_y_z_xyyy_yy, g_y_z_xyyy_yz, g_y_z_xyyy_zz, g_y_z_yyy_xx, g_y_z_yyy_xxx, g_y_z_yyy_xxy, g_y_z_yyy_xxz, g_y_z_yyy_xy, g_y_z_yyy_xyy, g_y_z_yyy_xyz, g_y_z_yyy_xz, g_y_z_yyy_xzz, g_y_z_yyy_yy, g_y_z_yyy_yz, g_y_z_yyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyyy_xx[k] = -g_y_z_yyy_xx[k] * ab_x + g_y_z_yyy_xxx[k];

                g_y_z_xyyy_xy[k] = -g_y_z_yyy_xy[k] * ab_x + g_y_z_yyy_xxy[k];

                g_y_z_xyyy_xz[k] = -g_y_z_yyy_xz[k] * ab_x + g_y_z_yyy_xxz[k];

                g_y_z_xyyy_yy[k] = -g_y_z_yyy_yy[k] * ab_x + g_y_z_yyy_xyy[k];

                g_y_z_xyyy_yz[k] = -g_y_z_yyy_yz[k] * ab_x + g_y_z_yyy_xyz[k];

                g_y_z_xyyy_zz[k] = -g_y_z_yyy_zz[k] * ab_x + g_y_z_yyy_xzz[k];
            }

            /// Set up 492-498 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyyz_xx = cbuffer.data(gd_geom_11_off + 492 * ccomps * dcomps);

            auto g_y_z_xyyz_xy = cbuffer.data(gd_geom_11_off + 493 * ccomps * dcomps);

            auto g_y_z_xyyz_xz = cbuffer.data(gd_geom_11_off + 494 * ccomps * dcomps);

            auto g_y_z_xyyz_yy = cbuffer.data(gd_geom_11_off + 495 * ccomps * dcomps);

            auto g_y_z_xyyz_yz = cbuffer.data(gd_geom_11_off + 496 * ccomps * dcomps);

            auto g_y_z_xyyz_zz = cbuffer.data(gd_geom_11_off + 497 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyyz_xx, g_y_z_xyyz_xy, g_y_z_xyyz_xz, g_y_z_xyyz_yy, g_y_z_xyyz_yz, g_y_z_xyyz_zz, g_y_z_yyz_xx, g_y_z_yyz_xxx, g_y_z_yyz_xxy, g_y_z_yyz_xxz, g_y_z_yyz_xy, g_y_z_yyz_xyy, g_y_z_yyz_xyz, g_y_z_yyz_xz, g_y_z_yyz_xzz, g_y_z_yyz_yy, g_y_z_yyz_yz, g_y_z_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyyz_xx[k] = -g_y_z_yyz_xx[k] * ab_x + g_y_z_yyz_xxx[k];

                g_y_z_xyyz_xy[k] = -g_y_z_yyz_xy[k] * ab_x + g_y_z_yyz_xxy[k];

                g_y_z_xyyz_xz[k] = -g_y_z_yyz_xz[k] * ab_x + g_y_z_yyz_xxz[k];

                g_y_z_xyyz_yy[k] = -g_y_z_yyz_yy[k] * ab_x + g_y_z_yyz_xyy[k];

                g_y_z_xyyz_yz[k] = -g_y_z_yyz_yz[k] * ab_x + g_y_z_yyz_xyz[k];

                g_y_z_xyyz_zz[k] = -g_y_z_yyz_zz[k] * ab_x + g_y_z_yyz_xzz[k];
            }

            /// Set up 498-504 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyzz_xx = cbuffer.data(gd_geom_11_off + 498 * ccomps * dcomps);

            auto g_y_z_xyzz_xy = cbuffer.data(gd_geom_11_off + 499 * ccomps * dcomps);

            auto g_y_z_xyzz_xz = cbuffer.data(gd_geom_11_off + 500 * ccomps * dcomps);

            auto g_y_z_xyzz_yy = cbuffer.data(gd_geom_11_off + 501 * ccomps * dcomps);

            auto g_y_z_xyzz_yz = cbuffer.data(gd_geom_11_off + 502 * ccomps * dcomps);

            auto g_y_z_xyzz_zz = cbuffer.data(gd_geom_11_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyzz_xx, g_y_z_xyzz_xy, g_y_z_xyzz_xz, g_y_z_xyzz_yy, g_y_z_xyzz_yz, g_y_z_xyzz_zz, g_y_z_yzz_xx, g_y_z_yzz_xxx, g_y_z_yzz_xxy, g_y_z_yzz_xxz, g_y_z_yzz_xy, g_y_z_yzz_xyy, g_y_z_yzz_xyz, g_y_z_yzz_xz, g_y_z_yzz_xzz, g_y_z_yzz_yy, g_y_z_yzz_yz, g_y_z_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyzz_xx[k] = -g_y_z_yzz_xx[k] * ab_x + g_y_z_yzz_xxx[k];

                g_y_z_xyzz_xy[k] = -g_y_z_yzz_xy[k] * ab_x + g_y_z_yzz_xxy[k];

                g_y_z_xyzz_xz[k] = -g_y_z_yzz_xz[k] * ab_x + g_y_z_yzz_xxz[k];

                g_y_z_xyzz_yy[k] = -g_y_z_yzz_yy[k] * ab_x + g_y_z_yzz_xyy[k];

                g_y_z_xyzz_yz[k] = -g_y_z_yzz_yz[k] * ab_x + g_y_z_yzz_xyz[k];

                g_y_z_xyzz_zz[k] = -g_y_z_yzz_zz[k] * ab_x + g_y_z_yzz_xzz[k];
            }

            /// Set up 504-510 components of targeted buffer : cbuffer.data(

            auto g_y_z_xzzz_xx = cbuffer.data(gd_geom_11_off + 504 * ccomps * dcomps);

            auto g_y_z_xzzz_xy = cbuffer.data(gd_geom_11_off + 505 * ccomps * dcomps);

            auto g_y_z_xzzz_xz = cbuffer.data(gd_geom_11_off + 506 * ccomps * dcomps);

            auto g_y_z_xzzz_yy = cbuffer.data(gd_geom_11_off + 507 * ccomps * dcomps);

            auto g_y_z_xzzz_yz = cbuffer.data(gd_geom_11_off + 508 * ccomps * dcomps);

            auto g_y_z_xzzz_zz = cbuffer.data(gd_geom_11_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xzzz_xx, g_y_z_xzzz_xy, g_y_z_xzzz_xz, g_y_z_xzzz_yy, g_y_z_xzzz_yz, g_y_z_xzzz_zz, g_y_z_zzz_xx, g_y_z_zzz_xxx, g_y_z_zzz_xxy, g_y_z_zzz_xxz, g_y_z_zzz_xy, g_y_z_zzz_xyy, g_y_z_zzz_xyz, g_y_z_zzz_xz, g_y_z_zzz_xzz, g_y_z_zzz_yy, g_y_z_zzz_yz, g_y_z_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xzzz_xx[k] = -g_y_z_zzz_xx[k] * ab_x + g_y_z_zzz_xxx[k];

                g_y_z_xzzz_xy[k] = -g_y_z_zzz_xy[k] * ab_x + g_y_z_zzz_xxy[k];

                g_y_z_xzzz_xz[k] = -g_y_z_zzz_xz[k] * ab_x + g_y_z_zzz_xxz[k];

                g_y_z_xzzz_yy[k] = -g_y_z_zzz_yy[k] * ab_x + g_y_z_zzz_xyy[k];

                g_y_z_xzzz_yz[k] = -g_y_z_zzz_yz[k] * ab_x + g_y_z_zzz_xyz[k];

                g_y_z_xzzz_zz[k] = -g_y_z_zzz_zz[k] * ab_x + g_y_z_zzz_xzz[k];
            }

            /// Set up 510-516 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyyy_xx = cbuffer.data(gd_geom_11_off + 510 * ccomps * dcomps);

            auto g_y_z_yyyy_xy = cbuffer.data(gd_geom_11_off + 511 * ccomps * dcomps);

            auto g_y_z_yyyy_xz = cbuffer.data(gd_geom_11_off + 512 * ccomps * dcomps);

            auto g_y_z_yyyy_yy = cbuffer.data(gd_geom_11_off + 513 * ccomps * dcomps);

            auto g_y_z_yyyy_yz = cbuffer.data(gd_geom_11_off + 514 * ccomps * dcomps);

            auto g_y_z_yyyy_zz = cbuffer.data(gd_geom_11_off + 515 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyy_xx, g_0_z_yyy_xy, g_0_z_yyy_xz, g_0_z_yyy_yy, g_0_z_yyy_yz, g_0_z_yyy_zz, g_y_z_yyy_xx, g_y_z_yyy_xxy, g_y_z_yyy_xy, g_y_z_yyy_xyy, g_y_z_yyy_xyz, g_y_z_yyy_xz, g_y_z_yyy_yy, g_y_z_yyy_yyy, g_y_z_yyy_yyz, g_y_z_yyy_yz, g_y_z_yyy_yzz, g_y_z_yyy_zz, g_y_z_yyyy_xx, g_y_z_yyyy_xy, g_y_z_yyyy_xz, g_y_z_yyyy_yy, g_y_z_yyyy_yz, g_y_z_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyyy_xx[k] = -g_0_z_yyy_xx[k] - g_y_z_yyy_xx[k] * ab_y + g_y_z_yyy_xxy[k];

                g_y_z_yyyy_xy[k] = -g_0_z_yyy_xy[k] - g_y_z_yyy_xy[k] * ab_y + g_y_z_yyy_xyy[k];

                g_y_z_yyyy_xz[k] = -g_0_z_yyy_xz[k] - g_y_z_yyy_xz[k] * ab_y + g_y_z_yyy_xyz[k];

                g_y_z_yyyy_yy[k] = -g_0_z_yyy_yy[k] - g_y_z_yyy_yy[k] * ab_y + g_y_z_yyy_yyy[k];

                g_y_z_yyyy_yz[k] = -g_0_z_yyy_yz[k] - g_y_z_yyy_yz[k] * ab_y + g_y_z_yyy_yyz[k];

                g_y_z_yyyy_zz[k] = -g_0_z_yyy_zz[k] - g_y_z_yyy_zz[k] * ab_y + g_y_z_yyy_yzz[k];
            }

            /// Set up 516-522 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyyz_xx = cbuffer.data(gd_geom_11_off + 516 * ccomps * dcomps);

            auto g_y_z_yyyz_xy = cbuffer.data(gd_geom_11_off + 517 * ccomps * dcomps);

            auto g_y_z_yyyz_xz = cbuffer.data(gd_geom_11_off + 518 * ccomps * dcomps);

            auto g_y_z_yyyz_yy = cbuffer.data(gd_geom_11_off + 519 * ccomps * dcomps);

            auto g_y_z_yyyz_yz = cbuffer.data(gd_geom_11_off + 520 * ccomps * dcomps);

            auto g_y_z_yyyz_zz = cbuffer.data(gd_geom_11_off + 521 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyz_xx, g_0_z_yyz_xy, g_0_z_yyz_xz, g_0_z_yyz_yy, g_0_z_yyz_yz, g_0_z_yyz_zz, g_y_z_yyyz_xx, g_y_z_yyyz_xy, g_y_z_yyyz_xz, g_y_z_yyyz_yy, g_y_z_yyyz_yz, g_y_z_yyyz_zz, g_y_z_yyz_xx, g_y_z_yyz_xxy, g_y_z_yyz_xy, g_y_z_yyz_xyy, g_y_z_yyz_xyz, g_y_z_yyz_xz, g_y_z_yyz_yy, g_y_z_yyz_yyy, g_y_z_yyz_yyz, g_y_z_yyz_yz, g_y_z_yyz_yzz, g_y_z_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyyz_xx[k] = -g_0_z_yyz_xx[k] - g_y_z_yyz_xx[k] * ab_y + g_y_z_yyz_xxy[k];

                g_y_z_yyyz_xy[k] = -g_0_z_yyz_xy[k] - g_y_z_yyz_xy[k] * ab_y + g_y_z_yyz_xyy[k];

                g_y_z_yyyz_xz[k] = -g_0_z_yyz_xz[k] - g_y_z_yyz_xz[k] * ab_y + g_y_z_yyz_xyz[k];

                g_y_z_yyyz_yy[k] = -g_0_z_yyz_yy[k] - g_y_z_yyz_yy[k] * ab_y + g_y_z_yyz_yyy[k];

                g_y_z_yyyz_yz[k] = -g_0_z_yyz_yz[k] - g_y_z_yyz_yz[k] * ab_y + g_y_z_yyz_yyz[k];

                g_y_z_yyyz_zz[k] = -g_0_z_yyz_zz[k] - g_y_z_yyz_zz[k] * ab_y + g_y_z_yyz_yzz[k];
            }

            /// Set up 522-528 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyzz_xx = cbuffer.data(gd_geom_11_off + 522 * ccomps * dcomps);

            auto g_y_z_yyzz_xy = cbuffer.data(gd_geom_11_off + 523 * ccomps * dcomps);

            auto g_y_z_yyzz_xz = cbuffer.data(gd_geom_11_off + 524 * ccomps * dcomps);

            auto g_y_z_yyzz_yy = cbuffer.data(gd_geom_11_off + 525 * ccomps * dcomps);

            auto g_y_z_yyzz_yz = cbuffer.data(gd_geom_11_off + 526 * ccomps * dcomps);

            auto g_y_z_yyzz_zz = cbuffer.data(gd_geom_11_off + 527 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzz_xx, g_0_z_yzz_xy, g_0_z_yzz_xz, g_0_z_yzz_yy, g_0_z_yzz_yz, g_0_z_yzz_zz, g_y_z_yyzz_xx, g_y_z_yyzz_xy, g_y_z_yyzz_xz, g_y_z_yyzz_yy, g_y_z_yyzz_yz, g_y_z_yyzz_zz, g_y_z_yzz_xx, g_y_z_yzz_xxy, g_y_z_yzz_xy, g_y_z_yzz_xyy, g_y_z_yzz_xyz, g_y_z_yzz_xz, g_y_z_yzz_yy, g_y_z_yzz_yyy, g_y_z_yzz_yyz, g_y_z_yzz_yz, g_y_z_yzz_yzz, g_y_z_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyzz_xx[k] = -g_0_z_yzz_xx[k] - g_y_z_yzz_xx[k] * ab_y + g_y_z_yzz_xxy[k];

                g_y_z_yyzz_xy[k] = -g_0_z_yzz_xy[k] - g_y_z_yzz_xy[k] * ab_y + g_y_z_yzz_xyy[k];

                g_y_z_yyzz_xz[k] = -g_0_z_yzz_xz[k] - g_y_z_yzz_xz[k] * ab_y + g_y_z_yzz_xyz[k];

                g_y_z_yyzz_yy[k] = -g_0_z_yzz_yy[k] - g_y_z_yzz_yy[k] * ab_y + g_y_z_yzz_yyy[k];

                g_y_z_yyzz_yz[k] = -g_0_z_yzz_yz[k] - g_y_z_yzz_yz[k] * ab_y + g_y_z_yzz_yyz[k];

                g_y_z_yyzz_zz[k] = -g_0_z_yzz_zz[k] - g_y_z_yzz_zz[k] * ab_y + g_y_z_yzz_yzz[k];
            }

            /// Set up 528-534 components of targeted buffer : cbuffer.data(

            auto g_y_z_yzzz_xx = cbuffer.data(gd_geom_11_off + 528 * ccomps * dcomps);

            auto g_y_z_yzzz_xy = cbuffer.data(gd_geom_11_off + 529 * ccomps * dcomps);

            auto g_y_z_yzzz_xz = cbuffer.data(gd_geom_11_off + 530 * ccomps * dcomps);

            auto g_y_z_yzzz_yy = cbuffer.data(gd_geom_11_off + 531 * ccomps * dcomps);

            auto g_y_z_yzzz_yz = cbuffer.data(gd_geom_11_off + 532 * ccomps * dcomps);

            auto g_y_z_yzzz_zz = cbuffer.data(gd_geom_11_off + 533 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzz_xx, g_0_z_zzz_xy, g_0_z_zzz_xz, g_0_z_zzz_yy, g_0_z_zzz_yz, g_0_z_zzz_zz, g_y_z_yzzz_xx, g_y_z_yzzz_xy, g_y_z_yzzz_xz, g_y_z_yzzz_yy, g_y_z_yzzz_yz, g_y_z_yzzz_zz, g_y_z_zzz_xx, g_y_z_zzz_xxy, g_y_z_zzz_xy, g_y_z_zzz_xyy, g_y_z_zzz_xyz, g_y_z_zzz_xz, g_y_z_zzz_yy, g_y_z_zzz_yyy, g_y_z_zzz_yyz, g_y_z_zzz_yz, g_y_z_zzz_yzz, g_y_z_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yzzz_xx[k] = -g_0_z_zzz_xx[k] - g_y_z_zzz_xx[k] * ab_y + g_y_z_zzz_xxy[k];

                g_y_z_yzzz_xy[k] = -g_0_z_zzz_xy[k] - g_y_z_zzz_xy[k] * ab_y + g_y_z_zzz_xyy[k];

                g_y_z_yzzz_xz[k] = -g_0_z_zzz_xz[k] - g_y_z_zzz_xz[k] * ab_y + g_y_z_zzz_xyz[k];

                g_y_z_yzzz_yy[k] = -g_0_z_zzz_yy[k] - g_y_z_zzz_yy[k] * ab_y + g_y_z_zzz_yyy[k];

                g_y_z_yzzz_yz[k] = -g_0_z_zzz_yz[k] - g_y_z_zzz_yz[k] * ab_y + g_y_z_zzz_yyz[k];

                g_y_z_yzzz_zz[k] = -g_0_z_zzz_zz[k] - g_y_z_zzz_zz[k] * ab_y + g_y_z_zzz_yzz[k];
            }

            /// Set up 534-540 components of targeted buffer : cbuffer.data(

            auto g_y_z_zzzz_xx = cbuffer.data(gd_geom_11_off + 534 * ccomps * dcomps);

            auto g_y_z_zzzz_xy = cbuffer.data(gd_geom_11_off + 535 * ccomps * dcomps);

            auto g_y_z_zzzz_xz = cbuffer.data(gd_geom_11_off + 536 * ccomps * dcomps);

            auto g_y_z_zzzz_yy = cbuffer.data(gd_geom_11_off + 537 * ccomps * dcomps);

            auto g_y_z_zzzz_yz = cbuffer.data(gd_geom_11_off + 538 * ccomps * dcomps);

            auto g_y_z_zzzz_zz = cbuffer.data(gd_geom_11_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzz_xx, g_y_0_zzz_xy, g_y_0_zzz_xz, g_y_0_zzz_yy, g_y_0_zzz_yz, g_y_0_zzz_zz, g_y_z_zzz_xx, g_y_z_zzz_xxz, g_y_z_zzz_xy, g_y_z_zzz_xyz, g_y_z_zzz_xz, g_y_z_zzz_xzz, g_y_z_zzz_yy, g_y_z_zzz_yyz, g_y_z_zzz_yz, g_y_z_zzz_yzz, g_y_z_zzz_zz, g_y_z_zzz_zzz, g_y_z_zzzz_xx, g_y_z_zzzz_xy, g_y_z_zzzz_xz, g_y_z_zzzz_yy, g_y_z_zzzz_yz, g_y_z_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_zzzz_xx[k] = g_y_0_zzz_xx[k] - g_y_z_zzz_xx[k] * ab_z + g_y_z_zzz_xxz[k];

                g_y_z_zzzz_xy[k] = g_y_0_zzz_xy[k] - g_y_z_zzz_xy[k] * ab_z + g_y_z_zzz_xyz[k];

                g_y_z_zzzz_xz[k] = g_y_0_zzz_xz[k] - g_y_z_zzz_xz[k] * ab_z + g_y_z_zzz_xzz[k];

                g_y_z_zzzz_yy[k] = g_y_0_zzz_yy[k] - g_y_z_zzz_yy[k] * ab_z + g_y_z_zzz_yyz[k];

                g_y_z_zzzz_yz[k] = g_y_0_zzz_yz[k] - g_y_z_zzz_yz[k] * ab_z + g_y_z_zzz_yzz[k];

                g_y_z_zzzz_zz[k] = g_y_0_zzz_zz[k] - g_y_z_zzz_zz[k] * ab_z + g_y_z_zzz_zzz[k];
            }

            /// Set up 540-546 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxxx_xx = cbuffer.data(gd_geom_11_off + 540 * ccomps * dcomps);

            auto g_z_x_xxxx_xy = cbuffer.data(gd_geom_11_off + 541 * ccomps * dcomps);

            auto g_z_x_xxxx_xz = cbuffer.data(gd_geom_11_off + 542 * ccomps * dcomps);

            auto g_z_x_xxxx_yy = cbuffer.data(gd_geom_11_off + 543 * ccomps * dcomps);

            auto g_z_x_xxxx_yz = cbuffer.data(gd_geom_11_off + 544 * ccomps * dcomps);

            auto g_z_x_xxxx_zz = cbuffer.data(gd_geom_11_off + 545 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxx_xx, g_z_0_xxx_xy, g_z_0_xxx_xz, g_z_0_xxx_yy, g_z_0_xxx_yz, g_z_0_xxx_zz, g_z_x_xxx_xx, g_z_x_xxx_xxx, g_z_x_xxx_xxy, g_z_x_xxx_xxz, g_z_x_xxx_xy, g_z_x_xxx_xyy, g_z_x_xxx_xyz, g_z_x_xxx_xz, g_z_x_xxx_xzz, g_z_x_xxx_yy, g_z_x_xxx_yz, g_z_x_xxx_zz, g_z_x_xxxx_xx, g_z_x_xxxx_xy, g_z_x_xxxx_xz, g_z_x_xxxx_yy, g_z_x_xxxx_yz, g_z_x_xxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxxx_xx[k] = g_z_0_xxx_xx[k] - g_z_x_xxx_xx[k] * ab_x + g_z_x_xxx_xxx[k];

                g_z_x_xxxx_xy[k] = g_z_0_xxx_xy[k] - g_z_x_xxx_xy[k] * ab_x + g_z_x_xxx_xxy[k];

                g_z_x_xxxx_xz[k] = g_z_0_xxx_xz[k] - g_z_x_xxx_xz[k] * ab_x + g_z_x_xxx_xxz[k];

                g_z_x_xxxx_yy[k] = g_z_0_xxx_yy[k] - g_z_x_xxx_yy[k] * ab_x + g_z_x_xxx_xyy[k];

                g_z_x_xxxx_yz[k] = g_z_0_xxx_yz[k] - g_z_x_xxx_yz[k] * ab_x + g_z_x_xxx_xyz[k];

                g_z_x_xxxx_zz[k] = g_z_0_xxx_zz[k] - g_z_x_xxx_zz[k] * ab_x + g_z_x_xxx_xzz[k];
            }

            /// Set up 546-552 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxxy_xx = cbuffer.data(gd_geom_11_off + 546 * ccomps * dcomps);

            auto g_z_x_xxxy_xy = cbuffer.data(gd_geom_11_off + 547 * ccomps * dcomps);

            auto g_z_x_xxxy_xz = cbuffer.data(gd_geom_11_off + 548 * ccomps * dcomps);

            auto g_z_x_xxxy_yy = cbuffer.data(gd_geom_11_off + 549 * ccomps * dcomps);

            auto g_z_x_xxxy_yz = cbuffer.data(gd_geom_11_off + 550 * ccomps * dcomps);

            auto g_z_x_xxxy_zz = cbuffer.data(gd_geom_11_off + 551 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xxx_xx, g_z_x_xxx_xxy, g_z_x_xxx_xy, g_z_x_xxx_xyy, g_z_x_xxx_xyz, g_z_x_xxx_xz, g_z_x_xxx_yy, g_z_x_xxx_yyy, g_z_x_xxx_yyz, g_z_x_xxx_yz, g_z_x_xxx_yzz, g_z_x_xxx_zz, g_z_x_xxxy_xx, g_z_x_xxxy_xy, g_z_x_xxxy_xz, g_z_x_xxxy_yy, g_z_x_xxxy_yz, g_z_x_xxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxxy_xx[k] = -g_z_x_xxx_xx[k] * ab_y + g_z_x_xxx_xxy[k];

                g_z_x_xxxy_xy[k] = -g_z_x_xxx_xy[k] * ab_y + g_z_x_xxx_xyy[k];

                g_z_x_xxxy_xz[k] = -g_z_x_xxx_xz[k] * ab_y + g_z_x_xxx_xyz[k];

                g_z_x_xxxy_yy[k] = -g_z_x_xxx_yy[k] * ab_y + g_z_x_xxx_yyy[k];

                g_z_x_xxxy_yz[k] = -g_z_x_xxx_yz[k] * ab_y + g_z_x_xxx_yyz[k];

                g_z_x_xxxy_zz[k] = -g_z_x_xxx_zz[k] * ab_y + g_z_x_xxx_yzz[k];
            }

            /// Set up 552-558 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxxz_xx = cbuffer.data(gd_geom_11_off + 552 * ccomps * dcomps);

            auto g_z_x_xxxz_xy = cbuffer.data(gd_geom_11_off + 553 * ccomps * dcomps);

            auto g_z_x_xxxz_xz = cbuffer.data(gd_geom_11_off + 554 * ccomps * dcomps);

            auto g_z_x_xxxz_yy = cbuffer.data(gd_geom_11_off + 555 * ccomps * dcomps);

            auto g_z_x_xxxz_yz = cbuffer.data(gd_geom_11_off + 556 * ccomps * dcomps);

            auto g_z_x_xxxz_zz = cbuffer.data(gd_geom_11_off + 557 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxz_xx, g_z_0_xxz_xy, g_z_0_xxz_xz, g_z_0_xxz_yy, g_z_0_xxz_yz, g_z_0_xxz_zz, g_z_x_xxxz_xx, g_z_x_xxxz_xy, g_z_x_xxxz_xz, g_z_x_xxxz_yy, g_z_x_xxxz_yz, g_z_x_xxxz_zz, g_z_x_xxz_xx, g_z_x_xxz_xxx, g_z_x_xxz_xxy, g_z_x_xxz_xxz, g_z_x_xxz_xy, g_z_x_xxz_xyy, g_z_x_xxz_xyz, g_z_x_xxz_xz, g_z_x_xxz_xzz, g_z_x_xxz_yy, g_z_x_xxz_yz, g_z_x_xxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxxz_xx[k] = g_z_0_xxz_xx[k] - g_z_x_xxz_xx[k] * ab_x + g_z_x_xxz_xxx[k];

                g_z_x_xxxz_xy[k] = g_z_0_xxz_xy[k] - g_z_x_xxz_xy[k] * ab_x + g_z_x_xxz_xxy[k];

                g_z_x_xxxz_xz[k] = g_z_0_xxz_xz[k] - g_z_x_xxz_xz[k] * ab_x + g_z_x_xxz_xxz[k];

                g_z_x_xxxz_yy[k] = g_z_0_xxz_yy[k] - g_z_x_xxz_yy[k] * ab_x + g_z_x_xxz_xyy[k];

                g_z_x_xxxz_yz[k] = g_z_0_xxz_yz[k] - g_z_x_xxz_yz[k] * ab_x + g_z_x_xxz_xyz[k];

                g_z_x_xxxz_zz[k] = g_z_0_xxz_zz[k] - g_z_x_xxz_zz[k] * ab_x + g_z_x_xxz_xzz[k];
            }

            /// Set up 558-564 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxyy_xx = cbuffer.data(gd_geom_11_off + 558 * ccomps * dcomps);

            auto g_z_x_xxyy_xy = cbuffer.data(gd_geom_11_off + 559 * ccomps * dcomps);

            auto g_z_x_xxyy_xz = cbuffer.data(gd_geom_11_off + 560 * ccomps * dcomps);

            auto g_z_x_xxyy_yy = cbuffer.data(gd_geom_11_off + 561 * ccomps * dcomps);

            auto g_z_x_xxyy_yz = cbuffer.data(gd_geom_11_off + 562 * ccomps * dcomps);

            auto g_z_x_xxyy_zz = cbuffer.data(gd_geom_11_off + 563 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xxy_xx, g_z_x_xxy_xxy, g_z_x_xxy_xy, g_z_x_xxy_xyy, g_z_x_xxy_xyz, g_z_x_xxy_xz, g_z_x_xxy_yy, g_z_x_xxy_yyy, g_z_x_xxy_yyz, g_z_x_xxy_yz, g_z_x_xxy_yzz, g_z_x_xxy_zz, g_z_x_xxyy_xx, g_z_x_xxyy_xy, g_z_x_xxyy_xz, g_z_x_xxyy_yy, g_z_x_xxyy_yz, g_z_x_xxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxyy_xx[k] = -g_z_x_xxy_xx[k] * ab_y + g_z_x_xxy_xxy[k];

                g_z_x_xxyy_xy[k] = -g_z_x_xxy_xy[k] * ab_y + g_z_x_xxy_xyy[k];

                g_z_x_xxyy_xz[k] = -g_z_x_xxy_xz[k] * ab_y + g_z_x_xxy_xyz[k];

                g_z_x_xxyy_yy[k] = -g_z_x_xxy_yy[k] * ab_y + g_z_x_xxy_yyy[k];

                g_z_x_xxyy_yz[k] = -g_z_x_xxy_yz[k] * ab_y + g_z_x_xxy_yyz[k];

                g_z_x_xxyy_zz[k] = -g_z_x_xxy_zz[k] * ab_y + g_z_x_xxy_yzz[k];
            }

            /// Set up 564-570 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxyz_xx = cbuffer.data(gd_geom_11_off + 564 * ccomps * dcomps);

            auto g_z_x_xxyz_xy = cbuffer.data(gd_geom_11_off + 565 * ccomps * dcomps);

            auto g_z_x_xxyz_xz = cbuffer.data(gd_geom_11_off + 566 * ccomps * dcomps);

            auto g_z_x_xxyz_yy = cbuffer.data(gd_geom_11_off + 567 * ccomps * dcomps);

            auto g_z_x_xxyz_yz = cbuffer.data(gd_geom_11_off + 568 * ccomps * dcomps);

            auto g_z_x_xxyz_zz = cbuffer.data(gd_geom_11_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xxyz_xx, g_z_x_xxyz_xy, g_z_x_xxyz_xz, g_z_x_xxyz_yy, g_z_x_xxyz_yz, g_z_x_xxyz_zz, g_z_x_xxz_xx, g_z_x_xxz_xxy, g_z_x_xxz_xy, g_z_x_xxz_xyy, g_z_x_xxz_xyz, g_z_x_xxz_xz, g_z_x_xxz_yy, g_z_x_xxz_yyy, g_z_x_xxz_yyz, g_z_x_xxz_yz, g_z_x_xxz_yzz, g_z_x_xxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxyz_xx[k] = -g_z_x_xxz_xx[k] * ab_y + g_z_x_xxz_xxy[k];

                g_z_x_xxyz_xy[k] = -g_z_x_xxz_xy[k] * ab_y + g_z_x_xxz_xyy[k];

                g_z_x_xxyz_xz[k] = -g_z_x_xxz_xz[k] * ab_y + g_z_x_xxz_xyz[k];

                g_z_x_xxyz_yy[k] = -g_z_x_xxz_yy[k] * ab_y + g_z_x_xxz_yyy[k];

                g_z_x_xxyz_yz[k] = -g_z_x_xxz_yz[k] * ab_y + g_z_x_xxz_yyz[k];

                g_z_x_xxyz_zz[k] = -g_z_x_xxz_zz[k] * ab_y + g_z_x_xxz_yzz[k];
            }

            /// Set up 570-576 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxzz_xx = cbuffer.data(gd_geom_11_off + 570 * ccomps * dcomps);

            auto g_z_x_xxzz_xy = cbuffer.data(gd_geom_11_off + 571 * ccomps * dcomps);

            auto g_z_x_xxzz_xz = cbuffer.data(gd_geom_11_off + 572 * ccomps * dcomps);

            auto g_z_x_xxzz_yy = cbuffer.data(gd_geom_11_off + 573 * ccomps * dcomps);

            auto g_z_x_xxzz_yz = cbuffer.data(gd_geom_11_off + 574 * ccomps * dcomps);

            auto g_z_x_xxzz_zz = cbuffer.data(gd_geom_11_off + 575 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xzz_xx, g_z_0_xzz_xy, g_z_0_xzz_xz, g_z_0_xzz_yy, g_z_0_xzz_yz, g_z_0_xzz_zz, g_z_x_xxzz_xx, g_z_x_xxzz_xy, g_z_x_xxzz_xz, g_z_x_xxzz_yy, g_z_x_xxzz_yz, g_z_x_xxzz_zz, g_z_x_xzz_xx, g_z_x_xzz_xxx, g_z_x_xzz_xxy, g_z_x_xzz_xxz, g_z_x_xzz_xy, g_z_x_xzz_xyy, g_z_x_xzz_xyz, g_z_x_xzz_xz, g_z_x_xzz_xzz, g_z_x_xzz_yy, g_z_x_xzz_yz, g_z_x_xzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxzz_xx[k] = g_z_0_xzz_xx[k] - g_z_x_xzz_xx[k] * ab_x + g_z_x_xzz_xxx[k];

                g_z_x_xxzz_xy[k] = g_z_0_xzz_xy[k] - g_z_x_xzz_xy[k] * ab_x + g_z_x_xzz_xxy[k];

                g_z_x_xxzz_xz[k] = g_z_0_xzz_xz[k] - g_z_x_xzz_xz[k] * ab_x + g_z_x_xzz_xxz[k];

                g_z_x_xxzz_yy[k] = g_z_0_xzz_yy[k] - g_z_x_xzz_yy[k] * ab_x + g_z_x_xzz_xyy[k];

                g_z_x_xxzz_yz[k] = g_z_0_xzz_yz[k] - g_z_x_xzz_yz[k] * ab_x + g_z_x_xzz_xyz[k];

                g_z_x_xxzz_zz[k] = g_z_0_xzz_zz[k] - g_z_x_xzz_zz[k] * ab_x + g_z_x_xzz_xzz[k];
            }

            /// Set up 576-582 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyyy_xx = cbuffer.data(gd_geom_11_off + 576 * ccomps * dcomps);

            auto g_z_x_xyyy_xy = cbuffer.data(gd_geom_11_off + 577 * ccomps * dcomps);

            auto g_z_x_xyyy_xz = cbuffer.data(gd_geom_11_off + 578 * ccomps * dcomps);

            auto g_z_x_xyyy_yy = cbuffer.data(gd_geom_11_off + 579 * ccomps * dcomps);

            auto g_z_x_xyyy_yz = cbuffer.data(gd_geom_11_off + 580 * ccomps * dcomps);

            auto g_z_x_xyyy_zz = cbuffer.data(gd_geom_11_off + 581 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xyy_xx, g_z_x_xyy_xxy, g_z_x_xyy_xy, g_z_x_xyy_xyy, g_z_x_xyy_xyz, g_z_x_xyy_xz, g_z_x_xyy_yy, g_z_x_xyy_yyy, g_z_x_xyy_yyz, g_z_x_xyy_yz, g_z_x_xyy_yzz, g_z_x_xyy_zz, g_z_x_xyyy_xx, g_z_x_xyyy_xy, g_z_x_xyyy_xz, g_z_x_xyyy_yy, g_z_x_xyyy_yz, g_z_x_xyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyyy_xx[k] = -g_z_x_xyy_xx[k] * ab_y + g_z_x_xyy_xxy[k];

                g_z_x_xyyy_xy[k] = -g_z_x_xyy_xy[k] * ab_y + g_z_x_xyy_xyy[k];

                g_z_x_xyyy_xz[k] = -g_z_x_xyy_xz[k] * ab_y + g_z_x_xyy_xyz[k];

                g_z_x_xyyy_yy[k] = -g_z_x_xyy_yy[k] * ab_y + g_z_x_xyy_yyy[k];

                g_z_x_xyyy_yz[k] = -g_z_x_xyy_yz[k] * ab_y + g_z_x_xyy_yyz[k];

                g_z_x_xyyy_zz[k] = -g_z_x_xyy_zz[k] * ab_y + g_z_x_xyy_yzz[k];
            }

            /// Set up 582-588 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyyz_xx = cbuffer.data(gd_geom_11_off + 582 * ccomps * dcomps);

            auto g_z_x_xyyz_xy = cbuffer.data(gd_geom_11_off + 583 * ccomps * dcomps);

            auto g_z_x_xyyz_xz = cbuffer.data(gd_geom_11_off + 584 * ccomps * dcomps);

            auto g_z_x_xyyz_yy = cbuffer.data(gd_geom_11_off + 585 * ccomps * dcomps);

            auto g_z_x_xyyz_yz = cbuffer.data(gd_geom_11_off + 586 * ccomps * dcomps);

            auto g_z_x_xyyz_zz = cbuffer.data(gd_geom_11_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xyyz_xx, g_z_x_xyyz_xy, g_z_x_xyyz_xz, g_z_x_xyyz_yy, g_z_x_xyyz_yz, g_z_x_xyyz_zz, g_z_x_xyz_xx, g_z_x_xyz_xxy, g_z_x_xyz_xy, g_z_x_xyz_xyy, g_z_x_xyz_xyz, g_z_x_xyz_xz, g_z_x_xyz_yy, g_z_x_xyz_yyy, g_z_x_xyz_yyz, g_z_x_xyz_yz, g_z_x_xyz_yzz, g_z_x_xyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyyz_xx[k] = -g_z_x_xyz_xx[k] * ab_y + g_z_x_xyz_xxy[k];

                g_z_x_xyyz_xy[k] = -g_z_x_xyz_xy[k] * ab_y + g_z_x_xyz_xyy[k];

                g_z_x_xyyz_xz[k] = -g_z_x_xyz_xz[k] * ab_y + g_z_x_xyz_xyz[k];

                g_z_x_xyyz_yy[k] = -g_z_x_xyz_yy[k] * ab_y + g_z_x_xyz_yyy[k];

                g_z_x_xyyz_yz[k] = -g_z_x_xyz_yz[k] * ab_y + g_z_x_xyz_yyz[k];

                g_z_x_xyyz_zz[k] = -g_z_x_xyz_zz[k] * ab_y + g_z_x_xyz_yzz[k];
            }

            /// Set up 588-594 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyzz_xx = cbuffer.data(gd_geom_11_off + 588 * ccomps * dcomps);

            auto g_z_x_xyzz_xy = cbuffer.data(gd_geom_11_off + 589 * ccomps * dcomps);

            auto g_z_x_xyzz_xz = cbuffer.data(gd_geom_11_off + 590 * ccomps * dcomps);

            auto g_z_x_xyzz_yy = cbuffer.data(gd_geom_11_off + 591 * ccomps * dcomps);

            auto g_z_x_xyzz_yz = cbuffer.data(gd_geom_11_off + 592 * ccomps * dcomps);

            auto g_z_x_xyzz_zz = cbuffer.data(gd_geom_11_off + 593 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xyzz_xx, g_z_x_xyzz_xy, g_z_x_xyzz_xz, g_z_x_xyzz_yy, g_z_x_xyzz_yz, g_z_x_xyzz_zz, g_z_x_xzz_xx, g_z_x_xzz_xxy, g_z_x_xzz_xy, g_z_x_xzz_xyy, g_z_x_xzz_xyz, g_z_x_xzz_xz, g_z_x_xzz_yy, g_z_x_xzz_yyy, g_z_x_xzz_yyz, g_z_x_xzz_yz, g_z_x_xzz_yzz, g_z_x_xzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyzz_xx[k] = -g_z_x_xzz_xx[k] * ab_y + g_z_x_xzz_xxy[k];

                g_z_x_xyzz_xy[k] = -g_z_x_xzz_xy[k] * ab_y + g_z_x_xzz_xyy[k];

                g_z_x_xyzz_xz[k] = -g_z_x_xzz_xz[k] * ab_y + g_z_x_xzz_xyz[k];

                g_z_x_xyzz_yy[k] = -g_z_x_xzz_yy[k] * ab_y + g_z_x_xzz_yyy[k];

                g_z_x_xyzz_yz[k] = -g_z_x_xzz_yz[k] * ab_y + g_z_x_xzz_yyz[k];

                g_z_x_xyzz_zz[k] = -g_z_x_xzz_zz[k] * ab_y + g_z_x_xzz_yzz[k];
            }

            /// Set up 594-600 components of targeted buffer : cbuffer.data(

            auto g_z_x_xzzz_xx = cbuffer.data(gd_geom_11_off + 594 * ccomps * dcomps);

            auto g_z_x_xzzz_xy = cbuffer.data(gd_geom_11_off + 595 * ccomps * dcomps);

            auto g_z_x_xzzz_xz = cbuffer.data(gd_geom_11_off + 596 * ccomps * dcomps);

            auto g_z_x_xzzz_yy = cbuffer.data(gd_geom_11_off + 597 * ccomps * dcomps);

            auto g_z_x_xzzz_yz = cbuffer.data(gd_geom_11_off + 598 * ccomps * dcomps);

            auto g_z_x_xzzz_zz = cbuffer.data(gd_geom_11_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzz_xx, g_z_0_zzz_xy, g_z_0_zzz_xz, g_z_0_zzz_yy, g_z_0_zzz_yz, g_z_0_zzz_zz, g_z_x_xzzz_xx, g_z_x_xzzz_xy, g_z_x_xzzz_xz, g_z_x_xzzz_yy, g_z_x_xzzz_yz, g_z_x_xzzz_zz, g_z_x_zzz_xx, g_z_x_zzz_xxx, g_z_x_zzz_xxy, g_z_x_zzz_xxz, g_z_x_zzz_xy, g_z_x_zzz_xyy, g_z_x_zzz_xyz, g_z_x_zzz_xz, g_z_x_zzz_xzz, g_z_x_zzz_yy, g_z_x_zzz_yz, g_z_x_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xzzz_xx[k] = g_z_0_zzz_xx[k] - g_z_x_zzz_xx[k] * ab_x + g_z_x_zzz_xxx[k];

                g_z_x_xzzz_xy[k] = g_z_0_zzz_xy[k] - g_z_x_zzz_xy[k] * ab_x + g_z_x_zzz_xxy[k];

                g_z_x_xzzz_xz[k] = g_z_0_zzz_xz[k] - g_z_x_zzz_xz[k] * ab_x + g_z_x_zzz_xxz[k];

                g_z_x_xzzz_yy[k] = g_z_0_zzz_yy[k] - g_z_x_zzz_yy[k] * ab_x + g_z_x_zzz_xyy[k];

                g_z_x_xzzz_yz[k] = g_z_0_zzz_yz[k] - g_z_x_zzz_yz[k] * ab_x + g_z_x_zzz_xyz[k];

                g_z_x_xzzz_zz[k] = g_z_0_zzz_zz[k] - g_z_x_zzz_zz[k] * ab_x + g_z_x_zzz_xzz[k];
            }

            /// Set up 600-606 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyyy_xx = cbuffer.data(gd_geom_11_off + 600 * ccomps * dcomps);

            auto g_z_x_yyyy_xy = cbuffer.data(gd_geom_11_off + 601 * ccomps * dcomps);

            auto g_z_x_yyyy_xz = cbuffer.data(gd_geom_11_off + 602 * ccomps * dcomps);

            auto g_z_x_yyyy_yy = cbuffer.data(gd_geom_11_off + 603 * ccomps * dcomps);

            auto g_z_x_yyyy_yz = cbuffer.data(gd_geom_11_off + 604 * ccomps * dcomps);

            auto g_z_x_yyyy_zz = cbuffer.data(gd_geom_11_off + 605 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yyy_xx, g_z_x_yyy_xxy, g_z_x_yyy_xy, g_z_x_yyy_xyy, g_z_x_yyy_xyz, g_z_x_yyy_xz, g_z_x_yyy_yy, g_z_x_yyy_yyy, g_z_x_yyy_yyz, g_z_x_yyy_yz, g_z_x_yyy_yzz, g_z_x_yyy_zz, g_z_x_yyyy_xx, g_z_x_yyyy_xy, g_z_x_yyyy_xz, g_z_x_yyyy_yy, g_z_x_yyyy_yz, g_z_x_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyyy_xx[k] = -g_z_x_yyy_xx[k] * ab_y + g_z_x_yyy_xxy[k];

                g_z_x_yyyy_xy[k] = -g_z_x_yyy_xy[k] * ab_y + g_z_x_yyy_xyy[k];

                g_z_x_yyyy_xz[k] = -g_z_x_yyy_xz[k] * ab_y + g_z_x_yyy_xyz[k];

                g_z_x_yyyy_yy[k] = -g_z_x_yyy_yy[k] * ab_y + g_z_x_yyy_yyy[k];

                g_z_x_yyyy_yz[k] = -g_z_x_yyy_yz[k] * ab_y + g_z_x_yyy_yyz[k];

                g_z_x_yyyy_zz[k] = -g_z_x_yyy_zz[k] * ab_y + g_z_x_yyy_yzz[k];
            }

            /// Set up 606-612 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyyz_xx = cbuffer.data(gd_geom_11_off + 606 * ccomps * dcomps);

            auto g_z_x_yyyz_xy = cbuffer.data(gd_geom_11_off + 607 * ccomps * dcomps);

            auto g_z_x_yyyz_xz = cbuffer.data(gd_geom_11_off + 608 * ccomps * dcomps);

            auto g_z_x_yyyz_yy = cbuffer.data(gd_geom_11_off + 609 * ccomps * dcomps);

            auto g_z_x_yyyz_yz = cbuffer.data(gd_geom_11_off + 610 * ccomps * dcomps);

            auto g_z_x_yyyz_zz = cbuffer.data(gd_geom_11_off + 611 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yyyz_xx, g_z_x_yyyz_xy, g_z_x_yyyz_xz, g_z_x_yyyz_yy, g_z_x_yyyz_yz, g_z_x_yyyz_zz, g_z_x_yyz_xx, g_z_x_yyz_xxy, g_z_x_yyz_xy, g_z_x_yyz_xyy, g_z_x_yyz_xyz, g_z_x_yyz_xz, g_z_x_yyz_yy, g_z_x_yyz_yyy, g_z_x_yyz_yyz, g_z_x_yyz_yz, g_z_x_yyz_yzz, g_z_x_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyyz_xx[k] = -g_z_x_yyz_xx[k] * ab_y + g_z_x_yyz_xxy[k];

                g_z_x_yyyz_xy[k] = -g_z_x_yyz_xy[k] * ab_y + g_z_x_yyz_xyy[k];

                g_z_x_yyyz_xz[k] = -g_z_x_yyz_xz[k] * ab_y + g_z_x_yyz_xyz[k];

                g_z_x_yyyz_yy[k] = -g_z_x_yyz_yy[k] * ab_y + g_z_x_yyz_yyy[k];

                g_z_x_yyyz_yz[k] = -g_z_x_yyz_yz[k] * ab_y + g_z_x_yyz_yyz[k];

                g_z_x_yyyz_zz[k] = -g_z_x_yyz_zz[k] * ab_y + g_z_x_yyz_yzz[k];
            }

            /// Set up 612-618 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyzz_xx = cbuffer.data(gd_geom_11_off + 612 * ccomps * dcomps);

            auto g_z_x_yyzz_xy = cbuffer.data(gd_geom_11_off + 613 * ccomps * dcomps);

            auto g_z_x_yyzz_xz = cbuffer.data(gd_geom_11_off + 614 * ccomps * dcomps);

            auto g_z_x_yyzz_yy = cbuffer.data(gd_geom_11_off + 615 * ccomps * dcomps);

            auto g_z_x_yyzz_yz = cbuffer.data(gd_geom_11_off + 616 * ccomps * dcomps);

            auto g_z_x_yyzz_zz = cbuffer.data(gd_geom_11_off + 617 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yyzz_xx, g_z_x_yyzz_xy, g_z_x_yyzz_xz, g_z_x_yyzz_yy, g_z_x_yyzz_yz, g_z_x_yyzz_zz, g_z_x_yzz_xx, g_z_x_yzz_xxy, g_z_x_yzz_xy, g_z_x_yzz_xyy, g_z_x_yzz_xyz, g_z_x_yzz_xz, g_z_x_yzz_yy, g_z_x_yzz_yyy, g_z_x_yzz_yyz, g_z_x_yzz_yz, g_z_x_yzz_yzz, g_z_x_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyzz_xx[k] = -g_z_x_yzz_xx[k] * ab_y + g_z_x_yzz_xxy[k];

                g_z_x_yyzz_xy[k] = -g_z_x_yzz_xy[k] * ab_y + g_z_x_yzz_xyy[k];

                g_z_x_yyzz_xz[k] = -g_z_x_yzz_xz[k] * ab_y + g_z_x_yzz_xyz[k];

                g_z_x_yyzz_yy[k] = -g_z_x_yzz_yy[k] * ab_y + g_z_x_yzz_yyy[k];

                g_z_x_yyzz_yz[k] = -g_z_x_yzz_yz[k] * ab_y + g_z_x_yzz_yyz[k];

                g_z_x_yyzz_zz[k] = -g_z_x_yzz_zz[k] * ab_y + g_z_x_yzz_yzz[k];
            }

            /// Set up 618-624 components of targeted buffer : cbuffer.data(

            auto g_z_x_yzzz_xx = cbuffer.data(gd_geom_11_off + 618 * ccomps * dcomps);

            auto g_z_x_yzzz_xy = cbuffer.data(gd_geom_11_off + 619 * ccomps * dcomps);

            auto g_z_x_yzzz_xz = cbuffer.data(gd_geom_11_off + 620 * ccomps * dcomps);

            auto g_z_x_yzzz_yy = cbuffer.data(gd_geom_11_off + 621 * ccomps * dcomps);

            auto g_z_x_yzzz_yz = cbuffer.data(gd_geom_11_off + 622 * ccomps * dcomps);

            auto g_z_x_yzzz_zz = cbuffer.data(gd_geom_11_off + 623 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yzzz_xx, g_z_x_yzzz_xy, g_z_x_yzzz_xz, g_z_x_yzzz_yy, g_z_x_yzzz_yz, g_z_x_yzzz_zz, g_z_x_zzz_xx, g_z_x_zzz_xxy, g_z_x_zzz_xy, g_z_x_zzz_xyy, g_z_x_zzz_xyz, g_z_x_zzz_xz, g_z_x_zzz_yy, g_z_x_zzz_yyy, g_z_x_zzz_yyz, g_z_x_zzz_yz, g_z_x_zzz_yzz, g_z_x_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yzzz_xx[k] = -g_z_x_zzz_xx[k] * ab_y + g_z_x_zzz_xxy[k];

                g_z_x_yzzz_xy[k] = -g_z_x_zzz_xy[k] * ab_y + g_z_x_zzz_xyy[k];

                g_z_x_yzzz_xz[k] = -g_z_x_zzz_xz[k] * ab_y + g_z_x_zzz_xyz[k];

                g_z_x_yzzz_yy[k] = -g_z_x_zzz_yy[k] * ab_y + g_z_x_zzz_yyy[k];

                g_z_x_yzzz_yz[k] = -g_z_x_zzz_yz[k] * ab_y + g_z_x_zzz_yyz[k];

                g_z_x_yzzz_zz[k] = -g_z_x_zzz_zz[k] * ab_y + g_z_x_zzz_yzz[k];
            }

            /// Set up 624-630 components of targeted buffer : cbuffer.data(

            auto g_z_x_zzzz_xx = cbuffer.data(gd_geom_11_off + 624 * ccomps * dcomps);

            auto g_z_x_zzzz_xy = cbuffer.data(gd_geom_11_off + 625 * ccomps * dcomps);

            auto g_z_x_zzzz_xz = cbuffer.data(gd_geom_11_off + 626 * ccomps * dcomps);

            auto g_z_x_zzzz_yy = cbuffer.data(gd_geom_11_off + 627 * ccomps * dcomps);

            auto g_z_x_zzzz_yz = cbuffer.data(gd_geom_11_off + 628 * ccomps * dcomps);

            auto g_z_x_zzzz_zz = cbuffer.data(gd_geom_11_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzz_xx, g_0_x_zzz_xy, g_0_x_zzz_xz, g_0_x_zzz_yy, g_0_x_zzz_yz, g_0_x_zzz_zz, g_z_x_zzz_xx, g_z_x_zzz_xxz, g_z_x_zzz_xy, g_z_x_zzz_xyz, g_z_x_zzz_xz, g_z_x_zzz_xzz, g_z_x_zzz_yy, g_z_x_zzz_yyz, g_z_x_zzz_yz, g_z_x_zzz_yzz, g_z_x_zzz_zz, g_z_x_zzz_zzz, g_z_x_zzzz_xx, g_z_x_zzzz_xy, g_z_x_zzzz_xz, g_z_x_zzzz_yy, g_z_x_zzzz_yz, g_z_x_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_zzzz_xx[k] = -g_0_x_zzz_xx[k] - g_z_x_zzz_xx[k] * ab_z + g_z_x_zzz_xxz[k];

                g_z_x_zzzz_xy[k] = -g_0_x_zzz_xy[k] - g_z_x_zzz_xy[k] * ab_z + g_z_x_zzz_xyz[k];

                g_z_x_zzzz_xz[k] = -g_0_x_zzz_xz[k] - g_z_x_zzz_xz[k] * ab_z + g_z_x_zzz_xzz[k];

                g_z_x_zzzz_yy[k] = -g_0_x_zzz_yy[k] - g_z_x_zzz_yy[k] * ab_z + g_z_x_zzz_yyz[k];

                g_z_x_zzzz_yz[k] = -g_0_x_zzz_yz[k] - g_z_x_zzz_yz[k] * ab_z + g_z_x_zzz_yzz[k];

                g_z_x_zzzz_zz[k] = -g_0_x_zzz_zz[k] - g_z_x_zzz_zz[k] * ab_z + g_z_x_zzz_zzz[k];
            }

            /// Set up 630-636 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxxx_xx = cbuffer.data(gd_geom_11_off + 630 * ccomps * dcomps);

            auto g_z_y_xxxx_xy = cbuffer.data(gd_geom_11_off + 631 * ccomps * dcomps);

            auto g_z_y_xxxx_xz = cbuffer.data(gd_geom_11_off + 632 * ccomps * dcomps);

            auto g_z_y_xxxx_yy = cbuffer.data(gd_geom_11_off + 633 * ccomps * dcomps);

            auto g_z_y_xxxx_yz = cbuffer.data(gd_geom_11_off + 634 * ccomps * dcomps);

            auto g_z_y_xxxx_zz = cbuffer.data(gd_geom_11_off + 635 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxx_xx, g_z_y_xxx_xxx, g_z_y_xxx_xxy, g_z_y_xxx_xxz, g_z_y_xxx_xy, g_z_y_xxx_xyy, g_z_y_xxx_xyz, g_z_y_xxx_xz, g_z_y_xxx_xzz, g_z_y_xxx_yy, g_z_y_xxx_yz, g_z_y_xxx_zz, g_z_y_xxxx_xx, g_z_y_xxxx_xy, g_z_y_xxxx_xz, g_z_y_xxxx_yy, g_z_y_xxxx_yz, g_z_y_xxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxxx_xx[k] = -g_z_y_xxx_xx[k] * ab_x + g_z_y_xxx_xxx[k];

                g_z_y_xxxx_xy[k] = -g_z_y_xxx_xy[k] * ab_x + g_z_y_xxx_xxy[k];

                g_z_y_xxxx_xz[k] = -g_z_y_xxx_xz[k] * ab_x + g_z_y_xxx_xxz[k];

                g_z_y_xxxx_yy[k] = -g_z_y_xxx_yy[k] * ab_x + g_z_y_xxx_xyy[k];

                g_z_y_xxxx_yz[k] = -g_z_y_xxx_yz[k] * ab_x + g_z_y_xxx_xyz[k];

                g_z_y_xxxx_zz[k] = -g_z_y_xxx_zz[k] * ab_x + g_z_y_xxx_xzz[k];
            }

            /// Set up 636-642 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxxy_xx = cbuffer.data(gd_geom_11_off + 636 * ccomps * dcomps);

            auto g_z_y_xxxy_xy = cbuffer.data(gd_geom_11_off + 637 * ccomps * dcomps);

            auto g_z_y_xxxy_xz = cbuffer.data(gd_geom_11_off + 638 * ccomps * dcomps);

            auto g_z_y_xxxy_yy = cbuffer.data(gd_geom_11_off + 639 * ccomps * dcomps);

            auto g_z_y_xxxy_yz = cbuffer.data(gd_geom_11_off + 640 * ccomps * dcomps);

            auto g_z_y_xxxy_zz = cbuffer.data(gd_geom_11_off + 641 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxxy_xx, g_z_y_xxxy_xy, g_z_y_xxxy_xz, g_z_y_xxxy_yy, g_z_y_xxxy_yz, g_z_y_xxxy_zz, g_z_y_xxy_xx, g_z_y_xxy_xxx, g_z_y_xxy_xxy, g_z_y_xxy_xxz, g_z_y_xxy_xy, g_z_y_xxy_xyy, g_z_y_xxy_xyz, g_z_y_xxy_xz, g_z_y_xxy_xzz, g_z_y_xxy_yy, g_z_y_xxy_yz, g_z_y_xxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxxy_xx[k] = -g_z_y_xxy_xx[k] * ab_x + g_z_y_xxy_xxx[k];

                g_z_y_xxxy_xy[k] = -g_z_y_xxy_xy[k] * ab_x + g_z_y_xxy_xxy[k];

                g_z_y_xxxy_xz[k] = -g_z_y_xxy_xz[k] * ab_x + g_z_y_xxy_xxz[k];

                g_z_y_xxxy_yy[k] = -g_z_y_xxy_yy[k] * ab_x + g_z_y_xxy_xyy[k];

                g_z_y_xxxy_yz[k] = -g_z_y_xxy_yz[k] * ab_x + g_z_y_xxy_xyz[k];

                g_z_y_xxxy_zz[k] = -g_z_y_xxy_zz[k] * ab_x + g_z_y_xxy_xzz[k];
            }

            /// Set up 642-648 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxxz_xx = cbuffer.data(gd_geom_11_off + 642 * ccomps * dcomps);

            auto g_z_y_xxxz_xy = cbuffer.data(gd_geom_11_off + 643 * ccomps * dcomps);

            auto g_z_y_xxxz_xz = cbuffer.data(gd_geom_11_off + 644 * ccomps * dcomps);

            auto g_z_y_xxxz_yy = cbuffer.data(gd_geom_11_off + 645 * ccomps * dcomps);

            auto g_z_y_xxxz_yz = cbuffer.data(gd_geom_11_off + 646 * ccomps * dcomps);

            auto g_z_y_xxxz_zz = cbuffer.data(gd_geom_11_off + 647 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxxz_xx, g_z_y_xxxz_xy, g_z_y_xxxz_xz, g_z_y_xxxz_yy, g_z_y_xxxz_yz, g_z_y_xxxz_zz, g_z_y_xxz_xx, g_z_y_xxz_xxx, g_z_y_xxz_xxy, g_z_y_xxz_xxz, g_z_y_xxz_xy, g_z_y_xxz_xyy, g_z_y_xxz_xyz, g_z_y_xxz_xz, g_z_y_xxz_xzz, g_z_y_xxz_yy, g_z_y_xxz_yz, g_z_y_xxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxxz_xx[k] = -g_z_y_xxz_xx[k] * ab_x + g_z_y_xxz_xxx[k];

                g_z_y_xxxz_xy[k] = -g_z_y_xxz_xy[k] * ab_x + g_z_y_xxz_xxy[k];

                g_z_y_xxxz_xz[k] = -g_z_y_xxz_xz[k] * ab_x + g_z_y_xxz_xxz[k];

                g_z_y_xxxz_yy[k] = -g_z_y_xxz_yy[k] * ab_x + g_z_y_xxz_xyy[k];

                g_z_y_xxxz_yz[k] = -g_z_y_xxz_yz[k] * ab_x + g_z_y_xxz_xyz[k];

                g_z_y_xxxz_zz[k] = -g_z_y_xxz_zz[k] * ab_x + g_z_y_xxz_xzz[k];
            }

            /// Set up 648-654 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxyy_xx = cbuffer.data(gd_geom_11_off + 648 * ccomps * dcomps);

            auto g_z_y_xxyy_xy = cbuffer.data(gd_geom_11_off + 649 * ccomps * dcomps);

            auto g_z_y_xxyy_xz = cbuffer.data(gd_geom_11_off + 650 * ccomps * dcomps);

            auto g_z_y_xxyy_yy = cbuffer.data(gd_geom_11_off + 651 * ccomps * dcomps);

            auto g_z_y_xxyy_yz = cbuffer.data(gd_geom_11_off + 652 * ccomps * dcomps);

            auto g_z_y_xxyy_zz = cbuffer.data(gd_geom_11_off + 653 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxyy_xx, g_z_y_xxyy_xy, g_z_y_xxyy_xz, g_z_y_xxyy_yy, g_z_y_xxyy_yz, g_z_y_xxyy_zz, g_z_y_xyy_xx, g_z_y_xyy_xxx, g_z_y_xyy_xxy, g_z_y_xyy_xxz, g_z_y_xyy_xy, g_z_y_xyy_xyy, g_z_y_xyy_xyz, g_z_y_xyy_xz, g_z_y_xyy_xzz, g_z_y_xyy_yy, g_z_y_xyy_yz, g_z_y_xyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxyy_xx[k] = -g_z_y_xyy_xx[k] * ab_x + g_z_y_xyy_xxx[k];

                g_z_y_xxyy_xy[k] = -g_z_y_xyy_xy[k] * ab_x + g_z_y_xyy_xxy[k];

                g_z_y_xxyy_xz[k] = -g_z_y_xyy_xz[k] * ab_x + g_z_y_xyy_xxz[k];

                g_z_y_xxyy_yy[k] = -g_z_y_xyy_yy[k] * ab_x + g_z_y_xyy_xyy[k];

                g_z_y_xxyy_yz[k] = -g_z_y_xyy_yz[k] * ab_x + g_z_y_xyy_xyz[k];

                g_z_y_xxyy_zz[k] = -g_z_y_xyy_zz[k] * ab_x + g_z_y_xyy_xzz[k];
            }

            /// Set up 654-660 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxyz_xx = cbuffer.data(gd_geom_11_off + 654 * ccomps * dcomps);

            auto g_z_y_xxyz_xy = cbuffer.data(gd_geom_11_off + 655 * ccomps * dcomps);

            auto g_z_y_xxyz_xz = cbuffer.data(gd_geom_11_off + 656 * ccomps * dcomps);

            auto g_z_y_xxyz_yy = cbuffer.data(gd_geom_11_off + 657 * ccomps * dcomps);

            auto g_z_y_xxyz_yz = cbuffer.data(gd_geom_11_off + 658 * ccomps * dcomps);

            auto g_z_y_xxyz_zz = cbuffer.data(gd_geom_11_off + 659 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxyz_xx, g_z_y_xxyz_xy, g_z_y_xxyz_xz, g_z_y_xxyz_yy, g_z_y_xxyz_yz, g_z_y_xxyz_zz, g_z_y_xyz_xx, g_z_y_xyz_xxx, g_z_y_xyz_xxy, g_z_y_xyz_xxz, g_z_y_xyz_xy, g_z_y_xyz_xyy, g_z_y_xyz_xyz, g_z_y_xyz_xz, g_z_y_xyz_xzz, g_z_y_xyz_yy, g_z_y_xyz_yz, g_z_y_xyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxyz_xx[k] = -g_z_y_xyz_xx[k] * ab_x + g_z_y_xyz_xxx[k];

                g_z_y_xxyz_xy[k] = -g_z_y_xyz_xy[k] * ab_x + g_z_y_xyz_xxy[k];

                g_z_y_xxyz_xz[k] = -g_z_y_xyz_xz[k] * ab_x + g_z_y_xyz_xxz[k];

                g_z_y_xxyz_yy[k] = -g_z_y_xyz_yy[k] * ab_x + g_z_y_xyz_xyy[k];

                g_z_y_xxyz_yz[k] = -g_z_y_xyz_yz[k] * ab_x + g_z_y_xyz_xyz[k];

                g_z_y_xxyz_zz[k] = -g_z_y_xyz_zz[k] * ab_x + g_z_y_xyz_xzz[k];
            }

            /// Set up 660-666 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxzz_xx = cbuffer.data(gd_geom_11_off + 660 * ccomps * dcomps);

            auto g_z_y_xxzz_xy = cbuffer.data(gd_geom_11_off + 661 * ccomps * dcomps);

            auto g_z_y_xxzz_xz = cbuffer.data(gd_geom_11_off + 662 * ccomps * dcomps);

            auto g_z_y_xxzz_yy = cbuffer.data(gd_geom_11_off + 663 * ccomps * dcomps);

            auto g_z_y_xxzz_yz = cbuffer.data(gd_geom_11_off + 664 * ccomps * dcomps);

            auto g_z_y_xxzz_zz = cbuffer.data(gd_geom_11_off + 665 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxzz_xx, g_z_y_xxzz_xy, g_z_y_xxzz_xz, g_z_y_xxzz_yy, g_z_y_xxzz_yz, g_z_y_xxzz_zz, g_z_y_xzz_xx, g_z_y_xzz_xxx, g_z_y_xzz_xxy, g_z_y_xzz_xxz, g_z_y_xzz_xy, g_z_y_xzz_xyy, g_z_y_xzz_xyz, g_z_y_xzz_xz, g_z_y_xzz_xzz, g_z_y_xzz_yy, g_z_y_xzz_yz, g_z_y_xzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxzz_xx[k] = -g_z_y_xzz_xx[k] * ab_x + g_z_y_xzz_xxx[k];

                g_z_y_xxzz_xy[k] = -g_z_y_xzz_xy[k] * ab_x + g_z_y_xzz_xxy[k];

                g_z_y_xxzz_xz[k] = -g_z_y_xzz_xz[k] * ab_x + g_z_y_xzz_xxz[k];

                g_z_y_xxzz_yy[k] = -g_z_y_xzz_yy[k] * ab_x + g_z_y_xzz_xyy[k];

                g_z_y_xxzz_yz[k] = -g_z_y_xzz_yz[k] * ab_x + g_z_y_xzz_xyz[k];

                g_z_y_xxzz_zz[k] = -g_z_y_xzz_zz[k] * ab_x + g_z_y_xzz_xzz[k];
            }

            /// Set up 666-672 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyyy_xx = cbuffer.data(gd_geom_11_off + 666 * ccomps * dcomps);

            auto g_z_y_xyyy_xy = cbuffer.data(gd_geom_11_off + 667 * ccomps * dcomps);

            auto g_z_y_xyyy_xz = cbuffer.data(gd_geom_11_off + 668 * ccomps * dcomps);

            auto g_z_y_xyyy_yy = cbuffer.data(gd_geom_11_off + 669 * ccomps * dcomps);

            auto g_z_y_xyyy_yz = cbuffer.data(gd_geom_11_off + 670 * ccomps * dcomps);

            auto g_z_y_xyyy_zz = cbuffer.data(gd_geom_11_off + 671 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyyy_xx, g_z_y_xyyy_xy, g_z_y_xyyy_xz, g_z_y_xyyy_yy, g_z_y_xyyy_yz, g_z_y_xyyy_zz, g_z_y_yyy_xx, g_z_y_yyy_xxx, g_z_y_yyy_xxy, g_z_y_yyy_xxz, g_z_y_yyy_xy, g_z_y_yyy_xyy, g_z_y_yyy_xyz, g_z_y_yyy_xz, g_z_y_yyy_xzz, g_z_y_yyy_yy, g_z_y_yyy_yz, g_z_y_yyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyyy_xx[k] = -g_z_y_yyy_xx[k] * ab_x + g_z_y_yyy_xxx[k];

                g_z_y_xyyy_xy[k] = -g_z_y_yyy_xy[k] * ab_x + g_z_y_yyy_xxy[k];

                g_z_y_xyyy_xz[k] = -g_z_y_yyy_xz[k] * ab_x + g_z_y_yyy_xxz[k];

                g_z_y_xyyy_yy[k] = -g_z_y_yyy_yy[k] * ab_x + g_z_y_yyy_xyy[k];

                g_z_y_xyyy_yz[k] = -g_z_y_yyy_yz[k] * ab_x + g_z_y_yyy_xyz[k];

                g_z_y_xyyy_zz[k] = -g_z_y_yyy_zz[k] * ab_x + g_z_y_yyy_xzz[k];
            }

            /// Set up 672-678 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyyz_xx = cbuffer.data(gd_geom_11_off + 672 * ccomps * dcomps);

            auto g_z_y_xyyz_xy = cbuffer.data(gd_geom_11_off + 673 * ccomps * dcomps);

            auto g_z_y_xyyz_xz = cbuffer.data(gd_geom_11_off + 674 * ccomps * dcomps);

            auto g_z_y_xyyz_yy = cbuffer.data(gd_geom_11_off + 675 * ccomps * dcomps);

            auto g_z_y_xyyz_yz = cbuffer.data(gd_geom_11_off + 676 * ccomps * dcomps);

            auto g_z_y_xyyz_zz = cbuffer.data(gd_geom_11_off + 677 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyyz_xx, g_z_y_xyyz_xy, g_z_y_xyyz_xz, g_z_y_xyyz_yy, g_z_y_xyyz_yz, g_z_y_xyyz_zz, g_z_y_yyz_xx, g_z_y_yyz_xxx, g_z_y_yyz_xxy, g_z_y_yyz_xxz, g_z_y_yyz_xy, g_z_y_yyz_xyy, g_z_y_yyz_xyz, g_z_y_yyz_xz, g_z_y_yyz_xzz, g_z_y_yyz_yy, g_z_y_yyz_yz, g_z_y_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyyz_xx[k] = -g_z_y_yyz_xx[k] * ab_x + g_z_y_yyz_xxx[k];

                g_z_y_xyyz_xy[k] = -g_z_y_yyz_xy[k] * ab_x + g_z_y_yyz_xxy[k];

                g_z_y_xyyz_xz[k] = -g_z_y_yyz_xz[k] * ab_x + g_z_y_yyz_xxz[k];

                g_z_y_xyyz_yy[k] = -g_z_y_yyz_yy[k] * ab_x + g_z_y_yyz_xyy[k];

                g_z_y_xyyz_yz[k] = -g_z_y_yyz_yz[k] * ab_x + g_z_y_yyz_xyz[k];

                g_z_y_xyyz_zz[k] = -g_z_y_yyz_zz[k] * ab_x + g_z_y_yyz_xzz[k];
            }

            /// Set up 678-684 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyzz_xx = cbuffer.data(gd_geom_11_off + 678 * ccomps * dcomps);

            auto g_z_y_xyzz_xy = cbuffer.data(gd_geom_11_off + 679 * ccomps * dcomps);

            auto g_z_y_xyzz_xz = cbuffer.data(gd_geom_11_off + 680 * ccomps * dcomps);

            auto g_z_y_xyzz_yy = cbuffer.data(gd_geom_11_off + 681 * ccomps * dcomps);

            auto g_z_y_xyzz_yz = cbuffer.data(gd_geom_11_off + 682 * ccomps * dcomps);

            auto g_z_y_xyzz_zz = cbuffer.data(gd_geom_11_off + 683 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyzz_xx, g_z_y_xyzz_xy, g_z_y_xyzz_xz, g_z_y_xyzz_yy, g_z_y_xyzz_yz, g_z_y_xyzz_zz, g_z_y_yzz_xx, g_z_y_yzz_xxx, g_z_y_yzz_xxy, g_z_y_yzz_xxz, g_z_y_yzz_xy, g_z_y_yzz_xyy, g_z_y_yzz_xyz, g_z_y_yzz_xz, g_z_y_yzz_xzz, g_z_y_yzz_yy, g_z_y_yzz_yz, g_z_y_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyzz_xx[k] = -g_z_y_yzz_xx[k] * ab_x + g_z_y_yzz_xxx[k];

                g_z_y_xyzz_xy[k] = -g_z_y_yzz_xy[k] * ab_x + g_z_y_yzz_xxy[k];

                g_z_y_xyzz_xz[k] = -g_z_y_yzz_xz[k] * ab_x + g_z_y_yzz_xxz[k];

                g_z_y_xyzz_yy[k] = -g_z_y_yzz_yy[k] * ab_x + g_z_y_yzz_xyy[k];

                g_z_y_xyzz_yz[k] = -g_z_y_yzz_yz[k] * ab_x + g_z_y_yzz_xyz[k];

                g_z_y_xyzz_zz[k] = -g_z_y_yzz_zz[k] * ab_x + g_z_y_yzz_xzz[k];
            }

            /// Set up 684-690 components of targeted buffer : cbuffer.data(

            auto g_z_y_xzzz_xx = cbuffer.data(gd_geom_11_off + 684 * ccomps * dcomps);

            auto g_z_y_xzzz_xy = cbuffer.data(gd_geom_11_off + 685 * ccomps * dcomps);

            auto g_z_y_xzzz_xz = cbuffer.data(gd_geom_11_off + 686 * ccomps * dcomps);

            auto g_z_y_xzzz_yy = cbuffer.data(gd_geom_11_off + 687 * ccomps * dcomps);

            auto g_z_y_xzzz_yz = cbuffer.data(gd_geom_11_off + 688 * ccomps * dcomps);

            auto g_z_y_xzzz_zz = cbuffer.data(gd_geom_11_off + 689 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xzzz_xx, g_z_y_xzzz_xy, g_z_y_xzzz_xz, g_z_y_xzzz_yy, g_z_y_xzzz_yz, g_z_y_xzzz_zz, g_z_y_zzz_xx, g_z_y_zzz_xxx, g_z_y_zzz_xxy, g_z_y_zzz_xxz, g_z_y_zzz_xy, g_z_y_zzz_xyy, g_z_y_zzz_xyz, g_z_y_zzz_xz, g_z_y_zzz_xzz, g_z_y_zzz_yy, g_z_y_zzz_yz, g_z_y_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xzzz_xx[k] = -g_z_y_zzz_xx[k] * ab_x + g_z_y_zzz_xxx[k];

                g_z_y_xzzz_xy[k] = -g_z_y_zzz_xy[k] * ab_x + g_z_y_zzz_xxy[k];

                g_z_y_xzzz_xz[k] = -g_z_y_zzz_xz[k] * ab_x + g_z_y_zzz_xxz[k];

                g_z_y_xzzz_yy[k] = -g_z_y_zzz_yy[k] * ab_x + g_z_y_zzz_xyy[k];

                g_z_y_xzzz_yz[k] = -g_z_y_zzz_yz[k] * ab_x + g_z_y_zzz_xyz[k];

                g_z_y_xzzz_zz[k] = -g_z_y_zzz_zz[k] * ab_x + g_z_y_zzz_xzz[k];
            }

            /// Set up 690-696 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyyy_xx = cbuffer.data(gd_geom_11_off + 690 * ccomps * dcomps);

            auto g_z_y_yyyy_xy = cbuffer.data(gd_geom_11_off + 691 * ccomps * dcomps);

            auto g_z_y_yyyy_xz = cbuffer.data(gd_geom_11_off + 692 * ccomps * dcomps);

            auto g_z_y_yyyy_yy = cbuffer.data(gd_geom_11_off + 693 * ccomps * dcomps);

            auto g_z_y_yyyy_yz = cbuffer.data(gd_geom_11_off + 694 * ccomps * dcomps);

            auto g_z_y_yyyy_zz = cbuffer.data(gd_geom_11_off + 695 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyy_xx, g_z_0_yyy_xy, g_z_0_yyy_xz, g_z_0_yyy_yy, g_z_0_yyy_yz, g_z_0_yyy_zz, g_z_y_yyy_xx, g_z_y_yyy_xxy, g_z_y_yyy_xy, g_z_y_yyy_xyy, g_z_y_yyy_xyz, g_z_y_yyy_xz, g_z_y_yyy_yy, g_z_y_yyy_yyy, g_z_y_yyy_yyz, g_z_y_yyy_yz, g_z_y_yyy_yzz, g_z_y_yyy_zz, g_z_y_yyyy_xx, g_z_y_yyyy_xy, g_z_y_yyyy_xz, g_z_y_yyyy_yy, g_z_y_yyyy_yz, g_z_y_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyyy_xx[k] = g_z_0_yyy_xx[k] - g_z_y_yyy_xx[k] * ab_y + g_z_y_yyy_xxy[k];

                g_z_y_yyyy_xy[k] = g_z_0_yyy_xy[k] - g_z_y_yyy_xy[k] * ab_y + g_z_y_yyy_xyy[k];

                g_z_y_yyyy_xz[k] = g_z_0_yyy_xz[k] - g_z_y_yyy_xz[k] * ab_y + g_z_y_yyy_xyz[k];

                g_z_y_yyyy_yy[k] = g_z_0_yyy_yy[k] - g_z_y_yyy_yy[k] * ab_y + g_z_y_yyy_yyy[k];

                g_z_y_yyyy_yz[k] = g_z_0_yyy_yz[k] - g_z_y_yyy_yz[k] * ab_y + g_z_y_yyy_yyz[k];

                g_z_y_yyyy_zz[k] = g_z_0_yyy_zz[k] - g_z_y_yyy_zz[k] * ab_y + g_z_y_yyy_yzz[k];
            }

            /// Set up 696-702 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyyz_xx = cbuffer.data(gd_geom_11_off + 696 * ccomps * dcomps);

            auto g_z_y_yyyz_xy = cbuffer.data(gd_geom_11_off + 697 * ccomps * dcomps);

            auto g_z_y_yyyz_xz = cbuffer.data(gd_geom_11_off + 698 * ccomps * dcomps);

            auto g_z_y_yyyz_yy = cbuffer.data(gd_geom_11_off + 699 * ccomps * dcomps);

            auto g_z_y_yyyz_yz = cbuffer.data(gd_geom_11_off + 700 * ccomps * dcomps);

            auto g_z_y_yyyz_zz = cbuffer.data(gd_geom_11_off + 701 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyz_xx, g_z_0_yyz_xy, g_z_0_yyz_xz, g_z_0_yyz_yy, g_z_0_yyz_yz, g_z_0_yyz_zz, g_z_y_yyyz_xx, g_z_y_yyyz_xy, g_z_y_yyyz_xz, g_z_y_yyyz_yy, g_z_y_yyyz_yz, g_z_y_yyyz_zz, g_z_y_yyz_xx, g_z_y_yyz_xxy, g_z_y_yyz_xy, g_z_y_yyz_xyy, g_z_y_yyz_xyz, g_z_y_yyz_xz, g_z_y_yyz_yy, g_z_y_yyz_yyy, g_z_y_yyz_yyz, g_z_y_yyz_yz, g_z_y_yyz_yzz, g_z_y_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyyz_xx[k] = g_z_0_yyz_xx[k] - g_z_y_yyz_xx[k] * ab_y + g_z_y_yyz_xxy[k];

                g_z_y_yyyz_xy[k] = g_z_0_yyz_xy[k] - g_z_y_yyz_xy[k] * ab_y + g_z_y_yyz_xyy[k];

                g_z_y_yyyz_xz[k] = g_z_0_yyz_xz[k] - g_z_y_yyz_xz[k] * ab_y + g_z_y_yyz_xyz[k];

                g_z_y_yyyz_yy[k] = g_z_0_yyz_yy[k] - g_z_y_yyz_yy[k] * ab_y + g_z_y_yyz_yyy[k];

                g_z_y_yyyz_yz[k] = g_z_0_yyz_yz[k] - g_z_y_yyz_yz[k] * ab_y + g_z_y_yyz_yyz[k];

                g_z_y_yyyz_zz[k] = g_z_0_yyz_zz[k] - g_z_y_yyz_zz[k] * ab_y + g_z_y_yyz_yzz[k];
            }

            /// Set up 702-708 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyzz_xx = cbuffer.data(gd_geom_11_off + 702 * ccomps * dcomps);

            auto g_z_y_yyzz_xy = cbuffer.data(gd_geom_11_off + 703 * ccomps * dcomps);

            auto g_z_y_yyzz_xz = cbuffer.data(gd_geom_11_off + 704 * ccomps * dcomps);

            auto g_z_y_yyzz_yy = cbuffer.data(gd_geom_11_off + 705 * ccomps * dcomps);

            auto g_z_y_yyzz_yz = cbuffer.data(gd_geom_11_off + 706 * ccomps * dcomps);

            auto g_z_y_yyzz_zz = cbuffer.data(gd_geom_11_off + 707 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yzz_xx, g_z_0_yzz_xy, g_z_0_yzz_xz, g_z_0_yzz_yy, g_z_0_yzz_yz, g_z_0_yzz_zz, g_z_y_yyzz_xx, g_z_y_yyzz_xy, g_z_y_yyzz_xz, g_z_y_yyzz_yy, g_z_y_yyzz_yz, g_z_y_yyzz_zz, g_z_y_yzz_xx, g_z_y_yzz_xxy, g_z_y_yzz_xy, g_z_y_yzz_xyy, g_z_y_yzz_xyz, g_z_y_yzz_xz, g_z_y_yzz_yy, g_z_y_yzz_yyy, g_z_y_yzz_yyz, g_z_y_yzz_yz, g_z_y_yzz_yzz, g_z_y_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyzz_xx[k] = g_z_0_yzz_xx[k] - g_z_y_yzz_xx[k] * ab_y + g_z_y_yzz_xxy[k];

                g_z_y_yyzz_xy[k] = g_z_0_yzz_xy[k] - g_z_y_yzz_xy[k] * ab_y + g_z_y_yzz_xyy[k];

                g_z_y_yyzz_xz[k] = g_z_0_yzz_xz[k] - g_z_y_yzz_xz[k] * ab_y + g_z_y_yzz_xyz[k];

                g_z_y_yyzz_yy[k] = g_z_0_yzz_yy[k] - g_z_y_yzz_yy[k] * ab_y + g_z_y_yzz_yyy[k];

                g_z_y_yyzz_yz[k] = g_z_0_yzz_yz[k] - g_z_y_yzz_yz[k] * ab_y + g_z_y_yzz_yyz[k];

                g_z_y_yyzz_zz[k] = g_z_0_yzz_zz[k] - g_z_y_yzz_zz[k] * ab_y + g_z_y_yzz_yzz[k];
            }

            /// Set up 708-714 components of targeted buffer : cbuffer.data(

            auto g_z_y_yzzz_xx = cbuffer.data(gd_geom_11_off + 708 * ccomps * dcomps);

            auto g_z_y_yzzz_xy = cbuffer.data(gd_geom_11_off + 709 * ccomps * dcomps);

            auto g_z_y_yzzz_xz = cbuffer.data(gd_geom_11_off + 710 * ccomps * dcomps);

            auto g_z_y_yzzz_yy = cbuffer.data(gd_geom_11_off + 711 * ccomps * dcomps);

            auto g_z_y_yzzz_yz = cbuffer.data(gd_geom_11_off + 712 * ccomps * dcomps);

            auto g_z_y_yzzz_zz = cbuffer.data(gd_geom_11_off + 713 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzz_xx, g_z_0_zzz_xy, g_z_0_zzz_xz, g_z_0_zzz_yy, g_z_0_zzz_yz, g_z_0_zzz_zz, g_z_y_yzzz_xx, g_z_y_yzzz_xy, g_z_y_yzzz_xz, g_z_y_yzzz_yy, g_z_y_yzzz_yz, g_z_y_yzzz_zz, g_z_y_zzz_xx, g_z_y_zzz_xxy, g_z_y_zzz_xy, g_z_y_zzz_xyy, g_z_y_zzz_xyz, g_z_y_zzz_xz, g_z_y_zzz_yy, g_z_y_zzz_yyy, g_z_y_zzz_yyz, g_z_y_zzz_yz, g_z_y_zzz_yzz, g_z_y_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yzzz_xx[k] = g_z_0_zzz_xx[k] - g_z_y_zzz_xx[k] * ab_y + g_z_y_zzz_xxy[k];

                g_z_y_yzzz_xy[k] = g_z_0_zzz_xy[k] - g_z_y_zzz_xy[k] * ab_y + g_z_y_zzz_xyy[k];

                g_z_y_yzzz_xz[k] = g_z_0_zzz_xz[k] - g_z_y_zzz_xz[k] * ab_y + g_z_y_zzz_xyz[k];

                g_z_y_yzzz_yy[k] = g_z_0_zzz_yy[k] - g_z_y_zzz_yy[k] * ab_y + g_z_y_zzz_yyy[k];

                g_z_y_yzzz_yz[k] = g_z_0_zzz_yz[k] - g_z_y_zzz_yz[k] * ab_y + g_z_y_zzz_yyz[k];

                g_z_y_yzzz_zz[k] = g_z_0_zzz_zz[k] - g_z_y_zzz_zz[k] * ab_y + g_z_y_zzz_yzz[k];
            }

            /// Set up 714-720 components of targeted buffer : cbuffer.data(

            auto g_z_y_zzzz_xx = cbuffer.data(gd_geom_11_off + 714 * ccomps * dcomps);

            auto g_z_y_zzzz_xy = cbuffer.data(gd_geom_11_off + 715 * ccomps * dcomps);

            auto g_z_y_zzzz_xz = cbuffer.data(gd_geom_11_off + 716 * ccomps * dcomps);

            auto g_z_y_zzzz_yy = cbuffer.data(gd_geom_11_off + 717 * ccomps * dcomps);

            auto g_z_y_zzzz_yz = cbuffer.data(gd_geom_11_off + 718 * ccomps * dcomps);

            auto g_z_y_zzzz_zz = cbuffer.data(gd_geom_11_off + 719 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzz_xx, g_0_y_zzz_xy, g_0_y_zzz_xz, g_0_y_zzz_yy, g_0_y_zzz_yz, g_0_y_zzz_zz, g_z_y_zzz_xx, g_z_y_zzz_xxz, g_z_y_zzz_xy, g_z_y_zzz_xyz, g_z_y_zzz_xz, g_z_y_zzz_xzz, g_z_y_zzz_yy, g_z_y_zzz_yyz, g_z_y_zzz_yz, g_z_y_zzz_yzz, g_z_y_zzz_zz, g_z_y_zzz_zzz, g_z_y_zzzz_xx, g_z_y_zzzz_xy, g_z_y_zzzz_xz, g_z_y_zzzz_yy, g_z_y_zzzz_yz, g_z_y_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_zzzz_xx[k] = -g_0_y_zzz_xx[k] - g_z_y_zzz_xx[k] * ab_z + g_z_y_zzz_xxz[k];

                g_z_y_zzzz_xy[k] = -g_0_y_zzz_xy[k] - g_z_y_zzz_xy[k] * ab_z + g_z_y_zzz_xyz[k];

                g_z_y_zzzz_xz[k] = -g_0_y_zzz_xz[k] - g_z_y_zzz_xz[k] * ab_z + g_z_y_zzz_xzz[k];

                g_z_y_zzzz_yy[k] = -g_0_y_zzz_yy[k] - g_z_y_zzz_yy[k] * ab_z + g_z_y_zzz_yyz[k];

                g_z_y_zzzz_yz[k] = -g_0_y_zzz_yz[k] - g_z_y_zzz_yz[k] * ab_z + g_z_y_zzz_yzz[k];

                g_z_y_zzzz_zz[k] = -g_0_y_zzz_zz[k] - g_z_y_zzz_zz[k] * ab_z + g_z_y_zzz_zzz[k];
            }

            /// Set up 720-726 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxxx_xx = cbuffer.data(gd_geom_11_off + 720 * ccomps * dcomps);

            auto g_z_z_xxxx_xy = cbuffer.data(gd_geom_11_off + 721 * ccomps * dcomps);

            auto g_z_z_xxxx_xz = cbuffer.data(gd_geom_11_off + 722 * ccomps * dcomps);

            auto g_z_z_xxxx_yy = cbuffer.data(gd_geom_11_off + 723 * ccomps * dcomps);

            auto g_z_z_xxxx_yz = cbuffer.data(gd_geom_11_off + 724 * ccomps * dcomps);

            auto g_z_z_xxxx_zz = cbuffer.data(gd_geom_11_off + 725 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxx_xx, g_z_z_xxx_xxx, g_z_z_xxx_xxy, g_z_z_xxx_xxz, g_z_z_xxx_xy, g_z_z_xxx_xyy, g_z_z_xxx_xyz, g_z_z_xxx_xz, g_z_z_xxx_xzz, g_z_z_xxx_yy, g_z_z_xxx_yz, g_z_z_xxx_zz, g_z_z_xxxx_xx, g_z_z_xxxx_xy, g_z_z_xxxx_xz, g_z_z_xxxx_yy, g_z_z_xxxx_yz, g_z_z_xxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxxx_xx[k] = -g_z_z_xxx_xx[k] * ab_x + g_z_z_xxx_xxx[k];

                g_z_z_xxxx_xy[k] = -g_z_z_xxx_xy[k] * ab_x + g_z_z_xxx_xxy[k];

                g_z_z_xxxx_xz[k] = -g_z_z_xxx_xz[k] * ab_x + g_z_z_xxx_xxz[k];

                g_z_z_xxxx_yy[k] = -g_z_z_xxx_yy[k] * ab_x + g_z_z_xxx_xyy[k];

                g_z_z_xxxx_yz[k] = -g_z_z_xxx_yz[k] * ab_x + g_z_z_xxx_xyz[k];

                g_z_z_xxxx_zz[k] = -g_z_z_xxx_zz[k] * ab_x + g_z_z_xxx_xzz[k];
            }

            /// Set up 726-732 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxxy_xx = cbuffer.data(gd_geom_11_off + 726 * ccomps * dcomps);

            auto g_z_z_xxxy_xy = cbuffer.data(gd_geom_11_off + 727 * ccomps * dcomps);

            auto g_z_z_xxxy_xz = cbuffer.data(gd_geom_11_off + 728 * ccomps * dcomps);

            auto g_z_z_xxxy_yy = cbuffer.data(gd_geom_11_off + 729 * ccomps * dcomps);

            auto g_z_z_xxxy_yz = cbuffer.data(gd_geom_11_off + 730 * ccomps * dcomps);

            auto g_z_z_xxxy_zz = cbuffer.data(gd_geom_11_off + 731 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxxy_xx, g_z_z_xxxy_xy, g_z_z_xxxy_xz, g_z_z_xxxy_yy, g_z_z_xxxy_yz, g_z_z_xxxy_zz, g_z_z_xxy_xx, g_z_z_xxy_xxx, g_z_z_xxy_xxy, g_z_z_xxy_xxz, g_z_z_xxy_xy, g_z_z_xxy_xyy, g_z_z_xxy_xyz, g_z_z_xxy_xz, g_z_z_xxy_xzz, g_z_z_xxy_yy, g_z_z_xxy_yz, g_z_z_xxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxxy_xx[k] = -g_z_z_xxy_xx[k] * ab_x + g_z_z_xxy_xxx[k];

                g_z_z_xxxy_xy[k] = -g_z_z_xxy_xy[k] * ab_x + g_z_z_xxy_xxy[k];

                g_z_z_xxxy_xz[k] = -g_z_z_xxy_xz[k] * ab_x + g_z_z_xxy_xxz[k];

                g_z_z_xxxy_yy[k] = -g_z_z_xxy_yy[k] * ab_x + g_z_z_xxy_xyy[k];

                g_z_z_xxxy_yz[k] = -g_z_z_xxy_yz[k] * ab_x + g_z_z_xxy_xyz[k];

                g_z_z_xxxy_zz[k] = -g_z_z_xxy_zz[k] * ab_x + g_z_z_xxy_xzz[k];
            }

            /// Set up 732-738 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxxz_xx = cbuffer.data(gd_geom_11_off + 732 * ccomps * dcomps);

            auto g_z_z_xxxz_xy = cbuffer.data(gd_geom_11_off + 733 * ccomps * dcomps);

            auto g_z_z_xxxz_xz = cbuffer.data(gd_geom_11_off + 734 * ccomps * dcomps);

            auto g_z_z_xxxz_yy = cbuffer.data(gd_geom_11_off + 735 * ccomps * dcomps);

            auto g_z_z_xxxz_yz = cbuffer.data(gd_geom_11_off + 736 * ccomps * dcomps);

            auto g_z_z_xxxz_zz = cbuffer.data(gd_geom_11_off + 737 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxxz_xx, g_z_z_xxxz_xy, g_z_z_xxxz_xz, g_z_z_xxxz_yy, g_z_z_xxxz_yz, g_z_z_xxxz_zz, g_z_z_xxz_xx, g_z_z_xxz_xxx, g_z_z_xxz_xxy, g_z_z_xxz_xxz, g_z_z_xxz_xy, g_z_z_xxz_xyy, g_z_z_xxz_xyz, g_z_z_xxz_xz, g_z_z_xxz_xzz, g_z_z_xxz_yy, g_z_z_xxz_yz, g_z_z_xxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxxz_xx[k] = -g_z_z_xxz_xx[k] * ab_x + g_z_z_xxz_xxx[k];

                g_z_z_xxxz_xy[k] = -g_z_z_xxz_xy[k] * ab_x + g_z_z_xxz_xxy[k];

                g_z_z_xxxz_xz[k] = -g_z_z_xxz_xz[k] * ab_x + g_z_z_xxz_xxz[k];

                g_z_z_xxxz_yy[k] = -g_z_z_xxz_yy[k] * ab_x + g_z_z_xxz_xyy[k];

                g_z_z_xxxz_yz[k] = -g_z_z_xxz_yz[k] * ab_x + g_z_z_xxz_xyz[k];

                g_z_z_xxxz_zz[k] = -g_z_z_xxz_zz[k] * ab_x + g_z_z_xxz_xzz[k];
            }

            /// Set up 738-744 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxyy_xx = cbuffer.data(gd_geom_11_off + 738 * ccomps * dcomps);

            auto g_z_z_xxyy_xy = cbuffer.data(gd_geom_11_off + 739 * ccomps * dcomps);

            auto g_z_z_xxyy_xz = cbuffer.data(gd_geom_11_off + 740 * ccomps * dcomps);

            auto g_z_z_xxyy_yy = cbuffer.data(gd_geom_11_off + 741 * ccomps * dcomps);

            auto g_z_z_xxyy_yz = cbuffer.data(gd_geom_11_off + 742 * ccomps * dcomps);

            auto g_z_z_xxyy_zz = cbuffer.data(gd_geom_11_off + 743 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxyy_xx, g_z_z_xxyy_xy, g_z_z_xxyy_xz, g_z_z_xxyy_yy, g_z_z_xxyy_yz, g_z_z_xxyy_zz, g_z_z_xyy_xx, g_z_z_xyy_xxx, g_z_z_xyy_xxy, g_z_z_xyy_xxz, g_z_z_xyy_xy, g_z_z_xyy_xyy, g_z_z_xyy_xyz, g_z_z_xyy_xz, g_z_z_xyy_xzz, g_z_z_xyy_yy, g_z_z_xyy_yz, g_z_z_xyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxyy_xx[k] = -g_z_z_xyy_xx[k] * ab_x + g_z_z_xyy_xxx[k];

                g_z_z_xxyy_xy[k] = -g_z_z_xyy_xy[k] * ab_x + g_z_z_xyy_xxy[k];

                g_z_z_xxyy_xz[k] = -g_z_z_xyy_xz[k] * ab_x + g_z_z_xyy_xxz[k];

                g_z_z_xxyy_yy[k] = -g_z_z_xyy_yy[k] * ab_x + g_z_z_xyy_xyy[k];

                g_z_z_xxyy_yz[k] = -g_z_z_xyy_yz[k] * ab_x + g_z_z_xyy_xyz[k];

                g_z_z_xxyy_zz[k] = -g_z_z_xyy_zz[k] * ab_x + g_z_z_xyy_xzz[k];
            }

            /// Set up 744-750 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxyz_xx = cbuffer.data(gd_geom_11_off + 744 * ccomps * dcomps);

            auto g_z_z_xxyz_xy = cbuffer.data(gd_geom_11_off + 745 * ccomps * dcomps);

            auto g_z_z_xxyz_xz = cbuffer.data(gd_geom_11_off + 746 * ccomps * dcomps);

            auto g_z_z_xxyz_yy = cbuffer.data(gd_geom_11_off + 747 * ccomps * dcomps);

            auto g_z_z_xxyz_yz = cbuffer.data(gd_geom_11_off + 748 * ccomps * dcomps);

            auto g_z_z_xxyz_zz = cbuffer.data(gd_geom_11_off + 749 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxyz_xx, g_z_z_xxyz_xy, g_z_z_xxyz_xz, g_z_z_xxyz_yy, g_z_z_xxyz_yz, g_z_z_xxyz_zz, g_z_z_xyz_xx, g_z_z_xyz_xxx, g_z_z_xyz_xxy, g_z_z_xyz_xxz, g_z_z_xyz_xy, g_z_z_xyz_xyy, g_z_z_xyz_xyz, g_z_z_xyz_xz, g_z_z_xyz_xzz, g_z_z_xyz_yy, g_z_z_xyz_yz, g_z_z_xyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxyz_xx[k] = -g_z_z_xyz_xx[k] * ab_x + g_z_z_xyz_xxx[k];

                g_z_z_xxyz_xy[k] = -g_z_z_xyz_xy[k] * ab_x + g_z_z_xyz_xxy[k];

                g_z_z_xxyz_xz[k] = -g_z_z_xyz_xz[k] * ab_x + g_z_z_xyz_xxz[k];

                g_z_z_xxyz_yy[k] = -g_z_z_xyz_yy[k] * ab_x + g_z_z_xyz_xyy[k];

                g_z_z_xxyz_yz[k] = -g_z_z_xyz_yz[k] * ab_x + g_z_z_xyz_xyz[k];

                g_z_z_xxyz_zz[k] = -g_z_z_xyz_zz[k] * ab_x + g_z_z_xyz_xzz[k];
            }

            /// Set up 750-756 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxzz_xx = cbuffer.data(gd_geom_11_off + 750 * ccomps * dcomps);

            auto g_z_z_xxzz_xy = cbuffer.data(gd_geom_11_off + 751 * ccomps * dcomps);

            auto g_z_z_xxzz_xz = cbuffer.data(gd_geom_11_off + 752 * ccomps * dcomps);

            auto g_z_z_xxzz_yy = cbuffer.data(gd_geom_11_off + 753 * ccomps * dcomps);

            auto g_z_z_xxzz_yz = cbuffer.data(gd_geom_11_off + 754 * ccomps * dcomps);

            auto g_z_z_xxzz_zz = cbuffer.data(gd_geom_11_off + 755 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxzz_xx, g_z_z_xxzz_xy, g_z_z_xxzz_xz, g_z_z_xxzz_yy, g_z_z_xxzz_yz, g_z_z_xxzz_zz, g_z_z_xzz_xx, g_z_z_xzz_xxx, g_z_z_xzz_xxy, g_z_z_xzz_xxz, g_z_z_xzz_xy, g_z_z_xzz_xyy, g_z_z_xzz_xyz, g_z_z_xzz_xz, g_z_z_xzz_xzz, g_z_z_xzz_yy, g_z_z_xzz_yz, g_z_z_xzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxzz_xx[k] = -g_z_z_xzz_xx[k] * ab_x + g_z_z_xzz_xxx[k];

                g_z_z_xxzz_xy[k] = -g_z_z_xzz_xy[k] * ab_x + g_z_z_xzz_xxy[k];

                g_z_z_xxzz_xz[k] = -g_z_z_xzz_xz[k] * ab_x + g_z_z_xzz_xxz[k];

                g_z_z_xxzz_yy[k] = -g_z_z_xzz_yy[k] * ab_x + g_z_z_xzz_xyy[k];

                g_z_z_xxzz_yz[k] = -g_z_z_xzz_yz[k] * ab_x + g_z_z_xzz_xyz[k];

                g_z_z_xxzz_zz[k] = -g_z_z_xzz_zz[k] * ab_x + g_z_z_xzz_xzz[k];
            }

            /// Set up 756-762 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyyy_xx = cbuffer.data(gd_geom_11_off + 756 * ccomps * dcomps);

            auto g_z_z_xyyy_xy = cbuffer.data(gd_geom_11_off + 757 * ccomps * dcomps);

            auto g_z_z_xyyy_xz = cbuffer.data(gd_geom_11_off + 758 * ccomps * dcomps);

            auto g_z_z_xyyy_yy = cbuffer.data(gd_geom_11_off + 759 * ccomps * dcomps);

            auto g_z_z_xyyy_yz = cbuffer.data(gd_geom_11_off + 760 * ccomps * dcomps);

            auto g_z_z_xyyy_zz = cbuffer.data(gd_geom_11_off + 761 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyyy_xx, g_z_z_xyyy_xy, g_z_z_xyyy_xz, g_z_z_xyyy_yy, g_z_z_xyyy_yz, g_z_z_xyyy_zz, g_z_z_yyy_xx, g_z_z_yyy_xxx, g_z_z_yyy_xxy, g_z_z_yyy_xxz, g_z_z_yyy_xy, g_z_z_yyy_xyy, g_z_z_yyy_xyz, g_z_z_yyy_xz, g_z_z_yyy_xzz, g_z_z_yyy_yy, g_z_z_yyy_yz, g_z_z_yyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyyy_xx[k] = -g_z_z_yyy_xx[k] * ab_x + g_z_z_yyy_xxx[k];

                g_z_z_xyyy_xy[k] = -g_z_z_yyy_xy[k] * ab_x + g_z_z_yyy_xxy[k];

                g_z_z_xyyy_xz[k] = -g_z_z_yyy_xz[k] * ab_x + g_z_z_yyy_xxz[k];

                g_z_z_xyyy_yy[k] = -g_z_z_yyy_yy[k] * ab_x + g_z_z_yyy_xyy[k];

                g_z_z_xyyy_yz[k] = -g_z_z_yyy_yz[k] * ab_x + g_z_z_yyy_xyz[k];

                g_z_z_xyyy_zz[k] = -g_z_z_yyy_zz[k] * ab_x + g_z_z_yyy_xzz[k];
            }

            /// Set up 762-768 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyyz_xx = cbuffer.data(gd_geom_11_off + 762 * ccomps * dcomps);

            auto g_z_z_xyyz_xy = cbuffer.data(gd_geom_11_off + 763 * ccomps * dcomps);

            auto g_z_z_xyyz_xz = cbuffer.data(gd_geom_11_off + 764 * ccomps * dcomps);

            auto g_z_z_xyyz_yy = cbuffer.data(gd_geom_11_off + 765 * ccomps * dcomps);

            auto g_z_z_xyyz_yz = cbuffer.data(gd_geom_11_off + 766 * ccomps * dcomps);

            auto g_z_z_xyyz_zz = cbuffer.data(gd_geom_11_off + 767 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyyz_xx, g_z_z_xyyz_xy, g_z_z_xyyz_xz, g_z_z_xyyz_yy, g_z_z_xyyz_yz, g_z_z_xyyz_zz, g_z_z_yyz_xx, g_z_z_yyz_xxx, g_z_z_yyz_xxy, g_z_z_yyz_xxz, g_z_z_yyz_xy, g_z_z_yyz_xyy, g_z_z_yyz_xyz, g_z_z_yyz_xz, g_z_z_yyz_xzz, g_z_z_yyz_yy, g_z_z_yyz_yz, g_z_z_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyyz_xx[k] = -g_z_z_yyz_xx[k] * ab_x + g_z_z_yyz_xxx[k];

                g_z_z_xyyz_xy[k] = -g_z_z_yyz_xy[k] * ab_x + g_z_z_yyz_xxy[k];

                g_z_z_xyyz_xz[k] = -g_z_z_yyz_xz[k] * ab_x + g_z_z_yyz_xxz[k];

                g_z_z_xyyz_yy[k] = -g_z_z_yyz_yy[k] * ab_x + g_z_z_yyz_xyy[k];

                g_z_z_xyyz_yz[k] = -g_z_z_yyz_yz[k] * ab_x + g_z_z_yyz_xyz[k];

                g_z_z_xyyz_zz[k] = -g_z_z_yyz_zz[k] * ab_x + g_z_z_yyz_xzz[k];
            }

            /// Set up 768-774 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyzz_xx = cbuffer.data(gd_geom_11_off + 768 * ccomps * dcomps);

            auto g_z_z_xyzz_xy = cbuffer.data(gd_geom_11_off + 769 * ccomps * dcomps);

            auto g_z_z_xyzz_xz = cbuffer.data(gd_geom_11_off + 770 * ccomps * dcomps);

            auto g_z_z_xyzz_yy = cbuffer.data(gd_geom_11_off + 771 * ccomps * dcomps);

            auto g_z_z_xyzz_yz = cbuffer.data(gd_geom_11_off + 772 * ccomps * dcomps);

            auto g_z_z_xyzz_zz = cbuffer.data(gd_geom_11_off + 773 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyzz_xx, g_z_z_xyzz_xy, g_z_z_xyzz_xz, g_z_z_xyzz_yy, g_z_z_xyzz_yz, g_z_z_xyzz_zz, g_z_z_yzz_xx, g_z_z_yzz_xxx, g_z_z_yzz_xxy, g_z_z_yzz_xxz, g_z_z_yzz_xy, g_z_z_yzz_xyy, g_z_z_yzz_xyz, g_z_z_yzz_xz, g_z_z_yzz_xzz, g_z_z_yzz_yy, g_z_z_yzz_yz, g_z_z_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyzz_xx[k] = -g_z_z_yzz_xx[k] * ab_x + g_z_z_yzz_xxx[k];

                g_z_z_xyzz_xy[k] = -g_z_z_yzz_xy[k] * ab_x + g_z_z_yzz_xxy[k];

                g_z_z_xyzz_xz[k] = -g_z_z_yzz_xz[k] * ab_x + g_z_z_yzz_xxz[k];

                g_z_z_xyzz_yy[k] = -g_z_z_yzz_yy[k] * ab_x + g_z_z_yzz_xyy[k];

                g_z_z_xyzz_yz[k] = -g_z_z_yzz_yz[k] * ab_x + g_z_z_yzz_xyz[k];

                g_z_z_xyzz_zz[k] = -g_z_z_yzz_zz[k] * ab_x + g_z_z_yzz_xzz[k];
            }

            /// Set up 774-780 components of targeted buffer : cbuffer.data(

            auto g_z_z_xzzz_xx = cbuffer.data(gd_geom_11_off + 774 * ccomps * dcomps);

            auto g_z_z_xzzz_xy = cbuffer.data(gd_geom_11_off + 775 * ccomps * dcomps);

            auto g_z_z_xzzz_xz = cbuffer.data(gd_geom_11_off + 776 * ccomps * dcomps);

            auto g_z_z_xzzz_yy = cbuffer.data(gd_geom_11_off + 777 * ccomps * dcomps);

            auto g_z_z_xzzz_yz = cbuffer.data(gd_geom_11_off + 778 * ccomps * dcomps);

            auto g_z_z_xzzz_zz = cbuffer.data(gd_geom_11_off + 779 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xzzz_xx, g_z_z_xzzz_xy, g_z_z_xzzz_xz, g_z_z_xzzz_yy, g_z_z_xzzz_yz, g_z_z_xzzz_zz, g_z_z_zzz_xx, g_z_z_zzz_xxx, g_z_z_zzz_xxy, g_z_z_zzz_xxz, g_z_z_zzz_xy, g_z_z_zzz_xyy, g_z_z_zzz_xyz, g_z_z_zzz_xz, g_z_z_zzz_xzz, g_z_z_zzz_yy, g_z_z_zzz_yz, g_z_z_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xzzz_xx[k] = -g_z_z_zzz_xx[k] * ab_x + g_z_z_zzz_xxx[k];

                g_z_z_xzzz_xy[k] = -g_z_z_zzz_xy[k] * ab_x + g_z_z_zzz_xxy[k];

                g_z_z_xzzz_xz[k] = -g_z_z_zzz_xz[k] * ab_x + g_z_z_zzz_xxz[k];

                g_z_z_xzzz_yy[k] = -g_z_z_zzz_yy[k] * ab_x + g_z_z_zzz_xyy[k];

                g_z_z_xzzz_yz[k] = -g_z_z_zzz_yz[k] * ab_x + g_z_z_zzz_xyz[k];

                g_z_z_xzzz_zz[k] = -g_z_z_zzz_zz[k] * ab_x + g_z_z_zzz_xzz[k];
            }

            /// Set up 780-786 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyyy_xx = cbuffer.data(gd_geom_11_off + 780 * ccomps * dcomps);

            auto g_z_z_yyyy_xy = cbuffer.data(gd_geom_11_off + 781 * ccomps * dcomps);

            auto g_z_z_yyyy_xz = cbuffer.data(gd_geom_11_off + 782 * ccomps * dcomps);

            auto g_z_z_yyyy_yy = cbuffer.data(gd_geom_11_off + 783 * ccomps * dcomps);

            auto g_z_z_yyyy_yz = cbuffer.data(gd_geom_11_off + 784 * ccomps * dcomps);

            auto g_z_z_yyyy_zz = cbuffer.data(gd_geom_11_off + 785 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yyy_xx, g_z_z_yyy_xxy, g_z_z_yyy_xy, g_z_z_yyy_xyy, g_z_z_yyy_xyz, g_z_z_yyy_xz, g_z_z_yyy_yy, g_z_z_yyy_yyy, g_z_z_yyy_yyz, g_z_z_yyy_yz, g_z_z_yyy_yzz, g_z_z_yyy_zz, g_z_z_yyyy_xx, g_z_z_yyyy_xy, g_z_z_yyyy_xz, g_z_z_yyyy_yy, g_z_z_yyyy_yz, g_z_z_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyyy_xx[k] = -g_z_z_yyy_xx[k] * ab_y + g_z_z_yyy_xxy[k];

                g_z_z_yyyy_xy[k] = -g_z_z_yyy_xy[k] * ab_y + g_z_z_yyy_xyy[k];

                g_z_z_yyyy_xz[k] = -g_z_z_yyy_xz[k] * ab_y + g_z_z_yyy_xyz[k];

                g_z_z_yyyy_yy[k] = -g_z_z_yyy_yy[k] * ab_y + g_z_z_yyy_yyy[k];

                g_z_z_yyyy_yz[k] = -g_z_z_yyy_yz[k] * ab_y + g_z_z_yyy_yyz[k];

                g_z_z_yyyy_zz[k] = -g_z_z_yyy_zz[k] * ab_y + g_z_z_yyy_yzz[k];
            }

            /// Set up 786-792 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyyz_xx = cbuffer.data(gd_geom_11_off + 786 * ccomps * dcomps);

            auto g_z_z_yyyz_xy = cbuffer.data(gd_geom_11_off + 787 * ccomps * dcomps);

            auto g_z_z_yyyz_xz = cbuffer.data(gd_geom_11_off + 788 * ccomps * dcomps);

            auto g_z_z_yyyz_yy = cbuffer.data(gd_geom_11_off + 789 * ccomps * dcomps);

            auto g_z_z_yyyz_yz = cbuffer.data(gd_geom_11_off + 790 * ccomps * dcomps);

            auto g_z_z_yyyz_zz = cbuffer.data(gd_geom_11_off + 791 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yyyz_xx, g_z_z_yyyz_xy, g_z_z_yyyz_xz, g_z_z_yyyz_yy, g_z_z_yyyz_yz, g_z_z_yyyz_zz, g_z_z_yyz_xx, g_z_z_yyz_xxy, g_z_z_yyz_xy, g_z_z_yyz_xyy, g_z_z_yyz_xyz, g_z_z_yyz_xz, g_z_z_yyz_yy, g_z_z_yyz_yyy, g_z_z_yyz_yyz, g_z_z_yyz_yz, g_z_z_yyz_yzz, g_z_z_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyyz_xx[k] = -g_z_z_yyz_xx[k] * ab_y + g_z_z_yyz_xxy[k];

                g_z_z_yyyz_xy[k] = -g_z_z_yyz_xy[k] * ab_y + g_z_z_yyz_xyy[k];

                g_z_z_yyyz_xz[k] = -g_z_z_yyz_xz[k] * ab_y + g_z_z_yyz_xyz[k];

                g_z_z_yyyz_yy[k] = -g_z_z_yyz_yy[k] * ab_y + g_z_z_yyz_yyy[k];

                g_z_z_yyyz_yz[k] = -g_z_z_yyz_yz[k] * ab_y + g_z_z_yyz_yyz[k];

                g_z_z_yyyz_zz[k] = -g_z_z_yyz_zz[k] * ab_y + g_z_z_yyz_yzz[k];
            }

            /// Set up 792-798 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyzz_xx = cbuffer.data(gd_geom_11_off + 792 * ccomps * dcomps);

            auto g_z_z_yyzz_xy = cbuffer.data(gd_geom_11_off + 793 * ccomps * dcomps);

            auto g_z_z_yyzz_xz = cbuffer.data(gd_geom_11_off + 794 * ccomps * dcomps);

            auto g_z_z_yyzz_yy = cbuffer.data(gd_geom_11_off + 795 * ccomps * dcomps);

            auto g_z_z_yyzz_yz = cbuffer.data(gd_geom_11_off + 796 * ccomps * dcomps);

            auto g_z_z_yyzz_zz = cbuffer.data(gd_geom_11_off + 797 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yyzz_xx, g_z_z_yyzz_xy, g_z_z_yyzz_xz, g_z_z_yyzz_yy, g_z_z_yyzz_yz, g_z_z_yyzz_zz, g_z_z_yzz_xx, g_z_z_yzz_xxy, g_z_z_yzz_xy, g_z_z_yzz_xyy, g_z_z_yzz_xyz, g_z_z_yzz_xz, g_z_z_yzz_yy, g_z_z_yzz_yyy, g_z_z_yzz_yyz, g_z_z_yzz_yz, g_z_z_yzz_yzz, g_z_z_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyzz_xx[k] = -g_z_z_yzz_xx[k] * ab_y + g_z_z_yzz_xxy[k];

                g_z_z_yyzz_xy[k] = -g_z_z_yzz_xy[k] * ab_y + g_z_z_yzz_xyy[k];

                g_z_z_yyzz_xz[k] = -g_z_z_yzz_xz[k] * ab_y + g_z_z_yzz_xyz[k];

                g_z_z_yyzz_yy[k] = -g_z_z_yzz_yy[k] * ab_y + g_z_z_yzz_yyy[k];

                g_z_z_yyzz_yz[k] = -g_z_z_yzz_yz[k] * ab_y + g_z_z_yzz_yyz[k];

                g_z_z_yyzz_zz[k] = -g_z_z_yzz_zz[k] * ab_y + g_z_z_yzz_yzz[k];
            }

            /// Set up 798-804 components of targeted buffer : cbuffer.data(

            auto g_z_z_yzzz_xx = cbuffer.data(gd_geom_11_off + 798 * ccomps * dcomps);

            auto g_z_z_yzzz_xy = cbuffer.data(gd_geom_11_off + 799 * ccomps * dcomps);

            auto g_z_z_yzzz_xz = cbuffer.data(gd_geom_11_off + 800 * ccomps * dcomps);

            auto g_z_z_yzzz_yy = cbuffer.data(gd_geom_11_off + 801 * ccomps * dcomps);

            auto g_z_z_yzzz_yz = cbuffer.data(gd_geom_11_off + 802 * ccomps * dcomps);

            auto g_z_z_yzzz_zz = cbuffer.data(gd_geom_11_off + 803 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yzzz_xx, g_z_z_yzzz_xy, g_z_z_yzzz_xz, g_z_z_yzzz_yy, g_z_z_yzzz_yz, g_z_z_yzzz_zz, g_z_z_zzz_xx, g_z_z_zzz_xxy, g_z_z_zzz_xy, g_z_z_zzz_xyy, g_z_z_zzz_xyz, g_z_z_zzz_xz, g_z_z_zzz_yy, g_z_z_zzz_yyy, g_z_z_zzz_yyz, g_z_z_zzz_yz, g_z_z_zzz_yzz, g_z_z_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yzzz_xx[k] = -g_z_z_zzz_xx[k] * ab_y + g_z_z_zzz_xxy[k];

                g_z_z_yzzz_xy[k] = -g_z_z_zzz_xy[k] * ab_y + g_z_z_zzz_xyy[k];

                g_z_z_yzzz_xz[k] = -g_z_z_zzz_xz[k] * ab_y + g_z_z_zzz_xyz[k];

                g_z_z_yzzz_yy[k] = -g_z_z_zzz_yy[k] * ab_y + g_z_z_zzz_yyy[k];

                g_z_z_yzzz_yz[k] = -g_z_z_zzz_yz[k] * ab_y + g_z_z_zzz_yyz[k];

                g_z_z_yzzz_zz[k] = -g_z_z_zzz_zz[k] * ab_y + g_z_z_zzz_yzz[k];
            }

            /// Set up 804-810 components of targeted buffer : cbuffer.data(

            auto g_z_z_zzzz_xx = cbuffer.data(gd_geom_11_off + 804 * ccomps * dcomps);

            auto g_z_z_zzzz_xy = cbuffer.data(gd_geom_11_off + 805 * ccomps * dcomps);

            auto g_z_z_zzzz_xz = cbuffer.data(gd_geom_11_off + 806 * ccomps * dcomps);

            auto g_z_z_zzzz_yy = cbuffer.data(gd_geom_11_off + 807 * ccomps * dcomps);

            auto g_z_z_zzzz_yz = cbuffer.data(gd_geom_11_off + 808 * ccomps * dcomps);

            auto g_z_z_zzzz_zz = cbuffer.data(gd_geom_11_off + 809 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzz_xx, g_0_z_zzz_xy, g_0_z_zzz_xz, g_0_z_zzz_yy, g_0_z_zzz_yz, g_0_z_zzz_zz, g_z_0_zzz_xx, g_z_0_zzz_xy, g_z_0_zzz_xz, g_z_0_zzz_yy, g_z_0_zzz_yz, g_z_0_zzz_zz, g_z_z_zzz_xx, g_z_z_zzz_xxz, g_z_z_zzz_xy, g_z_z_zzz_xyz, g_z_z_zzz_xz, g_z_z_zzz_xzz, g_z_z_zzz_yy, g_z_z_zzz_yyz, g_z_z_zzz_yz, g_z_z_zzz_yzz, g_z_z_zzz_zz, g_z_z_zzz_zzz, g_z_z_zzzz_xx, g_z_z_zzzz_xy, g_z_z_zzzz_xz, g_z_z_zzzz_yy, g_z_z_zzzz_yz, g_z_z_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_zzzz_xx[k] = -g_0_z_zzz_xx[k] + g_z_0_zzz_xx[k] - g_z_z_zzz_xx[k] * ab_z + g_z_z_zzz_xxz[k];

                g_z_z_zzzz_xy[k] = -g_0_z_zzz_xy[k] + g_z_0_zzz_xy[k] - g_z_z_zzz_xy[k] * ab_z + g_z_z_zzz_xyz[k];

                g_z_z_zzzz_xz[k] = -g_0_z_zzz_xz[k] + g_z_0_zzz_xz[k] - g_z_z_zzz_xz[k] * ab_z + g_z_z_zzz_xzz[k];

                g_z_z_zzzz_yy[k] = -g_0_z_zzz_yy[k] + g_z_0_zzz_yy[k] - g_z_z_zzz_yy[k] * ab_z + g_z_z_zzz_yyz[k];

                g_z_z_zzzz_yz[k] = -g_0_z_zzz_yz[k] + g_z_0_zzz_yz[k] - g_z_z_zzz_yz[k] * ab_z + g_z_z_zzz_yzz[k];

                g_z_z_zzzz_zz[k] = -g_0_z_zzz_zz[k] + g_z_0_zzz_zz[k] - g_z_z_zzz_zz[k] * ab_z + g_z_z_zzz_zzz[k];
            }
        }
    }
}

} // erirec namespace

