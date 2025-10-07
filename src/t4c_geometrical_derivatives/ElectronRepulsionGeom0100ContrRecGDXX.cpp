#include "ElectronRepulsionGeom0100ContrRecGDXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_gdxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_gdxx,
                                            const size_t idx_fdxx,
                                            const size_t idx_geom_01_fdxx,
                                            const size_t idx_geom_01_ffxx,
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

            const auto fd_off = idx_fdxx + i * dcomps + j;

            auto g_xxx_xx = cbuffer.data(fd_off + 0 * ccomps * dcomps);

            auto g_xxx_xy = cbuffer.data(fd_off + 1 * ccomps * dcomps);

            auto g_xxx_xz = cbuffer.data(fd_off + 2 * ccomps * dcomps);

            auto g_xxx_yy = cbuffer.data(fd_off + 3 * ccomps * dcomps);

            auto g_xxx_yz = cbuffer.data(fd_off + 4 * ccomps * dcomps);

            auto g_xxx_zz = cbuffer.data(fd_off + 5 * ccomps * dcomps);

            auto g_xxy_xx = cbuffer.data(fd_off + 6 * ccomps * dcomps);

            auto g_xxy_xy = cbuffer.data(fd_off + 7 * ccomps * dcomps);

            auto g_xxy_xz = cbuffer.data(fd_off + 8 * ccomps * dcomps);

            auto g_xxy_yy = cbuffer.data(fd_off + 9 * ccomps * dcomps);

            auto g_xxy_yz = cbuffer.data(fd_off + 10 * ccomps * dcomps);

            auto g_xxy_zz = cbuffer.data(fd_off + 11 * ccomps * dcomps);

            auto g_xxz_xx = cbuffer.data(fd_off + 12 * ccomps * dcomps);

            auto g_xxz_xy = cbuffer.data(fd_off + 13 * ccomps * dcomps);

            auto g_xxz_xz = cbuffer.data(fd_off + 14 * ccomps * dcomps);

            auto g_xxz_yy = cbuffer.data(fd_off + 15 * ccomps * dcomps);

            auto g_xxz_yz = cbuffer.data(fd_off + 16 * ccomps * dcomps);

            auto g_xxz_zz = cbuffer.data(fd_off + 17 * ccomps * dcomps);

            auto g_xyy_xx = cbuffer.data(fd_off + 18 * ccomps * dcomps);

            auto g_xyy_xy = cbuffer.data(fd_off + 19 * ccomps * dcomps);

            auto g_xyy_xz = cbuffer.data(fd_off + 20 * ccomps * dcomps);

            auto g_xyy_yy = cbuffer.data(fd_off + 21 * ccomps * dcomps);

            auto g_xyy_yz = cbuffer.data(fd_off + 22 * ccomps * dcomps);

            auto g_xyy_zz = cbuffer.data(fd_off + 23 * ccomps * dcomps);

            auto g_xyz_xx = cbuffer.data(fd_off + 24 * ccomps * dcomps);

            auto g_xyz_xy = cbuffer.data(fd_off + 25 * ccomps * dcomps);

            auto g_xyz_xz = cbuffer.data(fd_off + 26 * ccomps * dcomps);

            auto g_xyz_yy = cbuffer.data(fd_off + 27 * ccomps * dcomps);

            auto g_xyz_yz = cbuffer.data(fd_off + 28 * ccomps * dcomps);

            auto g_xyz_zz = cbuffer.data(fd_off + 29 * ccomps * dcomps);

            auto g_xzz_xx = cbuffer.data(fd_off + 30 * ccomps * dcomps);

            auto g_xzz_xy = cbuffer.data(fd_off + 31 * ccomps * dcomps);

            auto g_xzz_xz = cbuffer.data(fd_off + 32 * ccomps * dcomps);

            auto g_xzz_yy = cbuffer.data(fd_off + 33 * ccomps * dcomps);

            auto g_xzz_yz = cbuffer.data(fd_off + 34 * ccomps * dcomps);

            auto g_xzz_zz = cbuffer.data(fd_off + 35 * ccomps * dcomps);

            auto g_yyy_xx = cbuffer.data(fd_off + 36 * ccomps * dcomps);

            auto g_yyy_xy = cbuffer.data(fd_off + 37 * ccomps * dcomps);

            auto g_yyy_xz = cbuffer.data(fd_off + 38 * ccomps * dcomps);

            auto g_yyy_yy = cbuffer.data(fd_off + 39 * ccomps * dcomps);

            auto g_yyy_yz = cbuffer.data(fd_off + 40 * ccomps * dcomps);

            auto g_yyy_zz = cbuffer.data(fd_off + 41 * ccomps * dcomps);

            auto g_yyz_xx = cbuffer.data(fd_off + 42 * ccomps * dcomps);

            auto g_yyz_xy = cbuffer.data(fd_off + 43 * ccomps * dcomps);

            auto g_yyz_xz = cbuffer.data(fd_off + 44 * ccomps * dcomps);

            auto g_yyz_yy = cbuffer.data(fd_off + 45 * ccomps * dcomps);

            auto g_yyz_yz = cbuffer.data(fd_off + 46 * ccomps * dcomps);

            auto g_yyz_zz = cbuffer.data(fd_off + 47 * ccomps * dcomps);

            auto g_yzz_xx = cbuffer.data(fd_off + 48 * ccomps * dcomps);

            auto g_yzz_xy = cbuffer.data(fd_off + 49 * ccomps * dcomps);

            auto g_yzz_xz = cbuffer.data(fd_off + 50 * ccomps * dcomps);

            auto g_yzz_yy = cbuffer.data(fd_off + 51 * ccomps * dcomps);

            auto g_yzz_yz = cbuffer.data(fd_off + 52 * ccomps * dcomps);

            auto g_yzz_zz = cbuffer.data(fd_off + 53 * ccomps * dcomps);

            auto g_zzz_xx = cbuffer.data(fd_off + 54 * ccomps * dcomps);

            auto g_zzz_xy = cbuffer.data(fd_off + 55 * ccomps * dcomps);

            auto g_zzz_xz = cbuffer.data(fd_off + 56 * ccomps * dcomps);

            auto g_zzz_yy = cbuffer.data(fd_off + 57 * ccomps * dcomps);

            auto g_zzz_yz = cbuffer.data(fd_off + 58 * ccomps * dcomps);

            auto g_zzz_zz = cbuffer.data(fd_off + 59 * ccomps * dcomps);

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

            /// Set up components of auxilary buffer : FFSS

            const auto ff_geom_01_off = idx_geom_01_ffxx + i * dcomps + j;

            auto g_0_x_xxx_xxx = cbuffer.data(ff_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxx_xxy = cbuffer.data(ff_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxx_xxz = cbuffer.data(ff_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxx_xyy = cbuffer.data(ff_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxx_xyz = cbuffer.data(ff_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxx_xzz = cbuffer.data(ff_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxx_yyy = cbuffer.data(ff_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxx_yyz = cbuffer.data(ff_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxx_yzz = cbuffer.data(ff_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxx_zzz = cbuffer.data(ff_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxy_xxx = cbuffer.data(ff_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxy_xxy = cbuffer.data(ff_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxy_xxz = cbuffer.data(ff_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxy_xyy = cbuffer.data(ff_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxy_xyz = cbuffer.data(ff_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxy_xzz = cbuffer.data(ff_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxy_yyy = cbuffer.data(ff_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxy_yyz = cbuffer.data(ff_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxy_yzz = cbuffer.data(ff_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxy_zzz = cbuffer.data(ff_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxz_xxx = cbuffer.data(ff_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxz_xxy = cbuffer.data(ff_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxz_xxz = cbuffer.data(ff_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxz_xyy = cbuffer.data(ff_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxz_xyz = cbuffer.data(ff_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxz_xzz = cbuffer.data(ff_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxz_yyy = cbuffer.data(ff_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxz_yyz = cbuffer.data(ff_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxz_yzz = cbuffer.data(ff_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxz_zzz = cbuffer.data(ff_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xyy_xxx = cbuffer.data(ff_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xyy_xxy = cbuffer.data(ff_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xyy_xxz = cbuffer.data(ff_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xyy_xyy = cbuffer.data(ff_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xyy_xyz = cbuffer.data(ff_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xyy_xzz = cbuffer.data(ff_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xyy_yyy = cbuffer.data(ff_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xyy_yyz = cbuffer.data(ff_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xyy_yzz = cbuffer.data(ff_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xyy_zzz = cbuffer.data(ff_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xyz_xxx = cbuffer.data(ff_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xyz_xxy = cbuffer.data(ff_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xyz_xxz = cbuffer.data(ff_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xyz_xyy = cbuffer.data(ff_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xyz_xyz = cbuffer.data(ff_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xyz_xzz = cbuffer.data(ff_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xyz_yyy = cbuffer.data(ff_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xyz_yyz = cbuffer.data(ff_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xyz_yzz = cbuffer.data(ff_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xyz_zzz = cbuffer.data(ff_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xzz_xxx = cbuffer.data(ff_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xzz_xxy = cbuffer.data(ff_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xzz_xxz = cbuffer.data(ff_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xzz_xyy = cbuffer.data(ff_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xzz_xyz = cbuffer.data(ff_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xzz_xzz = cbuffer.data(ff_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xzz_yyy = cbuffer.data(ff_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xzz_yyz = cbuffer.data(ff_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xzz_yzz = cbuffer.data(ff_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xzz_zzz = cbuffer.data(ff_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_yyy_xxx = cbuffer.data(ff_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_yyy_xxy = cbuffer.data(ff_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_yyy_xxz = cbuffer.data(ff_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_yyy_xyy = cbuffer.data(ff_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_yyy_xyz = cbuffer.data(ff_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_yyy_xzz = cbuffer.data(ff_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_yyy_yyy = cbuffer.data(ff_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_yyy_yyz = cbuffer.data(ff_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_yyy_yzz = cbuffer.data(ff_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_yyy_zzz = cbuffer.data(ff_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_yyz_xxx = cbuffer.data(ff_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_yyz_xxy = cbuffer.data(ff_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_yyz_xxz = cbuffer.data(ff_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_yyz_xyy = cbuffer.data(ff_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_yyz_xyz = cbuffer.data(ff_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_yyz_xzz = cbuffer.data(ff_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_yyz_yyy = cbuffer.data(ff_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_yyz_yyz = cbuffer.data(ff_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_yyz_yzz = cbuffer.data(ff_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_yyz_zzz = cbuffer.data(ff_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_yzz_xxx = cbuffer.data(ff_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_yzz_xxy = cbuffer.data(ff_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_yzz_xxz = cbuffer.data(ff_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_yzz_xyy = cbuffer.data(ff_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_yzz_xyz = cbuffer.data(ff_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_yzz_xzz = cbuffer.data(ff_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_yzz_yyy = cbuffer.data(ff_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_yzz_yyz = cbuffer.data(ff_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_yzz_yzz = cbuffer.data(ff_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_yzz_zzz = cbuffer.data(ff_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_zzz_xxx = cbuffer.data(ff_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_zzz_xxy = cbuffer.data(ff_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_zzz_xxz = cbuffer.data(ff_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_zzz_xyy = cbuffer.data(ff_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_zzz_xyz = cbuffer.data(ff_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_zzz_xzz = cbuffer.data(ff_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_zzz_yyy = cbuffer.data(ff_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_zzz_yyz = cbuffer.data(ff_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_zzz_yzz = cbuffer.data(ff_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_zzz_zzz = cbuffer.data(ff_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_y_xxx_xxx = cbuffer.data(ff_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_y_xxx_xxy = cbuffer.data(ff_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_y_xxx_xxz = cbuffer.data(ff_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_y_xxx_xyy = cbuffer.data(ff_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_y_xxx_xyz = cbuffer.data(ff_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_y_xxx_xzz = cbuffer.data(ff_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_y_xxx_yyy = cbuffer.data(ff_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_y_xxx_yyz = cbuffer.data(ff_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_y_xxx_yzz = cbuffer.data(ff_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_xxx_zzz = cbuffer.data(ff_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_xxy_xxx = cbuffer.data(ff_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_y_xxy_xxy = cbuffer.data(ff_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_y_xxy_xxz = cbuffer.data(ff_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_xxy_xyy = cbuffer.data(ff_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_y_xxy_xyz = cbuffer.data(ff_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_xxy_xzz = cbuffer.data(ff_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_xxy_yyy = cbuffer.data(ff_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_y_xxy_yyz = cbuffer.data(ff_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_xxy_yzz = cbuffer.data(ff_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_xxy_zzz = cbuffer.data(ff_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_y_xxz_xxx = cbuffer.data(ff_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_y_xxz_xxy = cbuffer.data(ff_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_y_xxz_xxz = cbuffer.data(ff_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_y_xxz_xyy = cbuffer.data(ff_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_y_xxz_xyz = cbuffer.data(ff_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_y_xxz_xzz = cbuffer.data(ff_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_y_xxz_yyy = cbuffer.data(ff_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_y_xxz_yyz = cbuffer.data(ff_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_y_xxz_yzz = cbuffer.data(ff_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_y_xxz_zzz = cbuffer.data(ff_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_y_xyy_xxx = cbuffer.data(ff_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_y_xyy_xxy = cbuffer.data(ff_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_y_xyy_xxz = cbuffer.data(ff_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_y_xyy_xyy = cbuffer.data(ff_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_y_xyy_xyz = cbuffer.data(ff_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_y_xyy_xzz = cbuffer.data(ff_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_xyy_yyy = cbuffer.data(ff_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_xyy_yyz = cbuffer.data(ff_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_y_xyy_yzz = cbuffer.data(ff_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_xyy_zzz = cbuffer.data(ff_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_y_xyz_xxx = cbuffer.data(ff_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_y_xyz_xxy = cbuffer.data(ff_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_xyz_xxz = cbuffer.data(ff_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_xyz_xyy = cbuffer.data(ff_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_y_xyz_xyz = cbuffer.data(ff_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_xyz_xzz = cbuffer.data(ff_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_xyz_yyy = cbuffer.data(ff_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_y_xyz_yyz = cbuffer.data(ff_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_xyz_yzz = cbuffer.data(ff_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_xyz_zzz = cbuffer.data(ff_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_y_xzz_xxx = cbuffer.data(ff_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_xzz_xxy = cbuffer.data(ff_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_xzz_xxz = cbuffer.data(ff_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_xzz_xyy = cbuffer.data(ff_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_xzz_xyz = cbuffer.data(ff_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_xzz_xzz = cbuffer.data(ff_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_xzz_yyy = cbuffer.data(ff_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_xzz_yyz = cbuffer.data(ff_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_xzz_yzz = cbuffer.data(ff_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_xzz_zzz = cbuffer.data(ff_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_yyy_xxx = cbuffer.data(ff_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_yyy_xxy = cbuffer.data(ff_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_yyy_xxz = cbuffer.data(ff_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_yyy_xyy = cbuffer.data(ff_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_yyy_xyz = cbuffer.data(ff_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_yyy_xzz = cbuffer.data(ff_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_yyy_yyy = cbuffer.data(ff_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_yyy_yyz = cbuffer.data(ff_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_y_yyy_yzz = cbuffer.data(ff_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_yyy_zzz = cbuffer.data(ff_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_yyz_xxx = cbuffer.data(ff_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_yyz_xxy = cbuffer.data(ff_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_yyz_xxz = cbuffer.data(ff_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_yyz_xyy = cbuffer.data(ff_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_yyz_xyz = cbuffer.data(ff_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_yyz_xzz = cbuffer.data(ff_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_yyz_yyy = cbuffer.data(ff_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_yyz_yyz = cbuffer.data(ff_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_yyz_yzz = cbuffer.data(ff_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_yyz_zzz = cbuffer.data(ff_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_y_yzz_xxx = cbuffer.data(ff_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_y_yzz_xxy = cbuffer.data(ff_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_y_yzz_xxz = cbuffer.data(ff_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_y_yzz_xyy = cbuffer.data(ff_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_y_yzz_xyz = cbuffer.data(ff_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_y_yzz_xzz = cbuffer.data(ff_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_y_yzz_yyy = cbuffer.data(ff_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_y_yzz_yyz = cbuffer.data(ff_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_y_yzz_yzz = cbuffer.data(ff_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_y_yzz_zzz = cbuffer.data(ff_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_y_zzz_xxx = cbuffer.data(ff_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_y_zzz_xxy = cbuffer.data(ff_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_y_zzz_xxz = cbuffer.data(ff_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_y_zzz_xyy = cbuffer.data(ff_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_y_zzz_xyz = cbuffer.data(ff_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_y_zzz_xzz = cbuffer.data(ff_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_y_zzz_yyy = cbuffer.data(ff_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_y_zzz_yyz = cbuffer.data(ff_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_y_zzz_yzz = cbuffer.data(ff_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_y_zzz_zzz = cbuffer.data(ff_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_z_xxx_xxx = cbuffer.data(ff_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_z_xxx_xxy = cbuffer.data(ff_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_z_xxx_xxz = cbuffer.data(ff_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_z_xxx_xyy = cbuffer.data(ff_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_z_xxx_xyz = cbuffer.data(ff_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_z_xxx_xzz = cbuffer.data(ff_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_z_xxx_yyy = cbuffer.data(ff_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_z_xxx_yyz = cbuffer.data(ff_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_z_xxx_yzz = cbuffer.data(ff_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_z_xxx_zzz = cbuffer.data(ff_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_z_xxy_xxx = cbuffer.data(ff_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_z_xxy_xxy = cbuffer.data(ff_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_z_xxy_xxz = cbuffer.data(ff_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_z_xxy_xyy = cbuffer.data(ff_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_z_xxy_xyz = cbuffer.data(ff_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_z_xxy_xzz = cbuffer.data(ff_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_z_xxy_yyy = cbuffer.data(ff_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_z_xxy_yyz = cbuffer.data(ff_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_z_xxy_yzz = cbuffer.data(ff_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_z_xxy_zzz = cbuffer.data(ff_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_z_xxz_xxx = cbuffer.data(ff_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_z_xxz_xxy = cbuffer.data(ff_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_z_xxz_xxz = cbuffer.data(ff_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_z_xxz_xyy = cbuffer.data(ff_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_z_xxz_xyz = cbuffer.data(ff_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_z_xxz_xzz = cbuffer.data(ff_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_z_xxz_yyy = cbuffer.data(ff_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_z_xxz_yyz = cbuffer.data(ff_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_z_xxz_yzz = cbuffer.data(ff_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_z_xxz_zzz = cbuffer.data(ff_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_z_xyy_xxx = cbuffer.data(ff_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_z_xyy_xxy = cbuffer.data(ff_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_z_xyy_xxz = cbuffer.data(ff_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_z_xyy_xyy = cbuffer.data(ff_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_z_xyy_xyz = cbuffer.data(ff_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_z_xyy_xzz = cbuffer.data(ff_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_z_xyy_yyy = cbuffer.data(ff_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_z_xyy_yyz = cbuffer.data(ff_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_z_xyy_yzz = cbuffer.data(ff_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_z_xyy_zzz = cbuffer.data(ff_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_z_xyz_xxx = cbuffer.data(ff_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_z_xyz_xxy = cbuffer.data(ff_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_z_xyz_xxz = cbuffer.data(ff_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_z_xyz_xyy = cbuffer.data(ff_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_z_xyz_xyz = cbuffer.data(ff_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_z_xyz_xzz = cbuffer.data(ff_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_z_xyz_yyy = cbuffer.data(ff_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_z_xyz_yyz = cbuffer.data(ff_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_z_xyz_yzz = cbuffer.data(ff_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_z_xyz_zzz = cbuffer.data(ff_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_z_xzz_xxx = cbuffer.data(ff_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_z_xzz_xxy = cbuffer.data(ff_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_z_xzz_xxz = cbuffer.data(ff_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_z_xzz_xyy = cbuffer.data(ff_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_z_xzz_xyz = cbuffer.data(ff_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_z_xzz_xzz = cbuffer.data(ff_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_z_xzz_yyy = cbuffer.data(ff_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_z_xzz_yyz = cbuffer.data(ff_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_z_xzz_yzz = cbuffer.data(ff_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_z_xzz_zzz = cbuffer.data(ff_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_z_yyy_xxx = cbuffer.data(ff_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_z_yyy_xxy = cbuffer.data(ff_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_z_yyy_xxz = cbuffer.data(ff_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_z_yyy_xyy = cbuffer.data(ff_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_z_yyy_xyz = cbuffer.data(ff_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_z_yyy_xzz = cbuffer.data(ff_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_z_yyy_yyy = cbuffer.data(ff_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_z_yyy_yyz = cbuffer.data(ff_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_z_yyy_yzz = cbuffer.data(ff_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_z_yyy_zzz = cbuffer.data(ff_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_z_yyz_xxx = cbuffer.data(ff_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_z_yyz_xxy = cbuffer.data(ff_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_z_yyz_xxz = cbuffer.data(ff_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_z_yyz_xyy = cbuffer.data(ff_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_z_yyz_xyz = cbuffer.data(ff_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_z_yyz_xzz = cbuffer.data(ff_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_z_yyz_yyy = cbuffer.data(ff_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_z_yyz_yyz = cbuffer.data(ff_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_z_yyz_yzz = cbuffer.data(ff_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_z_yyz_zzz = cbuffer.data(ff_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_z_yzz_xxx = cbuffer.data(ff_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_z_yzz_xxy = cbuffer.data(ff_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_z_yzz_xxz = cbuffer.data(ff_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_z_yzz_xyy = cbuffer.data(ff_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_z_yzz_xyz = cbuffer.data(ff_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_z_yzz_xzz = cbuffer.data(ff_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_z_yzz_yyy = cbuffer.data(ff_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_z_yzz_yyz = cbuffer.data(ff_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_z_yzz_yzz = cbuffer.data(ff_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_z_yzz_zzz = cbuffer.data(ff_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_z_zzz_xxx = cbuffer.data(ff_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_z_zzz_xxy = cbuffer.data(ff_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_z_zzz_xxz = cbuffer.data(ff_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_z_zzz_xyy = cbuffer.data(ff_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_z_zzz_xyz = cbuffer.data(ff_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_z_zzz_xzz = cbuffer.data(ff_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_z_zzz_yyy = cbuffer.data(ff_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_z_zzz_yyz = cbuffer.data(ff_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_z_zzz_yzz = cbuffer.data(ff_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_z_zzz_zzz = cbuffer.data(ff_geom_01_off + 299 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_gdxx

            const auto gd_geom_01_off = idx_geom_01_gdxx + i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxx_xx = cbuffer.data(gd_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxx_xy = cbuffer.data(gd_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxx_xz = cbuffer.data(gd_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxx_yy = cbuffer.data(gd_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxx_yz = cbuffer.data(gd_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxx_zz = cbuffer.data(gd_geom_01_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxx_xx, g_0_x_xxx_xxx, g_0_x_xxx_xxy, g_0_x_xxx_xxz, g_0_x_xxx_xy, g_0_x_xxx_xyy, g_0_x_xxx_xyz, g_0_x_xxx_xz, g_0_x_xxx_xzz, g_0_x_xxx_yy, g_0_x_xxx_yz, g_0_x_xxx_zz, g_0_x_xxxx_xx, g_0_x_xxxx_xy, g_0_x_xxxx_xz, g_0_x_xxxx_yy, g_0_x_xxxx_yz, g_0_x_xxxx_zz, g_xxx_xx, g_xxx_xy, g_xxx_xz, g_xxx_yy, g_xxx_yz, g_xxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxx_xx[k] = g_xxx_xx[k] - g_0_x_xxx_xx[k] * ab_x + g_0_x_xxx_xxx[k];

                g_0_x_xxxx_xy[k] = g_xxx_xy[k] - g_0_x_xxx_xy[k] * ab_x + g_0_x_xxx_xxy[k];

                g_0_x_xxxx_xz[k] = g_xxx_xz[k] - g_0_x_xxx_xz[k] * ab_x + g_0_x_xxx_xxz[k];

                g_0_x_xxxx_yy[k] = g_xxx_yy[k] - g_0_x_xxx_yy[k] * ab_x + g_0_x_xxx_xyy[k];

                g_0_x_xxxx_yz[k] = g_xxx_yz[k] - g_0_x_xxx_yz[k] * ab_x + g_0_x_xxx_xyz[k];

                g_0_x_xxxx_zz[k] = g_xxx_zz[k] - g_0_x_xxx_zz[k] * ab_x + g_0_x_xxx_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxy_xx = cbuffer.data(gd_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxy_xy = cbuffer.data(gd_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxy_xz = cbuffer.data(gd_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxy_yy = cbuffer.data(gd_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxy_yz = cbuffer.data(gd_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxy_zz = cbuffer.data(gd_geom_01_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxx_xx, g_0_x_xxx_xxy, g_0_x_xxx_xy, g_0_x_xxx_xyy, g_0_x_xxx_xyz, g_0_x_xxx_xz, g_0_x_xxx_yy, g_0_x_xxx_yyy, g_0_x_xxx_yyz, g_0_x_xxx_yz, g_0_x_xxx_yzz, g_0_x_xxx_zz, g_0_x_xxxy_xx, g_0_x_xxxy_xy, g_0_x_xxxy_xz, g_0_x_xxxy_yy, g_0_x_xxxy_yz, g_0_x_xxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxy_xx[k] = -g_0_x_xxx_xx[k] * ab_y + g_0_x_xxx_xxy[k];

                g_0_x_xxxy_xy[k] = -g_0_x_xxx_xy[k] * ab_y + g_0_x_xxx_xyy[k];

                g_0_x_xxxy_xz[k] = -g_0_x_xxx_xz[k] * ab_y + g_0_x_xxx_xyz[k];

                g_0_x_xxxy_yy[k] = -g_0_x_xxx_yy[k] * ab_y + g_0_x_xxx_yyy[k];

                g_0_x_xxxy_yz[k] = -g_0_x_xxx_yz[k] * ab_y + g_0_x_xxx_yyz[k];

                g_0_x_xxxy_zz[k] = -g_0_x_xxx_zz[k] * ab_y + g_0_x_xxx_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxz_xx = cbuffer.data(gd_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxz_xy = cbuffer.data(gd_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxz_xz = cbuffer.data(gd_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxz_yy = cbuffer.data(gd_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxz_yz = cbuffer.data(gd_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxz_zz = cbuffer.data(gd_geom_01_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxx_xx, g_0_x_xxx_xxz, g_0_x_xxx_xy, g_0_x_xxx_xyz, g_0_x_xxx_xz, g_0_x_xxx_xzz, g_0_x_xxx_yy, g_0_x_xxx_yyz, g_0_x_xxx_yz, g_0_x_xxx_yzz, g_0_x_xxx_zz, g_0_x_xxx_zzz, g_0_x_xxxz_xx, g_0_x_xxxz_xy, g_0_x_xxxz_xz, g_0_x_xxxz_yy, g_0_x_xxxz_yz, g_0_x_xxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxz_xx[k] = -g_0_x_xxx_xx[k] * ab_z + g_0_x_xxx_xxz[k];

                g_0_x_xxxz_xy[k] = -g_0_x_xxx_xy[k] * ab_z + g_0_x_xxx_xyz[k];

                g_0_x_xxxz_xz[k] = -g_0_x_xxx_xz[k] * ab_z + g_0_x_xxx_xzz[k];

                g_0_x_xxxz_yy[k] = -g_0_x_xxx_yy[k] * ab_z + g_0_x_xxx_yyz[k];

                g_0_x_xxxz_yz[k] = -g_0_x_xxx_yz[k] * ab_z + g_0_x_xxx_yzz[k];

                g_0_x_xxxz_zz[k] = -g_0_x_xxx_zz[k] * ab_z + g_0_x_xxx_zzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyy_xx = cbuffer.data(gd_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxyy_xy = cbuffer.data(gd_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxyy_xz = cbuffer.data(gd_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxyy_yy = cbuffer.data(gd_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxyy_yz = cbuffer.data(gd_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxyy_zz = cbuffer.data(gd_geom_01_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxy_xx, g_0_x_xxy_xxy, g_0_x_xxy_xy, g_0_x_xxy_xyy, g_0_x_xxy_xyz, g_0_x_xxy_xz, g_0_x_xxy_yy, g_0_x_xxy_yyy, g_0_x_xxy_yyz, g_0_x_xxy_yz, g_0_x_xxy_yzz, g_0_x_xxy_zz, g_0_x_xxyy_xx, g_0_x_xxyy_xy, g_0_x_xxyy_xz, g_0_x_xxyy_yy, g_0_x_xxyy_yz, g_0_x_xxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyy_xx[k] = -g_0_x_xxy_xx[k] * ab_y + g_0_x_xxy_xxy[k];

                g_0_x_xxyy_xy[k] = -g_0_x_xxy_xy[k] * ab_y + g_0_x_xxy_xyy[k];

                g_0_x_xxyy_xz[k] = -g_0_x_xxy_xz[k] * ab_y + g_0_x_xxy_xyz[k];

                g_0_x_xxyy_yy[k] = -g_0_x_xxy_yy[k] * ab_y + g_0_x_xxy_yyy[k];

                g_0_x_xxyy_yz[k] = -g_0_x_xxy_yz[k] * ab_y + g_0_x_xxy_yyz[k];

                g_0_x_xxyy_zz[k] = -g_0_x_xxy_zz[k] * ab_y + g_0_x_xxy_yzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyz_xx = cbuffer.data(gd_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxyz_xy = cbuffer.data(gd_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxyz_xz = cbuffer.data(gd_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxyz_yy = cbuffer.data(gd_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxyz_yz = cbuffer.data(gd_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxyz_zz = cbuffer.data(gd_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyz_xx, g_0_x_xxyz_xy, g_0_x_xxyz_xz, g_0_x_xxyz_yy, g_0_x_xxyz_yz, g_0_x_xxyz_zz, g_0_x_xxz_xx, g_0_x_xxz_xxy, g_0_x_xxz_xy, g_0_x_xxz_xyy, g_0_x_xxz_xyz, g_0_x_xxz_xz, g_0_x_xxz_yy, g_0_x_xxz_yyy, g_0_x_xxz_yyz, g_0_x_xxz_yz, g_0_x_xxz_yzz, g_0_x_xxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyz_xx[k] = -g_0_x_xxz_xx[k] * ab_y + g_0_x_xxz_xxy[k];

                g_0_x_xxyz_xy[k] = -g_0_x_xxz_xy[k] * ab_y + g_0_x_xxz_xyy[k];

                g_0_x_xxyz_xz[k] = -g_0_x_xxz_xz[k] * ab_y + g_0_x_xxz_xyz[k];

                g_0_x_xxyz_yy[k] = -g_0_x_xxz_yy[k] * ab_y + g_0_x_xxz_yyy[k];

                g_0_x_xxyz_yz[k] = -g_0_x_xxz_yz[k] * ab_y + g_0_x_xxz_yyz[k];

                g_0_x_xxyz_zz[k] = -g_0_x_xxz_zz[k] * ab_y + g_0_x_xxz_yzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxzz_xx = cbuffer.data(gd_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxzz_xy = cbuffer.data(gd_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxzz_xz = cbuffer.data(gd_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxzz_yy = cbuffer.data(gd_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxzz_yz = cbuffer.data(gd_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxzz_zz = cbuffer.data(gd_geom_01_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxz_xx, g_0_x_xxz_xxz, g_0_x_xxz_xy, g_0_x_xxz_xyz, g_0_x_xxz_xz, g_0_x_xxz_xzz, g_0_x_xxz_yy, g_0_x_xxz_yyz, g_0_x_xxz_yz, g_0_x_xxz_yzz, g_0_x_xxz_zz, g_0_x_xxz_zzz, g_0_x_xxzz_xx, g_0_x_xxzz_xy, g_0_x_xxzz_xz, g_0_x_xxzz_yy, g_0_x_xxzz_yz, g_0_x_xxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxzz_xx[k] = -g_0_x_xxz_xx[k] * ab_z + g_0_x_xxz_xxz[k];

                g_0_x_xxzz_xy[k] = -g_0_x_xxz_xy[k] * ab_z + g_0_x_xxz_xyz[k];

                g_0_x_xxzz_xz[k] = -g_0_x_xxz_xz[k] * ab_z + g_0_x_xxz_xzz[k];

                g_0_x_xxzz_yy[k] = -g_0_x_xxz_yy[k] * ab_z + g_0_x_xxz_yyz[k];

                g_0_x_xxzz_yz[k] = -g_0_x_xxz_yz[k] * ab_z + g_0_x_xxz_yzz[k];

                g_0_x_xxzz_zz[k] = -g_0_x_xxz_zz[k] * ab_z + g_0_x_xxz_zzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyy_xx = cbuffer.data(gd_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xyyy_xy = cbuffer.data(gd_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xyyy_xz = cbuffer.data(gd_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xyyy_yy = cbuffer.data(gd_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xyyy_yz = cbuffer.data(gd_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xyyy_zz = cbuffer.data(gd_geom_01_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyy_xx, g_0_x_xyy_xxy, g_0_x_xyy_xy, g_0_x_xyy_xyy, g_0_x_xyy_xyz, g_0_x_xyy_xz, g_0_x_xyy_yy, g_0_x_xyy_yyy, g_0_x_xyy_yyz, g_0_x_xyy_yz, g_0_x_xyy_yzz, g_0_x_xyy_zz, g_0_x_xyyy_xx, g_0_x_xyyy_xy, g_0_x_xyyy_xz, g_0_x_xyyy_yy, g_0_x_xyyy_yz, g_0_x_xyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyy_xx[k] = -g_0_x_xyy_xx[k] * ab_y + g_0_x_xyy_xxy[k];

                g_0_x_xyyy_xy[k] = -g_0_x_xyy_xy[k] * ab_y + g_0_x_xyy_xyy[k];

                g_0_x_xyyy_xz[k] = -g_0_x_xyy_xz[k] * ab_y + g_0_x_xyy_xyz[k];

                g_0_x_xyyy_yy[k] = -g_0_x_xyy_yy[k] * ab_y + g_0_x_xyy_yyy[k];

                g_0_x_xyyy_yz[k] = -g_0_x_xyy_yz[k] * ab_y + g_0_x_xyy_yyz[k];

                g_0_x_xyyy_zz[k] = -g_0_x_xyy_zz[k] * ab_y + g_0_x_xyy_yzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyz_xx = cbuffer.data(gd_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xyyz_xy = cbuffer.data(gd_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xyyz_xz = cbuffer.data(gd_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xyyz_yy = cbuffer.data(gd_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xyyz_yz = cbuffer.data(gd_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xyyz_zz = cbuffer.data(gd_geom_01_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyz_xx, g_0_x_xyyz_xy, g_0_x_xyyz_xz, g_0_x_xyyz_yy, g_0_x_xyyz_yz, g_0_x_xyyz_zz, g_0_x_xyz_xx, g_0_x_xyz_xxy, g_0_x_xyz_xy, g_0_x_xyz_xyy, g_0_x_xyz_xyz, g_0_x_xyz_xz, g_0_x_xyz_yy, g_0_x_xyz_yyy, g_0_x_xyz_yyz, g_0_x_xyz_yz, g_0_x_xyz_yzz, g_0_x_xyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyz_xx[k] = -g_0_x_xyz_xx[k] * ab_y + g_0_x_xyz_xxy[k];

                g_0_x_xyyz_xy[k] = -g_0_x_xyz_xy[k] * ab_y + g_0_x_xyz_xyy[k];

                g_0_x_xyyz_xz[k] = -g_0_x_xyz_xz[k] * ab_y + g_0_x_xyz_xyz[k];

                g_0_x_xyyz_yy[k] = -g_0_x_xyz_yy[k] * ab_y + g_0_x_xyz_yyy[k];

                g_0_x_xyyz_yz[k] = -g_0_x_xyz_yz[k] * ab_y + g_0_x_xyz_yyz[k];

                g_0_x_xyyz_zz[k] = -g_0_x_xyz_zz[k] * ab_y + g_0_x_xyz_yzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyzz_xx = cbuffer.data(gd_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xyzz_xy = cbuffer.data(gd_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xyzz_xz = cbuffer.data(gd_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xyzz_yy = cbuffer.data(gd_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xyzz_yz = cbuffer.data(gd_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xyzz_zz = cbuffer.data(gd_geom_01_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyzz_xx, g_0_x_xyzz_xy, g_0_x_xyzz_xz, g_0_x_xyzz_yy, g_0_x_xyzz_yz, g_0_x_xyzz_zz, g_0_x_xzz_xx, g_0_x_xzz_xxy, g_0_x_xzz_xy, g_0_x_xzz_xyy, g_0_x_xzz_xyz, g_0_x_xzz_xz, g_0_x_xzz_yy, g_0_x_xzz_yyy, g_0_x_xzz_yyz, g_0_x_xzz_yz, g_0_x_xzz_yzz, g_0_x_xzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyzz_xx[k] = -g_0_x_xzz_xx[k] * ab_y + g_0_x_xzz_xxy[k];

                g_0_x_xyzz_xy[k] = -g_0_x_xzz_xy[k] * ab_y + g_0_x_xzz_xyy[k];

                g_0_x_xyzz_xz[k] = -g_0_x_xzz_xz[k] * ab_y + g_0_x_xzz_xyz[k];

                g_0_x_xyzz_yy[k] = -g_0_x_xzz_yy[k] * ab_y + g_0_x_xzz_yyy[k];

                g_0_x_xyzz_yz[k] = -g_0_x_xzz_yz[k] * ab_y + g_0_x_xzz_yyz[k];

                g_0_x_xyzz_zz[k] = -g_0_x_xzz_zz[k] * ab_y + g_0_x_xzz_yzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_0_x_xzzz_xx = cbuffer.data(gd_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xzzz_xy = cbuffer.data(gd_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xzzz_xz = cbuffer.data(gd_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xzzz_yy = cbuffer.data(gd_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xzzz_yz = cbuffer.data(gd_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xzzz_zz = cbuffer.data(gd_geom_01_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xzz_xx, g_0_x_xzz_xxz, g_0_x_xzz_xy, g_0_x_xzz_xyz, g_0_x_xzz_xz, g_0_x_xzz_xzz, g_0_x_xzz_yy, g_0_x_xzz_yyz, g_0_x_xzz_yz, g_0_x_xzz_yzz, g_0_x_xzz_zz, g_0_x_xzz_zzz, g_0_x_xzzz_xx, g_0_x_xzzz_xy, g_0_x_xzzz_xz, g_0_x_xzzz_yy, g_0_x_xzzz_yz, g_0_x_xzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzzz_xx[k] = -g_0_x_xzz_xx[k] * ab_z + g_0_x_xzz_xxz[k];

                g_0_x_xzzz_xy[k] = -g_0_x_xzz_xy[k] * ab_z + g_0_x_xzz_xyz[k];

                g_0_x_xzzz_xz[k] = -g_0_x_xzz_xz[k] * ab_z + g_0_x_xzz_xzz[k];

                g_0_x_xzzz_yy[k] = -g_0_x_xzz_yy[k] * ab_z + g_0_x_xzz_yyz[k];

                g_0_x_xzzz_yz[k] = -g_0_x_xzz_yz[k] * ab_z + g_0_x_xzz_yzz[k];

                g_0_x_xzzz_zz[k] = -g_0_x_xzz_zz[k] * ab_z + g_0_x_xzz_zzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyy_xx = cbuffer.data(gd_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_yyyy_xy = cbuffer.data(gd_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_yyyy_xz = cbuffer.data(gd_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_yyyy_yy = cbuffer.data(gd_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_yyyy_yz = cbuffer.data(gd_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_yyyy_zz = cbuffer.data(gd_geom_01_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyy_xx, g_0_x_yyy_xxy, g_0_x_yyy_xy, g_0_x_yyy_xyy, g_0_x_yyy_xyz, g_0_x_yyy_xz, g_0_x_yyy_yy, g_0_x_yyy_yyy, g_0_x_yyy_yyz, g_0_x_yyy_yz, g_0_x_yyy_yzz, g_0_x_yyy_zz, g_0_x_yyyy_xx, g_0_x_yyyy_xy, g_0_x_yyyy_xz, g_0_x_yyyy_yy, g_0_x_yyyy_yz, g_0_x_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyy_xx[k] = -g_0_x_yyy_xx[k] * ab_y + g_0_x_yyy_xxy[k];

                g_0_x_yyyy_xy[k] = -g_0_x_yyy_xy[k] * ab_y + g_0_x_yyy_xyy[k];

                g_0_x_yyyy_xz[k] = -g_0_x_yyy_xz[k] * ab_y + g_0_x_yyy_xyz[k];

                g_0_x_yyyy_yy[k] = -g_0_x_yyy_yy[k] * ab_y + g_0_x_yyy_yyy[k];

                g_0_x_yyyy_yz[k] = -g_0_x_yyy_yz[k] * ab_y + g_0_x_yyy_yyz[k];

                g_0_x_yyyy_zz[k] = -g_0_x_yyy_zz[k] * ab_y + g_0_x_yyy_yzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyz_xx = cbuffer.data(gd_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_yyyz_xy = cbuffer.data(gd_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_yyyz_xz = cbuffer.data(gd_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_yyyz_yy = cbuffer.data(gd_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_yyyz_yz = cbuffer.data(gd_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_yyyz_zz = cbuffer.data(gd_geom_01_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyz_xx, g_0_x_yyyz_xy, g_0_x_yyyz_xz, g_0_x_yyyz_yy, g_0_x_yyyz_yz, g_0_x_yyyz_zz, g_0_x_yyz_xx, g_0_x_yyz_xxy, g_0_x_yyz_xy, g_0_x_yyz_xyy, g_0_x_yyz_xyz, g_0_x_yyz_xz, g_0_x_yyz_yy, g_0_x_yyz_yyy, g_0_x_yyz_yyz, g_0_x_yyz_yz, g_0_x_yyz_yzz, g_0_x_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyz_xx[k] = -g_0_x_yyz_xx[k] * ab_y + g_0_x_yyz_xxy[k];

                g_0_x_yyyz_xy[k] = -g_0_x_yyz_xy[k] * ab_y + g_0_x_yyz_xyy[k];

                g_0_x_yyyz_xz[k] = -g_0_x_yyz_xz[k] * ab_y + g_0_x_yyz_xyz[k];

                g_0_x_yyyz_yy[k] = -g_0_x_yyz_yy[k] * ab_y + g_0_x_yyz_yyy[k];

                g_0_x_yyyz_yz[k] = -g_0_x_yyz_yz[k] * ab_y + g_0_x_yyz_yyz[k];

                g_0_x_yyyz_zz[k] = -g_0_x_yyz_zz[k] * ab_y + g_0_x_yyz_yzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyzz_xx = cbuffer.data(gd_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_yyzz_xy = cbuffer.data(gd_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_yyzz_xz = cbuffer.data(gd_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_yyzz_yy = cbuffer.data(gd_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_yyzz_yz = cbuffer.data(gd_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_yyzz_zz = cbuffer.data(gd_geom_01_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyzz_xx, g_0_x_yyzz_xy, g_0_x_yyzz_xz, g_0_x_yyzz_yy, g_0_x_yyzz_yz, g_0_x_yyzz_zz, g_0_x_yzz_xx, g_0_x_yzz_xxy, g_0_x_yzz_xy, g_0_x_yzz_xyy, g_0_x_yzz_xyz, g_0_x_yzz_xz, g_0_x_yzz_yy, g_0_x_yzz_yyy, g_0_x_yzz_yyz, g_0_x_yzz_yz, g_0_x_yzz_yzz, g_0_x_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyzz_xx[k] = -g_0_x_yzz_xx[k] * ab_y + g_0_x_yzz_xxy[k];

                g_0_x_yyzz_xy[k] = -g_0_x_yzz_xy[k] * ab_y + g_0_x_yzz_xyy[k];

                g_0_x_yyzz_xz[k] = -g_0_x_yzz_xz[k] * ab_y + g_0_x_yzz_xyz[k];

                g_0_x_yyzz_yy[k] = -g_0_x_yzz_yy[k] * ab_y + g_0_x_yzz_yyy[k];

                g_0_x_yyzz_yz[k] = -g_0_x_yzz_yz[k] * ab_y + g_0_x_yzz_yyz[k];

                g_0_x_yyzz_zz[k] = -g_0_x_yzz_zz[k] * ab_y + g_0_x_yzz_yzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_0_x_yzzz_xx = cbuffer.data(gd_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_yzzz_xy = cbuffer.data(gd_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_yzzz_xz = cbuffer.data(gd_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_yzzz_yy = cbuffer.data(gd_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_yzzz_yz = cbuffer.data(gd_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_yzzz_zz = cbuffer.data(gd_geom_01_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yzzz_xx, g_0_x_yzzz_xy, g_0_x_yzzz_xz, g_0_x_yzzz_yy, g_0_x_yzzz_yz, g_0_x_yzzz_zz, g_0_x_zzz_xx, g_0_x_zzz_xxy, g_0_x_zzz_xy, g_0_x_zzz_xyy, g_0_x_zzz_xyz, g_0_x_zzz_xz, g_0_x_zzz_yy, g_0_x_zzz_yyy, g_0_x_zzz_yyz, g_0_x_zzz_yz, g_0_x_zzz_yzz, g_0_x_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzzz_xx[k] = -g_0_x_zzz_xx[k] * ab_y + g_0_x_zzz_xxy[k];

                g_0_x_yzzz_xy[k] = -g_0_x_zzz_xy[k] * ab_y + g_0_x_zzz_xyy[k];

                g_0_x_yzzz_xz[k] = -g_0_x_zzz_xz[k] * ab_y + g_0_x_zzz_xyz[k];

                g_0_x_yzzz_yy[k] = -g_0_x_zzz_yy[k] * ab_y + g_0_x_zzz_yyy[k];

                g_0_x_yzzz_yz[k] = -g_0_x_zzz_yz[k] * ab_y + g_0_x_zzz_yyz[k];

                g_0_x_yzzz_zz[k] = -g_0_x_zzz_zz[k] * ab_y + g_0_x_zzz_yzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_0_x_zzzz_xx = cbuffer.data(gd_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_zzzz_xy = cbuffer.data(gd_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_zzzz_xz = cbuffer.data(gd_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_zzzz_yy = cbuffer.data(gd_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_zzzz_yz = cbuffer.data(gd_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_zzzz_zz = cbuffer.data(gd_geom_01_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzz_xx, g_0_x_zzz_xxz, g_0_x_zzz_xy, g_0_x_zzz_xyz, g_0_x_zzz_xz, g_0_x_zzz_xzz, g_0_x_zzz_yy, g_0_x_zzz_yyz, g_0_x_zzz_yz, g_0_x_zzz_yzz, g_0_x_zzz_zz, g_0_x_zzz_zzz, g_0_x_zzzz_xx, g_0_x_zzzz_xy, g_0_x_zzzz_xz, g_0_x_zzzz_yy, g_0_x_zzzz_yz, g_0_x_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzzz_xx[k] = -g_0_x_zzz_xx[k] * ab_z + g_0_x_zzz_xxz[k];

                g_0_x_zzzz_xy[k] = -g_0_x_zzz_xy[k] * ab_z + g_0_x_zzz_xyz[k];

                g_0_x_zzzz_xz[k] = -g_0_x_zzz_xz[k] * ab_z + g_0_x_zzz_xzz[k];

                g_0_x_zzzz_yy[k] = -g_0_x_zzz_yy[k] * ab_z + g_0_x_zzz_yyz[k];

                g_0_x_zzzz_yz[k] = -g_0_x_zzz_yz[k] * ab_z + g_0_x_zzz_yzz[k];

                g_0_x_zzzz_zz[k] = -g_0_x_zzz_zz[k] * ab_z + g_0_x_zzz_zzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxx_xx = cbuffer.data(gd_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_y_xxxx_xy = cbuffer.data(gd_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_y_xxxx_xz = cbuffer.data(gd_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_y_xxxx_yy = cbuffer.data(gd_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_y_xxxx_yz = cbuffer.data(gd_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_y_xxxx_zz = cbuffer.data(gd_geom_01_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxx_xx, g_0_y_xxx_xxx, g_0_y_xxx_xxy, g_0_y_xxx_xxz, g_0_y_xxx_xy, g_0_y_xxx_xyy, g_0_y_xxx_xyz, g_0_y_xxx_xz, g_0_y_xxx_xzz, g_0_y_xxx_yy, g_0_y_xxx_yz, g_0_y_xxx_zz, g_0_y_xxxx_xx, g_0_y_xxxx_xy, g_0_y_xxxx_xz, g_0_y_xxxx_yy, g_0_y_xxxx_yz, g_0_y_xxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxx_xx[k] = -g_0_y_xxx_xx[k] * ab_x + g_0_y_xxx_xxx[k];

                g_0_y_xxxx_xy[k] = -g_0_y_xxx_xy[k] * ab_x + g_0_y_xxx_xxy[k];

                g_0_y_xxxx_xz[k] = -g_0_y_xxx_xz[k] * ab_x + g_0_y_xxx_xxz[k];

                g_0_y_xxxx_yy[k] = -g_0_y_xxx_yy[k] * ab_x + g_0_y_xxx_xyy[k];

                g_0_y_xxxx_yz[k] = -g_0_y_xxx_yz[k] * ab_x + g_0_y_xxx_xyz[k];

                g_0_y_xxxx_zz[k] = -g_0_y_xxx_zz[k] * ab_x + g_0_y_xxx_xzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxy_xx = cbuffer.data(gd_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_y_xxxy_xy = cbuffer.data(gd_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_y_xxxy_xz = cbuffer.data(gd_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_y_xxxy_yy = cbuffer.data(gd_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_y_xxxy_yz = cbuffer.data(gd_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_y_xxxy_zz = cbuffer.data(gd_geom_01_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxy_xx, g_0_y_xxxy_xy, g_0_y_xxxy_xz, g_0_y_xxxy_yy, g_0_y_xxxy_yz, g_0_y_xxxy_zz, g_0_y_xxy_xx, g_0_y_xxy_xxx, g_0_y_xxy_xxy, g_0_y_xxy_xxz, g_0_y_xxy_xy, g_0_y_xxy_xyy, g_0_y_xxy_xyz, g_0_y_xxy_xz, g_0_y_xxy_xzz, g_0_y_xxy_yy, g_0_y_xxy_yz, g_0_y_xxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxy_xx[k] = -g_0_y_xxy_xx[k] * ab_x + g_0_y_xxy_xxx[k];

                g_0_y_xxxy_xy[k] = -g_0_y_xxy_xy[k] * ab_x + g_0_y_xxy_xxy[k];

                g_0_y_xxxy_xz[k] = -g_0_y_xxy_xz[k] * ab_x + g_0_y_xxy_xxz[k];

                g_0_y_xxxy_yy[k] = -g_0_y_xxy_yy[k] * ab_x + g_0_y_xxy_xyy[k];

                g_0_y_xxxy_yz[k] = -g_0_y_xxy_yz[k] * ab_x + g_0_y_xxy_xyz[k];

                g_0_y_xxxy_zz[k] = -g_0_y_xxy_zz[k] * ab_x + g_0_y_xxy_xzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxz_xx = cbuffer.data(gd_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_y_xxxz_xy = cbuffer.data(gd_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_y_xxxz_xz = cbuffer.data(gd_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_y_xxxz_yy = cbuffer.data(gd_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_y_xxxz_yz = cbuffer.data(gd_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_y_xxxz_zz = cbuffer.data(gd_geom_01_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxz_xx, g_0_y_xxxz_xy, g_0_y_xxxz_xz, g_0_y_xxxz_yy, g_0_y_xxxz_yz, g_0_y_xxxz_zz, g_0_y_xxz_xx, g_0_y_xxz_xxx, g_0_y_xxz_xxy, g_0_y_xxz_xxz, g_0_y_xxz_xy, g_0_y_xxz_xyy, g_0_y_xxz_xyz, g_0_y_xxz_xz, g_0_y_xxz_xzz, g_0_y_xxz_yy, g_0_y_xxz_yz, g_0_y_xxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxz_xx[k] = -g_0_y_xxz_xx[k] * ab_x + g_0_y_xxz_xxx[k];

                g_0_y_xxxz_xy[k] = -g_0_y_xxz_xy[k] * ab_x + g_0_y_xxz_xxy[k];

                g_0_y_xxxz_xz[k] = -g_0_y_xxz_xz[k] * ab_x + g_0_y_xxz_xxz[k];

                g_0_y_xxxz_yy[k] = -g_0_y_xxz_yy[k] * ab_x + g_0_y_xxz_xyy[k];

                g_0_y_xxxz_yz[k] = -g_0_y_xxz_yz[k] * ab_x + g_0_y_xxz_xyz[k];

                g_0_y_xxxz_zz[k] = -g_0_y_xxz_zz[k] * ab_x + g_0_y_xxz_xzz[k];
            }

            /// Set up 108-114 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyy_xx = cbuffer.data(gd_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_xxyy_xy = cbuffer.data(gd_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_xxyy_xz = cbuffer.data(gd_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_y_xxyy_yy = cbuffer.data(gd_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_y_xxyy_yz = cbuffer.data(gd_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_xxyy_zz = cbuffer.data(gd_geom_01_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyy_xx, g_0_y_xxyy_xy, g_0_y_xxyy_xz, g_0_y_xxyy_yy, g_0_y_xxyy_yz, g_0_y_xxyy_zz, g_0_y_xyy_xx, g_0_y_xyy_xxx, g_0_y_xyy_xxy, g_0_y_xyy_xxz, g_0_y_xyy_xy, g_0_y_xyy_xyy, g_0_y_xyy_xyz, g_0_y_xyy_xz, g_0_y_xyy_xzz, g_0_y_xyy_yy, g_0_y_xyy_yz, g_0_y_xyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyy_xx[k] = -g_0_y_xyy_xx[k] * ab_x + g_0_y_xyy_xxx[k];

                g_0_y_xxyy_xy[k] = -g_0_y_xyy_xy[k] * ab_x + g_0_y_xyy_xxy[k];

                g_0_y_xxyy_xz[k] = -g_0_y_xyy_xz[k] * ab_x + g_0_y_xyy_xxz[k];

                g_0_y_xxyy_yy[k] = -g_0_y_xyy_yy[k] * ab_x + g_0_y_xyy_xyy[k];

                g_0_y_xxyy_yz[k] = -g_0_y_xyy_yz[k] * ab_x + g_0_y_xyy_xyz[k];

                g_0_y_xxyy_zz[k] = -g_0_y_xyy_zz[k] * ab_x + g_0_y_xyy_xzz[k];
            }

            /// Set up 114-120 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyz_xx = cbuffer.data(gd_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_xxyz_xy = cbuffer.data(gd_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_xxyz_xz = cbuffer.data(gd_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_y_xxyz_yy = cbuffer.data(gd_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_xxyz_yz = cbuffer.data(gd_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_xxyz_zz = cbuffer.data(gd_geom_01_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyz_xx, g_0_y_xxyz_xy, g_0_y_xxyz_xz, g_0_y_xxyz_yy, g_0_y_xxyz_yz, g_0_y_xxyz_zz, g_0_y_xyz_xx, g_0_y_xyz_xxx, g_0_y_xyz_xxy, g_0_y_xyz_xxz, g_0_y_xyz_xy, g_0_y_xyz_xyy, g_0_y_xyz_xyz, g_0_y_xyz_xz, g_0_y_xyz_xzz, g_0_y_xyz_yy, g_0_y_xyz_yz, g_0_y_xyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyz_xx[k] = -g_0_y_xyz_xx[k] * ab_x + g_0_y_xyz_xxx[k];

                g_0_y_xxyz_xy[k] = -g_0_y_xyz_xy[k] * ab_x + g_0_y_xyz_xxy[k];

                g_0_y_xxyz_xz[k] = -g_0_y_xyz_xz[k] * ab_x + g_0_y_xyz_xxz[k];

                g_0_y_xxyz_yy[k] = -g_0_y_xyz_yy[k] * ab_x + g_0_y_xyz_xyy[k];

                g_0_y_xxyz_yz[k] = -g_0_y_xyz_yz[k] * ab_x + g_0_y_xyz_xyz[k];

                g_0_y_xxyz_zz[k] = -g_0_y_xyz_zz[k] * ab_x + g_0_y_xyz_xzz[k];
            }

            /// Set up 120-126 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxzz_xx = cbuffer.data(gd_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_y_xxzz_xy = cbuffer.data(gd_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_y_xxzz_xz = cbuffer.data(gd_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_y_xxzz_yy = cbuffer.data(gd_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_y_xxzz_yz = cbuffer.data(gd_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_y_xxzz_zz = cbuffer.data(gd_geom_01_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxzz_xx, g_0_y_xxzz_xy, g_0_y_xxzz_xz, g_0_y_xxzz_yy, g_0_y_xxzz_yz, g_0_y_xxzz_zz, g_0_y_xzz_xx, g_0_y_xzz_xxx, g_0_y_xzz_xxy, g_0_y_xzz_xxz, g_0_y_xzz_xy, g_0_y_xzz_xyy, g_0_y_xzz_xyz, g_0_y_xzz_xz, g_0_y_xzz_xzz, g_0_y_xzz_yy, g_0_y_xzz_yz, g_0_y_xzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxzz_xx[k] = -g_0_y_xzz_xx[k] * ab_x + g_0_y_xzz_xxx[k];

                g_0_y_xxzz_xy[k] = -g_0_y_xzz_xy[k] * ab_x + g_0_y_xzz_xxy[k];

                g_0_y_xxzz_xz[k] = -g_0_y_xzz_xz[k] * ab_x + g_0_y_xzz_xxz[k];

                g_0_y_xxzz_yy[k] = -g_0_y_xzz_yy[k] * ab_x + g_0_y_xzz_xyy[k];

                g_0_y_xxzz_yz[k] = -g_0_y_xzz_yz[k] * ab_x + g_0_y_xzz_xyz[k];

                g_0_y_xxzz_zz[k] = -g_0_y_xzz_zz[k] * ab_x + g_0_y_xzz_xzz[k];
            }

            /// Set up 126-132 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyy_xx = cbuffer.data(gd_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_y_xyyy_xy = cbuffer.data(gd_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_y_xyyy_xz = cbuffer.data(gd_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_y_xyyy_yy = cbuffer.data(gd_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_y_xyyy_yz = cbuffer.data(gd_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_y_xyyy_zz = cbuffer.data(gd_geom_01_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyy_xx, g_0_y_xyyy_xy, g_0_y_xyyy_xz, g_0_y_xyyy_yy, g_0_y_xyyy_yz, g_0_y_xyyy_zz, g_0_y_yyy_xx, g_0_y_yyy_xxx, g_0_y_yyy_xxy, g_0_y_yyy_xxz, g_0_y_yyy_xy, g_0_y_yyy_xyy, g_0_y_yyy_xyz, g_0_y_yyy_xz, g_0_y_yyy_xzz, g_0_y_yyy_yy, g_0_y_yyy_yz, g_0_y_yyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyy_xx[k] = -g_0_y_yyy_xx[k] * ab_x + g_0_y_yyy_xxx[k];

                g_0_y_xyyy_xy[k] = -g_0_y_yyy_xy[k] * ab_x + g_0_y_yyy_xxy[k];

                g_0_y_xyyy_xz[k] = -g_0_y_yyy_xz[k] * ab_x + g_0_y_yyy_xxz[k];

                g_0_y_xyyy_yy[k] = -g_0_y_yyy_yy[k] * ab_x + g_0_y_yyy_xyy[k];

                g_0_y_xyyy_yz[k] = -g_0_y_yyy_yz[k] * ab_x + g_0_y_yyy_xyz[k];

                g_0_y_xyyy_zz[k] = -g_0_y_yyy_zz[k] * ab_x + g_0_y_yyy_xzz[k];
            }

            /// Set up 132-138 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyz_xx = cbuffer.data(gd_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_y_xyyz_xy = cbuffer.data(gd_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_y_xyyz_xz = cbuffer.data(gd_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_y_xyyz_yy = cbuffer.data(gd_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_xyyz_yz = cbuffer.data(gd_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_xyyz_zz = cbuffer.data(gd_geom_01_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyz_xx, g_0_y_xyyz_xy, g_0_y_xyyz_xz, g_0_y_xyyz_yy, g_0_y_xyyz_yz, g_0_y_xyyz_zz, g_0_y_yyz_xx, g_0_y_yyz_xxx, g_0_y_yyz_xxy, g_0_y_yyz_xxz, g_0_y_yyz_xy, g_0_y_yyz_xyy, g_0_y_yyz_xyz, g_0_y_yyz_xz, g_0_y_yyz_xzz, g_0_y_yyz_yy, g_0_y_yyz_yz, g_0_y_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyz_xx[k] = -g_0_y_yyz_xx[k] * ab_x + g_0_y_yyz_xxx[k];

                g_0_y_xyyz_xy[k] = -g_0_y_yyz_xy[k] * ab_x + g_0_y_yyz_xxy[k];

                g_0_y_xyyz_xz[k] = -g_0_y_yyz_xz[k] * ab_x + g_0_y_yyz_xxz[k];

                g_0_y_xyyz_yy[k] = -g_0_y_yyz_yy[k] * ab_x + g_0_y_yyz_xyy[k];

                g_0_y_xyyz_yz[k] = -g_0_y_yyz_yz[k] * ab_x + g_0_y_yyz_xyz[k];

                g_0_y_xyyz_zz[k] = -g_0_y_yyz_zz[k] * ab_x + g_0_y_yyz_xzz[k];
            }

            /// Set up 138-144 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyzz_xx = cbuffer.data(gd_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_xyzz_xy = cbuffer.data(gd_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_y_xyzz_xz = cbuffer.data(gd_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_y_xyzz_yy = cbuffer.data(gd_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_xyzz_yz = cbuffer.data(gd_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_xyzz_zz = cbuffer.data(gd_geom_01_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyzz_xx, g_0_y_xyzz_xy, g_0_y_xyzz_xz, g_0_y_xyzz_yy, g_0_y_xyzz_yz, g_0_y_xyzz_zz, g_0_y_yzz_xx, g_0_y_yzz_xxx, g_0_y_yzz_xxy, g_0_y_yzz_xxz, g_0_y_yzz_xy, g_0_y_yzz_xyy, g_0_y_yzz_xyz, g_0_y_yzz_xz, g_0_y_yzz_xzz, g_0_y_yzz_yy, g_0_y_yzz_yz, g_0_y_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyzz_xx[k] = -g_0_y_yzz_xx[k] * ab_x + g_0_y_yzz_xxx[k];

                g_0_y_xyzz_xy[k] = -g_0_y_yzz_xy[k] * ab_x + g_0_y_yzz_xxy[k];

                g_0_y_xyzz_xz[k] = -g_0_y_yzz_xz[k] * ab_x + g_0_y_yzz_xxz[k];

                g_0_y_xyzz_yy[k] = -g_0_y_yzz_yy[k] * ab_x + g_0_y_yzz_xyy[k];

                g_0_y_xyzz_yz[k] = -g_0_y_yzz_yz[k] * ab_x + g_0_y_yzz_xyz[k];

                g_0_y_xyzz_zz[k] = -g_0_y_yzz_zz[k] * ab_x + g_0_y_yzz_xzz[k];
            }

            /// Set up 144-150 components of targeted buffer : cbuffer.data(

            auto g_0_y_xzzz_xx = cbuffer.data(gd_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_xzzz_xy = cbuffer.data(gd_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_xzzz_xz = cbuffer.data(gd_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_y_xzzz_yy = cbuffer.data(gd_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_xzzz_yz = cbuffer.data(gd_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_xzzz_zz = cbuffer.data(gd_geom_01_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xzzz_xx, g_0_y_xzzz_xy, g_0_y_xzzz_xz, g_0_y_xzzz_yy, g_0_y_xzzz_yz, g_0_y_xzzz_zz, g_0_y_zzz_xx, g_0_y_zzz_xxx, g_0_y_zzz_xxy, g_0_y_zzz_xxz, g_0_y_zzz_xy, g_0_y_zzz_xyy, g_0_y_zzz_xyz, g_0_y_zzz_xz, g_0_y_zzz_xzz, g_0_y_zzz_yy, g_0_y_zzz_yz, g_0_y_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzzz_xx[k] = -g_0_y_zzz_xx[k] * ab_x + g_0_y_zzz_xxx[k];

                g_0_y_xzzz_xy[k] = -g_0_y_zzz_xy[k] * ab_x + g_0_y_zzz_xxy[k];

                g_0_y_xzzz_xz[k] = -g_0_y_zzz_xz[k] * ab_x + g_0_y_zzz_xxz[k];

                g_0_y_xzzz_yy[k] = -g_0_y_zzz_yy[k] * ab_x + g_0_y_zzz_xyy[k];

                g_0_y_xzzz_yz[k] = -g_0_y_zzz_yz[k] * ab_x + g_0_y_zzz_xyz[k];

                g_0_y_xzzz_zz[k] = -g_0_y_zzz_zz[k] * ab_x + g_0_y_zzz_xzz[k];
            }

            /// Set up 150-156 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyy_xx = cbuffer.data(gd_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_yyyy_xy = cbuffer.data(gd_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_yyyy_xz = cbuffer.data(gd_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_yyyy_yy = cbuffer.data(gd_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_yyyy_yz = cbuffer.data(gd_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_yyyy_zz = cbuffer.data(gd_geom_01_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyy_xx, g_0_y_yyy_xxy, g_0_y_yyy_xy, g_0_y_yyy_xyy, g_0_y_yyy_xyz, g_0_y_yyy_xz, g_0_y_yyy_yy, g_0_y_yyy_yyy, g_0_y_yyy_yyz, g_0_y_yyy_yz, g_0_y_yyy_yzz, g_0_y_yyy_zz, g_0_y_yyyy_xx, g_0_y_yyyy_xy, g_0_y_yyyy_xz, g_0_y_yyyy_yy, g_0_y_yyyy_yz, g_0_y_yyyy_zz, g_yyy_xx, g_yyy_xy, g_yyy_xz, g_yyy_yy, g_yyy_yz, g_yyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyy_xx[k] = g_yyy_xx[k] - g_0_y_yyy_xx[k] * ab_y + g_0_y_yyy_xxy[k];

                g_0_y_yyyy_xy[k] = g_yyy_xy[k] - g_0_y_yyy_xy[k] * ab_y + g_0_y_yyy_xyy[k];

                g_0_y_yyyy_xz[k] = g_yyy_xz[k] - g_0_y_yyy_xz[k] * ab_y + g_0_y_yyy_xyz[k];

                g_0_y_yyyy_yy[k] = g_yyy_yy[k] - g_0_y_yyy_yy[k] * ab_y + g_0_y_yyy_yyy[k];

                g_0_y_yyyy_yz[k] = g_yyy_yz[k] - g_0_y_yyy_yz[k] * ab_y + g_0_y_yyy_yyz[k];

                g_0_y_yyyy_zz[k] = g_yyy_zz[k] - g_0_y_yyy_zz[k] * ab_y + g_0_y_yyy_yzz[k];
            }

            /// Set up 156-162 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyz_xx = cbuffer.data(gd_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_yyyz_xy = cbuffer.data(gd_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_yyyz_xz = cbuffer.data(gd_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_yyyz_yy = cbuffer.data(gd_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_yyyz_yz = cbuffer.data(gd_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_yyyz_zz = cbuffer.data(gd_geom_01_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyy_xx, g_0_y_yyy_xxz, g_0_y_yyy_xy, g_0_y_yyy_xyz, g_0_y_yyy_xz, g_0_y_yyy_xzz, g_0_y_yyy_yy, g_0_y_yyy_yyz, g_0_y_yyy_yz, g_0_y_yyy_yzz, g_0_y_yyy_zz, g_0_y_yyy_zzz, g_0_y_yyyz_xx, g_0_y_yyyz_xy, g_0_y_yyyz_xz, g_0_y_yyyz_yy, g_0_y_yyyz_yz, g_0_y_yyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyz_xx[k] = -g_0_y_yyy_xx[k] * ab_z + g_0_y_yyy_xxz[k];

                g_0_y_yyyz_xy[k] = -g_0_y_yyy_xy[k] * ab_z + g_0_y_yyy_xyz[k];

                g_0_y_yyyz_xz[k] = -g_0_y_yyy_xz[k] * ab_z + g_0_y_yyy_xzz[k];

                g_0_y_yyyz_yy[k] = -g_0_y_yyy_yy[k] * ab_z + g_0_y_yyy_yyz[k];

                g_0_y_yyyz_yz[k] = -g_0_y_yyy_yz[k] * ab_z + g_0_y_yyy_yzz[k];

                g_0_y_yyyz_zz[k] = -g_0_y_yyy_zz[k] * ab_z + g_0_y_yyy_zzz[k];
            }

            /// Set up 162-168 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyzz_xx = cbuffer.data(gd_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_yyzz_xy = cbuffer.data(gd_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_yyzz_xz = cbuffer.data(gd_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_yyzz_yy = cbuffer.data(gd_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_yyzz_yz = cbuffer.data(gd_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_yyzz_zz = cbuffer.data(gd_geom_01_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyz_xx, g_0_y_yyz_xxz, g_0_y_yyz_xy, g_0_y_yyz_xyz, g_0_y_yyz_xz, g_0_y_yyz_xzz, g_0_y_yyz_yy, g_0_y_yyz_yyz, g_0_y_yyz_yz, g_0_y_yyz_yzz, g_0_y_yyz_zz, g_0_y_yyz_zzz, g_0_y_yyzz_xx, g_0_y_yyzz_xy, g_0_y_yyzz_xz, g_0_y_yyzz_yy, g_0_y_yyzz_yz, g_0_y_yyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyzz_xx[k] = -g_0_y_yyz_xx[k] * ab_z + g_0_y_yyz_xxz[k];

                g_0_y_yyzz_xy[k] = -g_0_y_yyz_xy[k] * ab_z + g_0_y_yyz_xyz[k];

                g_0_y_yyzz_xz[k] = -g_0_y_yyz_xz[k] * ab_z + g_0_y_yyz_xzz[k];

                g_0_y_yyzz_yy[k] = -g_0_y_yyz_yy[k] * ab_z + g_0_y_yyz_yyz[k];

                g_0_y_yyzz_yz[k] = -g_0_y_yyz_yz[k] * ab_z + g_0_y_yyz_yzz[k];

                g_0_y_yyzz_zz[k] = -g_0_y_yyz_zz[k] * ab_z + g_0_y_yyz_zzz[k];
            }

            /// Set up 168-174 components of targeted buffer : cbuffer.data(

            auto g_0_y_yzzz_xx = cbuffer.data(gd_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_yzzz_xy = cbuffer.data(gd_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_yzzz_xz = cbuffer.data(gd_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_yzzz_yy = cbuffer.data(gd_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_yzzz_yz = cbuffer.data(gd_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_yzzz_zz = cbuffer.data(gd_geom_01_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yzz_xx, g_0_y_yzz_xxz, g_0_y_yzz_xy, g_0_y_yzz_xyz, g_0_y_yzz_xz, g_0_y_yzz_xzz, g_0_y_yzz_yy, g_0_y_yzz_yyz, g_0_y_yzz_yz, g_0_y_yzz_yzz, g_0_y_yzz_zz, g_0_y_yzz_zzz, g_0_y_yzzz_xx, g_0_y_yzzz_xy, g_0_y_yzzz_xz, g_0_y_yzzz_yy, g_0_y_yzzz_yz, g_0_y_yzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzzz_xx[k] = -g_0_y_yzz_xx[k] * ab_z + g_0_y_yzz_xxz[k];

                g_0_y_yzzz_xy[k] = -g_0_y_yzz_xy[k] * ab_z + g_0_y_yzz_xyz[k];

                g_0_y_yzzz_xz[k] = -g_0_y_yzz_xz[k] * ab_z + g_0_y_yzz_xzz[k];

                g_0_y_yzzz_yy[k] = -g_0_y_yzz_yy[k] * ab_z + g_0_y_yzz_yyz[k];

                g_0_y_yzzz_yz[k] = -g_0_y_yzz_yz[k] * ab_z + g_0_y_yzz_yzz[k];

                g_0_y_yzzz_zz[k] = -g_0_y_yzz_zz[k] * ab_z + g_0_y_yzz_zzz[k];
            }

            /// Set up 174-180 components of targeted buffer : cbuffer.data(

            auto g_0_y_zzzz_xx = cbuffer.data(gd_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_zzzz_xy = cbuffer.data(gd_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_zzzz_xz = cbuffer.data(gd_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_zzzz_yy = cbuffer.data(gd_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_zzzz_yz = cbuffer.data(gd_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_zzzz_zz = cbuffer.data(gd_geom_01_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzz_xx, g_0_y_zzz_xxz, g_0_y_zzz_xy, g_0_y_zzz_xyz, g_0_y_zzz_xz, g_0_y_zzz_xzz, g_0_y_zzz_yy, g_0_y_zzz_yyz, g_0_y_zzz_yz, g_0_y_zzz_yzz, g_0_y_zzz_zz, g_0_y_zzz_zzz, g_0_y_zzzz_xx, g_0_y_zzzz_xy, g_0_y_zzzz_xz, g_0_y_zzzz_yy, g_0_y_zzzz_yz, g_0_y_zzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzzz_xx[k] = -g_0_y_zzz_xx[k] * ab_z + g_0_y_zzz_xxz[k];

                g_0_y_zzzz_xy[k] = -g_0_y_zzz_xy[k] * ab_z + g_0_y_zzz_xyz[k];

                g_0_y_zzzz_xz[k] = -g_0_y_zzz_xz[k] * ab_z + g_0_y_zzz_xzz[k];

                g_0_y_zzzz_yy[k] = -g_0_y_zzz_yy[k] * ab_z + g_0_y_zzz_yyz[k];

                g_0_y_zzzz_yz[k] = -g_0_y_zzz_yz[k] * ab_z + g_0_y_zzz_yzz[k];

                g_0_y_zzzz_zz[k] = -g_0_y_zzz_zz[k] * ab_z + g_0_y_zzz_zzz[k];
            }

            /// Set up 180-186 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxx_xx = cbuffer.data(gd_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_z_xxxx_xy = cbuffer.data(gd_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_z_xxxx_xz = cbuffer.data(gd_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_z_xxxx_yy = cbuffer.data(gd_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_z_xxxx_yz = cbuffer.data(gd_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_z_xxxx_zz = cbuffer.data(gd_geom_01_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxx_xx, g_0_z_xxx_xxx, g_0_z_xxx_xxy, g_0_z_xxx_xxz, g_0_z_xxx_xy, g_0_z_xxx_xyy, g_0_z_xxx_xyz, g_0_z_xxx_xz, g_0_z_xxx_xzz, g_0_z_xxx_yy, g_0_z_xxx_yz, g_0_z_xxx_zz, g_0_z_xxxx_xx, g_0_z_xxxx_xy, g_0_z_xxxx_xz, g_0_z_xxxx_yy, g_0_z_xxxx_yz, g_0_z_xxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxx_xx[k] = -g_0_z_xxx_xx[k] * ab_x + g_0_z_xxx_xxx[k];

                g_0_z_xxxx_xy[k] = -g_0_z_xxx_xy[k] * ab_x + g_0_z_xxx_xxy[k];

                g_0_z_xxxx_xz[k] = -g_0_z_xxx_xz[k] * ab_x + g_0_z_xxx_xxz[k];

                g_0_z_xxxx_yy[k] = -g_0_z_xxx_yy[k] * ab_x + g_0_z_xxx_xyy[k];

                g_0_z_xxxx_yz[k] = -g_0_z_xxx_yz[k] * ab_x + g_0_z_xxx_xyz[k];

                g_0_z_xxxx_zz[k] = -g_0_z_xxx_zz[k] * ab_x + g_0_z_xxx_xzz[k];
            }

            /// Set up 186-192 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxy_xx = cbuffer.data(gd_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_z_xxxy_xy = cbuffer.data(gd_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_z_xxxy_xz = cbuffer.data(gd_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_z_xxxy_yy = cbuffer.data(gd_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_z_xxxy_yz = cbuffer.data(gd_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_z_xxxy_zz = cbuffer.data(gd_geom_01_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxy_xx, g_0_z_xxxy_xy, g_0_z_xxxy_xz, g_0_z_xxxy_yy, g_0_z_xxxy_yz, g_0_z_xxxy_zz, g_0_z_xxy_xx, g_0_z_xxy_xxx, g_0_z_xxy_xxy, g_0_z_xxy_xxz, g_0_z_xxy_xy, g_0_z_xxy_xyy, g_0_z_xxy_xyz, g_0_z_xxy_xz, g_0_z_xxy_xzz, g_0_z_xxy_yy, g_0_z_xxy_yz, g_0_z_xxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxy_xx[k] = -g_0_z_xxy_xx[k] * ab_x + g_0_z_xxy_xxx[k];

                g_0_z_xxxy_xy[k] = -g_0_z_xxy_xy[k] * ab_x + g_0_z_xxy_xxy[k];

                g_0_z_xxxy_xz[k] = -g_0_z_xxy_xz[k] * ab_x + g_0_z_xxy_xxz[k];

                g_0_z_xxxy_yy[k] = -g_0_z_xxy_yy[k] * ab_x + g_0_z_xxy_xyy[k];

                g_0_z_xxxy_yz[k] = -g_0_z_xxy_yz[k] * ab_x + g_0_z_xxy_xyz[k];

                g_0_z_xxxy_zz[k] = -g_0_z_xxy_zz[k] * ab_x + g_0_z_xxy_xzz[k];
            }

            /// Set up 192-198 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxz_xx = cbuffer.data(gd_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_z_xxxz_xy = cbuffer.data(gd_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_z_xxxz_xz = cbuffer.data(gd_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_z_xxxz_yy = cbuffer.data(gd_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_z_xxxz_yz = cbuffer.data(gd_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_z_xxxz_zz = cbuffer.data(gd_geom_01_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxz_xx, g_0_z_xxxz_xy, g_0_z_xxxz_xz, g_0_z_xxxz_yy, g_0_z_xxxz_yz, g_0_z_xxxz_zz, g_0_z_xxz_xx, g_0_z_xxz_xxx, g_0_z_xxz_xxy, g_0_z_xxz_xxz, g_0_z_xxz_xy, g_0_z_xxz_xyy, g_0_z_xxz_xyz, g_0_z_xxz_xz, g_0_z_xxz_xzz, g_0_z_xxz_yy, g_0_z_xxz_yz, g_0_z_xxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxz_xx[k] = -g_0_z_xxz_xx[k] * ab_x + g_0_z_xxz_xxx[k];

                g_0_z_xxxz_xy[k] = -g_0_z_xxz_xy[k] * ab_x + g_0_z_xxz_xxy[k];

                g_0_z_xxxz_xz[k] = -g_0_z_xxz_xz[k] * ab_x + g_0_z_xxz_xxz[k];

                g_0_z_xxxz_yy[k] = -g_0_z_xxz_yy[k] * ab_x + g_0_z_xxz_xyy[k];

                g_0_z_xxxz_yz[k] = -g_0_z_xxz_yz[k] * ab_x + g_0_z_xxz_xyz[k];

                g_0_z_xxxz_zz[k] = -g_0_z_xxz_zz[k] * ab_x + g_0_z_xxz_xzz[k];
            }

            /// Set up 198-204 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyy_xx = cbuffer.data(gd_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_z_xxyy_xy = cbuffer.data(gd_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_z_xxyy_xz = cbuffer.data(gd_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_z_xxyy_yy = cbuffer.data(gd_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_z_xxyy_yz = cbuffer.data(gd_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_z_xxyy_zz = cbuffer.data(gd_geom_01_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyy_xx, g_0_z_xxyy_xy, g_0_z_xxyy_xz, g_0_z_xxyy_yy, g_0_z_xxyy_yz, g_0_z_xxyy_zz, g_0_z_xyy_xx, g_0_z_xyy_xxx, g_0_z_xyy_xxy, g_0_z_xyy_xxz, g_0_z_xyy_xy, g_0_z_xyy_xyy, g_0_z_xyy_xyz, g_0_z_xyy_xz, g_0_z_xyy_xzz, g_0_z_xyy_yy, g_0_z_xyy_yz, g_0_z_xyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyy_xx[k] = -g_0_z_xyy_xx[k] * ab_x + g_0_z_xyy_xxx[k];

                g_0_z_xxyy_xy[k] = -g_0_z_xyy_xy[k] * ab_x + g_0_z_xyy_xxy[k];

                g_0_z_xxyy_xz[k] = -g_0_z_xyy_xz[k] * ab_x + g_0_z_xyy_xxz[k];

                g_0_z_xxyy_yy[k] = -g_0_z_xyy_yy[k] * ab_x + g_0_z_xyy_xyy[k];

                g_0_z_xxyy_yz[k] = -g_0_z_xyy_yz[k] * ab_x + g_0_z_xyy_xyz[k];

                g_0_z_xxyy_zz[k] = -g_0_z_xyy_zz[k] * ab_x + g_0_z_xyy_xzz[k];
            }

            /// Set up 204-210 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyz_xx = cbuffer.data(gd_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_z_xxyz_xy = cbuffer.data(gd_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_z_xxyz_xz = cbuffer.data(gd_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_z_xxyz_yy = cbuffer.data(gd_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_z_xxyz_yz = cbuffer.data(gd_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_z_xxyz_zz = cbuffer.data(gd_geom_01_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyz_xx, g_0_z_xxyz_xy, g_0_z_xxyz_xz, g_0_z_xxyz_yy, g_0_z_xxyz_yz, g_0_z_xxyz_zz, g_0_z_xyz_xx, g_0_z_xyz_xxx, g_0_z_xyz_xxy, g_0_z_xyz_xxz, g_0_z_xyz_xy, g_0_z_xyz_xyy, g_0_z_xyz_xyz, g_0_z_xyz_xz, g_0_z_xyz_xzz, g_0_z_xyz_yy, g_0_z_xyz_yz, g_0_z_xyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyz_xx[k] = -g_0_z_xyz_xx[k] * ab_x + g_0_z_xyz_xxx[k];

                g_0_z_xxyz_xy[k] = -g_0_z_xyz_xy[k] * ab_x + g_0_z_xyz_xxy[k];

                g_0_z_xxyz_xz[k] = -g_0_z_xyz_xz[k] * ab_x + g_0_z_xyz_xxz[k];

                g_0_z_xxyz_yy[k] = -g_0_z_xyz_yy[k] * ab_x + g_0_z_xyz_xyy[k];

                g_0_z_xxyz_yz[k] = -g_0_z_xyz_yz[k] * ab_x + g_0_z_xyz_xyz[k];

                g_0_z_xxyz_zz[k] = -g_0_z_xyz_zz[k] * ab_x + g_0_z_xyz_xzz[k];
            }

            /// Set up 210-216 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxzz_xx = cbuffer.data(gd_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_z_xxzz_xy = cbuffer.data(gd_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_z_xxzz_xz = cbuffer.data(gd_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_z_xxzz_yy = cbuffer.data(gd_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_z_xxzz_yz = cbuffer.data(gd_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_z_xxzz_zz = cbuffer.data(gd_geom_01_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxzz_xx, g_0_z_xxzz_xy, g_0_z_xxzz_xz, g_0_z_xxzz_yy, g_0_z_xxzz_yz, g_0_z_xxzz_zz, g_0_z_xzz_xx, g_0_z_xzz_xxx, g_0_z_xzz_xxy, g_0_z_xzz_xxz, g_0_z_xzz_xy, g_0_z_xzz_xyy, g_0_z_xzz_xyz, g_0_z_xzz_xz, g_0_z_xzz_xzz, g_0_z_xzz_yy, g_0_z_xzz_yz, g_0_z_xzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxzz_xx[k] = -g_0_z_xzz_xx[k] * ab_x + g_0_z_xzz_xxx[k];

                g_0_z_xxzz_xy[k] = -g_0_z_xzz_xy[k] * ab_x + g_0_z_xzz_xxy[k];

                g_0_z_xxzz_xz[k] = -g_0_z_xzz_xz[k] * ab_x + g_0_z_xzz_xxz[k];

                g_0_z_xxzz_yy[k] = -g_0_z_xzz_yy[k] * ab_x + g_0_z_xzz_xyy[k];

                g_0_z_xxzz_yz[k] = -g_0_z_xzz_yz[k] * ab_x + g_0_z_xzz_xyz[k];

                g_0_z_xxzz_zz[k] = -g_0_z_xzz_zz[k] * ab_x + g_0_z_xzz_xzz[k];
            }

            /// Set up 216-222 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyy_xx = cbuffer.data(gd_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_z_xyyy_xy = cbuffer.data(gd_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_z_xyyy_xz = cbuffer.data(gd_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_z_xyyy_yy = cbuffer.data(gd_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_z_xyyy_yz = cbuffer.data(gd_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_z_xyyy_zz = cbuffer.data(gd_geom_01_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyy_xx, g_0_z_xyyy_xy, g_0_z_xyyy_xz, g_0_z_xyyy_yy, g_0_z_xyyy_yz, g_0_z_xyyy_zz, g_0_z_yyy_xx, g_0_z_yyy_xxx, g_0_z_yyy_xxy, g_0_z_yyy_xxz, g_0_z_yyy_xy, g_0_z_yyy_xyy, g_0_z_yyy_xyz, g_0_z_yyy_xz, g_0_z_yyy_xzz, g_0_z_yyy_yy, g_0_z_yyy_yz, g_0_z_yyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyy_xx[k] = -g_0_z_yyy_xx[k] * ab_x + g_0_z_yyy_xxx[k];

                g_0_z_xyyy_xy[k] = -g_0_z_yyy_xy[k] * ab_x + g_0_z_yyy_xxy[k];

                g_0_z_xyyy_xz[k] = -g_0_z_yyy_xz[k] * ab_x + g_0_z_yyy_xxz[k];

                g_0_z_xyyy_yy[k] = -g_0_z_yyy_yy[k] * ab_x + g_0_z_yyy_xyy[k];

                g_0_z_xyyy_yz[k] = -g_0_z_yyy_yz[k] * ab_x + g_0_z_yyy_xyz[k];

                g_0_z_xyyy_zz[k] = -g_0_z_yyy_zz[k] * ab_x + g_0_z_yyy_xzz[k];
            }

            /// Set up 222-228 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyz_xx = cbuffer.data(gd_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_z_xyyz_xy = cbuffer.data(gd_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_z_xyyz_xz = cbuffer.data(gd_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_z_xyyz_yy = cbuffer.data(gd_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_z_xyyz_yz = cbuffer.data(gd_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_z_xyyz_zz = cbuffer.data(gd_geom_01_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyz_xx, g_0_z_xyyz_xy, g_0_z_xyyz_xz, g_0_z_xyyz_yy, g_0_z_xyyz_yz, g_0_z_xyyz_zz, g_0_z_yyz_xx, g_0_z_yyz_xxx, g_0_z_yyz_xxy, g_0_z_yyz_xxz, g_0_z_yyz_xy, g_0_z_yyz_xyy, g_0_z_yyz_xyz, g_0_z_yyz_xz, g_0_z_yyz_xzz, g_0_z_yyz_yy, g_0_z_yyz_yz, g_0_z_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyz_xx[k] = -g_0_z_yyz_xx[k] * ab_x + g_0_z_yyz_xxx[k];

                g_0_z_xyyz_xy[k] = -g_0_z_yyz_xy[k] * ab_x + g_0_z_yyz_xxy[k];

                g_0_z_xyyz_xz[k] = -g_0_z_yyz_xz[k] * ab_x + g_0_z_yyz_xxz[k];

                g_0_z_xyyz_yy[k] = -g_0_z_yyz_yy[k] * ab_x + g_0_z_yyz_xyy[k];

                g_0_z_xyyz_yz[k] = -g_0_z_yyz_yz[k] * ab_x + g_0_z_yyz_xyz[k];

                g_0_z_xyyz_zz[k] = -g_0_z_yyz_zz[k] * ab_x + g_0_z_yyz_xzz[k];
            }

            /// Set up 228-234 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyzz_xx = cbuffer.data(gd_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_z_xyzz_xy = cbuffer.data(gd_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_z_xyzz_xz = cbuffer.data(gd_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_z_xyzz_yy = cbuffer.data(gd_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_z_xyzz_yz = cbuffer.data(gd_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_z_xyzz_zz = cbuffer.data(gd_geom_01_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyzz_xx, g_0_z_xyzz_xy, g_0_z_xyzz_xz, g_0_z_xyzz_yy, g_0_z_xyzz_yz, g_0_z_xyzz_zz, g_0_z_yzz_xx, g_0_z_yzz_xxx, g_0_z_yzz_xxy, g_0_z_yzz_xxz, g_0_z_yzz_xy, g_0_z_yzz_xyy, g_0_z_yzz_xyz, g_0_z_yzz_xz, g_0_z_yzz_xzz, g_0_z_yzz_yy, g_0_z_yzz_yz, g_0_z_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyzz_xx[k] = -g_0_z_yzz_xx[k] * ab_x + g_0_z_yzz_xxx[k];

                g_0_z_xyzz_xy[k] = -g_0_z_yzz_xy[k] * ab_x + g_0_z_yzz_xxy[k];

                g_0_z_xyzz_xz[k] = -g_0_z_yzz_xz[k] * ab_x + g_0_z_yzz_xxz[k];

                g_0_z_xyzz_yy[k] = -g_0_z_yzz_yy[k] * ab_x + g_0_z_yzz_xyy[k];

                g_0_z_xyzz_yz[k] = -g_0_z_yzz_yz[k] * ab_x + g_0_z_yzz_xyz[k];

                g_0_z_xyzz_zz[k] = -g_0_z_yzz_zz[k] * ab_x + g_0_z_yzz_xzz[k];
            }

            /// Set up 234-240 components of targeted buffer : cbuffer.data(

            auto g_0_z_xzzz_xx = cbuffer.data(gd_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_z_xzzz_xy = cbuffer.data(gd_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_z_xzzz_xz = cbuffer.data(gd_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_z_xzzz_yy = cbuffer.data(gd_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_z_xzzz_yz = cbuffer.data(gd_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_z_xzzz_zz = cbuffer.data(gd_geom_01_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzzz_xx, g_0_z_xzzz_xy, g_0_z_xzzz_xz, g_0_z_xzzz_yy, g_0_z_xzzz_yz, g_0_z_xzzz_zz, g_0_z_zzz_xx, g_0_z_zzz_xxx, g_0_z_zzz_xxy, g_0_z_zzz_xxz, g_0_z_zzz_xy, g_0_z_zzz_xyy, g_0_z_zzz_xyz, g_0_z_zzz_xz, g_0_z_zzz_xzz, g_0_z_zzz_yy, g_0_z_zzz_yz, g_0_z_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzzz_xx[k] = -g_0_z_zzz_xx[k] * ab_x + g_0_z_zzz_xxx[k];

                g_0_z_xzzz_xy[k] = -g_0_z_zzz_xy[k] * ab_x + g_0_z_zzz_xxy[k];

                g_0_z_xzzz_xz[k] = -g_0_z_zzz_xz[k] * ab_x + g_0_z_zzz_xxz[k];

                g_0_z_xzzz_yy[k] = -g_0_z_zzz_yy[k] * ab_x + g_0_z_zzz_xyy[k];

                g_0_z_xzzz_yz[k] = -g_0_z_zzz_yz[k] * ab_x + g_0_z_zzz_xyz[k];

                g_0_z_xzzz_zz[k] = -g_0_z_zzz_zz[k] * ab_x + g_0_z_zzz_xzz[k];
            }

            /// Set up 240-246 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyy_xx = cbuffer.data(gd_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_z_yyyy_xy = cbuffer.data(gd_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_z_yyyy_xz = cbuffer.data(gd_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_z_yyyy_yy = cbuffer.data(gd_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_z_yyyy_yz = cbuffer.data(gd_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_z_yyyy_zz = cbuffer.data(gd_geom_01_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyy_xx, g_0_z_yyy_xxy, g_0_z_yyy_xy, g_0_z_yyy_xyy, g_0_z_yyy_xyz, g_0_z_yyy_xz, g_0_z_yyy_yy, g_0_z_yyy_yyy, g_0_z_yyy_yyz, g_0_z_yyy_yz, g_0_z_yyy_yzz, g_0_z_yyy_zz, g_0_z_yyyy_xx, g_0_z_yyyy_xy, g_0_z_yyyy_xz, g_0_z_yyyy_yy, g_0_z_yyyy_yz, g_0_z_yyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyy_xx[k] = -g_0_z_yyy_xx[k] * ab_y + g_0_z_yyy_xxy[k];

                g_0_z_yyyy_xy[k] = -g_0_z_yyy_xy[k] * ab_y + g_0_z_yyy_xyy[k];

                g_0_z_yyyy_xz[k] = -g_0_z_yyy_xz[k] * ab_y + g_0_z_yyy_xyz[k];

                g_0_z_yyyy_yy[k] = -g_0_z_yyy_yy[k] * ab_y + g_0_z_yyy_yyy[k];

                g_0_z_yyyy_yz[k] = -g_0_z_yyy_yz[k] * ab_y + g_0_z_yyy_yyz[k];

                g_0_z_yyyy_zz[k] = -g_0_z_yyy_zz[k] * ab_y + g_0_z_yyy_yzz[k];
            }

            /// Set up 246-252 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyz_xx = cbuffer.data(gd_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_z_yyyz_xy = cbuffer.data(gd_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_z_yyyz_xz = cbuffer.data(gd_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_z_yyyz_yy = cbuffer.data(gd_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_z_yyyz_yz = cbuffer.data(gd_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_z_yyyz_zz = cbuffer.data(gd_geom_01_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyz_xx, g_0_z_yyyz_xy, g_0_z_yyyz_xz, g_0_z_yyyz_yy, g_0_z_yyyz_yz, g_0_z_yyyz_zz, g_0_z_yyz_xx, g_0_z_yyz_xxy, g_0_z_yyz_xy, g_0_z_yyz_xyy, g_0_z_yyz_xyz, g_0_z_yyz_xz, g_0_z_yyz_yy, g_0_z_yyz_yyy, g_0_z_yyz_yyz, g_0_z_yyz_yz, g_0_z_yyz_yzz, g_0_z_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyz_xx[k] = -g_0_z_yyz_xx[k] * ab_y + g_0_z_yyz_xxy[k];

                g_0_z_yyyz_xy[k] = -g_0_z_yyz_xy[k] * ab_y + g_0_z_yyz_xyy[k];

                g_0_z_yyyz_xz[k] = -g_0_z_yyz_xz[k] * ab_y + g_0_z_yyz_xyz[k];

                g_0_z_yyyz_yy[k] = -g_0_z_yyz_yy[k] * ab_y + g_0_z_yyz_yyy[k];

                g_0_z_yyyz_yz[k] = -g_0_z_yyz_yz[k] * ab_y + g_0_z_yyz_yyz[k];

                g_0_z_yyyz_zz[k] = -g_0_z_yyz_zz[k] * ab_y + g_0_z_yyz_yzz[k];
            }

            /// Set up 252-258 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyzz_xx = cbuffer.data(gd_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_z_yyzz_xy = cbuffer.data(gd_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_z_yyzz_xz = cbuffer.data(gd_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_z_yyzz_yy = cbuffer.data(gd_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_z_yyzz_yz = cbuffer.data(gd_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_z_yyzz_zz = cbuffer.data(gd_geom_01_off + 257 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyzz_xx, g_0_z_yyzz_xy, g_0_z_yyzz_xz, g_0_z_yyzz_yy, g_0_z_yyzz_yz, g_0_z_yyzz_zz, g_0_z_yzz_xx, g_0_z_yzz_xxy, g_0_z_yzz_xy, g_0_z_yzz_xyy, g_0_z_yzz_xyz, g_0_z_yzz_xz, g_0_z_yzz_yy, g_0_z_yzz_yyy, g_0_z_yzz_yyz, g_0_z_yzz_yz, g_0_z_yzz_yzz, g_0_z_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyzz_xx[k] = -g_0_z_yzz_xx[k] * ab_y + g_0_z_yzz_xxy[k];

                g_0_z_yyzz_xy[k] = -g_0_z_yzz_xy[k] * ab_y + g_0_z_yzz_xyy[k];

                g_0_z_yyzz_xz[k] = -g_0_z_yzz_xz[k] * ab_y + g_0_z_yzz_xyz[k];

                g_0_z_yyzz_yy[k] = -g_0_z_yzz_yy[k] * ab_y + g_0_z_yzz_yyy[k];

                g_0_z_yyzz_yz[k] = -g_0_z_yzz_yz[k] * ab_y + g_0_z_yzz_yyz[k];

                g_0_z_yyzz_zz[k] = -g_0_z_yzz_zz[k] * ab_y + g_0_z_yzz_yzz[k];
            }

            /// Set up 258-264 components of targeted buffer : cbuffer.data(

            auto g_0_z_yzzz_xx = cbuffer.data(gd_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_z_yzzz_xy = cbuffer.data(gd_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_z_yzzz_xz = cbuffer.data(gd_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_z_yzzz_yy = cbuffer.data(gd_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_z_yzzz_yz = cbuffer.data(gd_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_z_yzzz_zz = cbuffer.data(gd_geom_01_off + 263 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzzz_xx, g_0_z_yzzz_xy, g_0_z_yzzz_xz, g_0_z_yzzz_yy, g_0_z_yzzz_yz, g_0_z_yzzz_zz, g_0_z_zzz_xx, g_0_z_zzz_xxy, g_0_z_zzz_xy, g_0_z_zzz_xyy, g_0_z_zzz_xyz, g_0_z_zzz_xz, g_0_z_zzz_yy, g_0_z_zzz_yyy, g_0_z_zzz_yyz, g_0_z_zzz_yz, g_0_z_zzz_yzz, g_0_z_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzzz_xx[k] = -g_0_z_zzz_xx[k] * ab_y + g_0_z_zzz_xxy[k];

                g_0_z_yzzz_xy[k] = -g_0_z_zzz_xy[k] * ab_y + g_0_z_zzz_xyy[k];

                g_0_z_yzzz_xz[k] = -g_0_z_zzz_xz[k] * ab_y + g_0_z_zzz_xyz[k];

                g_0_z_yzzz_yy[k] = -g_0_z_zzz_yy[k] * ab_y + g_0_z_zzz_yyy[k];

                g_0_z_yzzz_yz[k] = -g_0_z_zzz_yz[k] * ab_y + g_0_z_zzz_yyz[k];

                g_0_z_yzzz_zz[k] = -g_0_z_zzz_zz[k] * ab_y + g_0_z_zzz_yzz[k];
            }

            /// Set up 264-270 components of targeted buffer : cbuffer.data(

            auto g_0_z_zzzz_xx = cbuffer.data(gd_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_z_zzzz_xy = cbuffer.data(gd_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_z_zzzz_xz = cbuffer.data(gd_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_z_zzzz_yy = cbuffer.data(gd_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_z_zzzz_yz = cbuffer.data(gd_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_z_zzzz_zz = cbuffer.data(gd_geom_01_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzz_xx, g_0_z_zzz_xxz, g_0_z_zzz_xy, g_0_z_zzz_xyz, g_0_z_zzz_xz, g_0_z_zzz_xzz, g_0_z_zzz_yy, g_0_z_zzz_yyz, g_0_z_zzz_yz, g_0_z_zzz_yzz, g_0_z_zzz_zz, g_0_z_zzz_zzz, g_0_z_zzzz_xx, g_0_z_zzzz_xy, g_0_z_zzzz_xz, g_0_z_zzzz_yy, g_0_z_zzzz_yz, g_0_z_zzzz_zz, g_zzz_xx, g_zzz_xy, g_zzz_xz, g_zzz_yy, g_zzz_yz, g_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzzz_xx[k] = g_zzz_xx[k] - g_0_z_zzz_xx[k] * ab_z + g_0_z_zzz_xxz[k];

                g_0_z_zzzz_xy[k] = g_zzz_xy[k] - g_0_z_zzz_xy[k] * ab_z + g_0_z_zzz_xyz[k];

                g_0_z_zzzz_xz[k] = g_zzz_xz[k] - g_0_z_zzz_xz[k] * ab_z + g_0_z_zzz_xzz[k];

                g_0_z_zzzz_yy[k] = g_zzz_yy[k] - g_0_z_zzz_yy[k] * ab_z + g_0_z_zzz_yyz[k];

                g_0_z_zzzz_yz[k] = g_zzz_yz[k] - g_0_z_zzz_yz[k] * ab_z + g_0_z_zzz_yzz[k];

                g_0_z_zzzz_zz[k] = g_zzz_zz[k] - g_0_z_zzz_zz[k] * ab_z + g_0_z_zzz_zzz[k];
            }
        }
    }
}

} // erirec namespace

