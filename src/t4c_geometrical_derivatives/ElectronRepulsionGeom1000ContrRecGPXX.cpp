#include "ElectronRepulsionGeom1000ContrRecGPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_gpxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_gpxx,
                                            const size_t idx_fpxx,
                                            const size_t idx_geom_10_fpxx,
                                            const size_t idx_geom_10_fdxx,
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
            /// Set up components of auxilary buffer : FPSS

            const auto fp_off = idx_fpxx + i * dcomps + j;

            auto g_xxx_x = cbuffer.data(fp_off + 0 * ccomps * dcomps);

            auto g_xxx_y = cbuffer.data(fp_off + 1 * ccomps * dcomps);

            auto g_xxx_z = cbuffer.data(fp_off + 2 * ccomps * dcomps);

            auto g_yyy_x = cbuffer.data(fp_off + 18 * ccomps * dcomps);

            auto g_yyy_y = cbuffer.data(fp_off + 19 * ccomps * dcomps);

            auto g_yyy_z = cbuffer.data(fp_off + 20 * ccomps * dcomps);

            auto g_zzz_x = cbuffer.data(fp_off + 27 * ccomps * dcomps);

            auto g_zzz_y = cbuffer.data(fp_off + 28 * ccomps * dcomps);

            auto g_zzz_z = cbuffer.data(fp_off + 29 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FPSS

            const auto fp_geom_10_off = idx_geom_10_fpxx + i * dcomps + j;

            auto g_x_0_xxx_x = cbuffer.data(fp_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxx_y = cbuffer.data(fp_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxx_z = cbuffer.data(fp_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxy_x = cbuffer.data(fp_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxy_y = cbuffer.data(fp_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxy_z = cbuffer.data(fp_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxz_x = cbuffer.data(fp_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxz_y = cbuffer.data(fp_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxz_z = cbuffer.data(fp_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xyy_x = cbuffer.data(fp_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xyy_y = cbuffer.data(fp_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xyy_z = cbuffer.data(fp_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xyz_x = cbuffer.data(fp_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xyz_y = cbuffer.data(fp_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xyz_z = cbuffer.data(fp_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xzz_x = cbuffer.data(fp_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xzz_y = cbuffer.data(fp_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xzz_z = cbuffer.data(fp_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_yyy_x = cbuffer.data(fp_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_yyy_y = cbuffer.data(fp_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_yyy_z = cbuffer.data(fp_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_yyz_x = cbuffer.data(fp_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_yyz_y = cbuffer.data(fp_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_yyz_z = cbuffer.data(fp_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_yzz_x = cbuffer.data(fp_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_yzz_y = cbuffer.data(fp_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_yzz_z = cbuffer.data(fp_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_zzz_x = cbuffer.data(fp_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_zzz_y = cbuffer.data(fp_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_zzz_z = cbuffer.data(fp_geom_10_off + 29 * ccomps * dcomps);

            auto g_y_0_xxx_x = cbuffer.data(fp_geom_10_off + 30 * ccomps * dcomps);

            auto g_y_0_xxx_y = cbuffer.data(fp_geom_10_off + 31 * ccomps * dcomps);

            auto g_y_0_xxx_z = cbuffer.data(fp_geom_10_off + 32 * ccomps * dcomps);

            auto g_y_0_xxy_x = cbuffer.data(fp_geom_10_off + 33 * ccomps * dcomps);

            auto g_y_0_xxy_y = cbuffer.data(fp_geom_10_off + 34 * ccomps * dcomps);

            auto g_y_0_xxy_z = cbuffer.data(fp_geom_10_off + 35 * ccomps * dcomps);

            auto g_y_0_xxz_x = cbuffer.data(fp_geom_10_off + 36 * ccomps * dcomps);

            auto g_y_0_xxz_y = cbuffer.data(fp_geom_10_off + 37 * ccomps * dcomps);

            auto g_y_0_xxz_z = cbuffer.data(fp_geom_10_off + 38 * ccomps * dcomps);

            auto g_y_0_xyy_x = cbuffer.data(fp_geom_10_off + 39 * ccomps * dcomps);

            auto g_y_0_xyy_y = cbuffer.data(fp_geom_10_off + 40 * ccomps * dcomps);

            auto g_y_0_xyy_z = cbuffer.data(fp_geom_10_off + 41 * ccomps * dcomps);

            auto g_y_0_xyz_x = cbuffer.data(fp_geom_10_off + 42 * ccomps * dcomps);

            auto g_y_0_xyz_y = cbuffer.data(fp_geom_10_off + 43 * ccomps * dcomps);

            auto g_y_0_xyz_z = cbuffer.data(fp_geom_10_off + 44 * ccomps * dcomps);

            auto g_y_0_xzz_x = cbuffer.data(fp_geom_10_off + 45 * ccomps * dcomps);

            auto g_y_0_xzz_y = cbuffer.data(fp_geom_10_off + 46 * ccomps * dcomps);

            auto g_y_0_xzz_z = cbuffer.data(fp_geom_10_off + 47 * ccomps * dcomps);

            auto g_y_0_yyy_x = cbuffer.data(fp_geom_10_off + 48 * ccomps * dcomps);

            auto g_y_0_yyy_y = cbuffer.data(fp_geom_10_off + 49 * ccomps * dcomps);

            auto g_y_0_yyy_z = cbuffer.data(fp_geom_10_off + 50 * ccomps * dcomps);

            auto g_y_0_yyz_x = cbuffer.data(fp_geom_10_off + 51 * ccomps * dcomps);

            auto g_y_0_yyz_y = cbuffer.data(fp_geom_10_off + 52 * ccomps * dcomps);

            auto g_y_0_yyz_z = cbuffer.data(fp_geom_10_off + 53 * ccomps * dcomps);

            auto g_y_0_yzz_x = cbuffer.data(fp_geom_10_off + 54 * ccomps * dcomps);

            auto g_y_0_yzz_y = cbuffer.data(fp_geom_10_off + 55 * ccomps * dcomps);

            auto g_y_0_yzz_z = cbuffer.data(fp_geom_10_off + 56 * ccomps * dcomps);

            auto g_y_0_zzz_x = cbuffer.data(fp_geom_10_off + 57 * ccomps * dcomps);

            auto g_y_0_zzz_y = cbuffer.data(fp_geom_10_off + 58 * ccomps * dcomps);

            auto g_y_0_zzz_z = cbuffer.data(fp_geom_10_off + 59 * ccomps * dcomps);

            auto g_z_0_xxx_x = cbuffer.data(fp_geom_10_off + 60 * ccomps * dcomps);

            auto g_z_0_xxx_y = cbuffer.data(fp_geom_10_off + 61 * ccomps * dcomps);

            auto g_z_0_xxx_z = cbuffer.data(fp_geom_10_off + 62 * ccomps * dcomps);

            auto g_z_0_xxy_x = cbuffer.data(fp_geom_10_off + 63 * ccomps * dcomps);

            auto g_z_0_xxy_y = cbuffer.data(fp_geom_10_off + 64 * ccomps * dcomps);

            auto g_z_0_xxy_z = cbuffer.data(fp_geom_10_off + 65 * ccomps * dcomps);

            auto g_z_0_xxz_x = cbuffer.data(fp_geom_10_off + 66 * ccomps * dcomps);

            auto g_z_0_xxz_y = cbuffer.data(fp_geom_10_off + 67 * ccomps * dcomps);

            auto g_z_0_xxz_z = cbuffer.data(fp_geom_10_off + 68 * ccomps * dcomps);

            auto g_z_0_xyy_x = cbuffer.data(fp_geom_10_off + 69 * ccomps * dcomps);

            auto g_z_0_xyy_y = cbuffer.data(fp_geom_10_off + 70 * ccomps * dcomps);

            auto g_z_0_xyy_z = cbuffer.data(fp_geom_10_off + 71 * ccomps * dcomps);

            auto g_z_0_xyz_x = cbuffer.data(fp_geom_10_off + 72 * ccomps * dcomps);

            auto g_z_0_xyz_y = cbuffer.data(fp_geom_10_off + 73 * ccomps * dcomps);

            auto g_z_0_xyz_z = cbuffer.data(fp_geom_10_off + 74 * ccomps * dcomps);

            auto g_z_0_xzz_x = cbuffer.data(fp_geom_10_off + 75 * ccomps * dcomps);

            auto g_z_0_xzz_y = cbuffer.data(fp_geom_10_off + 76 * ccomps * dcomps);

            auto g_z_0_xzz_z = cbuffer.data(fp_geom_10_off + 77 * ccomps * dcomps);

            auto g_z_0_yyy_x = cbuffer.data(fp_geom_10_off + 78 * ccomps * dcomps);

            auto g_z_0_yyy_y = cbuffer.data(fp_geom_10_off + 79 * ccomps * dcomps);

            auto g_z_0_yyy_z = cbuffer.data(fp_geom_10_off + 80 * ccomps * dcomps);

            auto g_z_0_yyz_x = cbuffer.data(fp_geom_10_off + 81 * ccomps * dcomps);

            auto g_z_0_yyz_y = cbuffer.data(fp_geom_10_off + 82 * ccomps * dcomps);

            auto g_z_0_yyz_z = cbuffer.data(fp_geom_10_off + 83 * ccomps * dcomps);

            auto g_z_0_yzz_x = cbuffer.data(fp_geom_10_off + 84 * ccomps * dcomps);

            auto g_z_0_yzz_y = cbuffer.data(fp_geom_10_off + 85 * ccomps * dcomps);

            auto g_z_0_yzz_z = cbuffer.data(fp_geom_10_off + 86 * ccomps * dcomps);

            auto g_z_0_zzz_x = cbuffer.data(fp_geom_10_off + 87 * ccomps * dcomps);

            auto g_z_0_zzz_y = cbuffer.data(fp_geom_10_off + 88 * ccomps * dcomps);

            auto g_z_0_zzz_z = cbuffer.data(fp_geom_10_off + 89 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FDSS

            const auto fd_geom_10_off = idx_geom_10_fdxx + i * dcomps + j;

            auto g_x_0_xxx_xx = cbuffer.data(fd_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxx_xy = cbuffer.data(fd_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxx_xz = cbuffer.data(fd_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxx_yy = cbuffer.data(fd_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxx_yz = cbuffer.data(fd_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxx_zz = cbuffer.data(fd_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxy_xy = cbuffer.data(fd_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxy_yy = cbuffer.data(fd_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxy_yz = cbuffer.data(fd_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxz_xy = cbuffer.data(fd_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxz_xz = cbuffer.data(fd_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxz_yy = cbuffer.data(fd_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxz_yz = cbuffer.data(fd_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxz_zz = cbuffer.data(fd_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xyy_xy = cbuffer.data(fd_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xyy_yy = cbuffer.data(fd_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xyy_yz = cbuffer.data(fd_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xyz_xy = cbuffer.data(fd_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xyz_yy = cbuffer.data(fd_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xyz_yz = cbuffer.data(fd_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xzz_xy = cbuffer.data(fd_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xzz_xz = cbuffer.data(fd_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xzz_yy = cbuffer.data(fd_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xzz_yz = cbuffer.data(fd_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xzz_zz = cbuffer.data(fd_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_yyy_xy = cbuffer.data(fd_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_yyy_yy = cbuffer.data(fd_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_yyy_yz = cbuffer.data(fd_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_yyz_xy = cbuffer.data(fd_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_yyz_yy = cbuffer.data(fd_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_yyz_yz = cbuffer.data(fd_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_yzz_xy = cbuffer.data(fd_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_yzz_yy = cbuffer.data(fd_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_yzz_yz = cbuffer.data(fd_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_zzz_xy = cbuffer.data(fd_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_zzz_xz = cbuffer.data(fd_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_zzz_yy = cbuffer.data(fd_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_zzz_yz = cbuffer.data(fd_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_zzz_zz = cbuffer.data(fd_geom_10_off + 59 * ccomps * dcomps);

            auto g_y_0_xxx_xx = cbuffer.data(fd_geom_10_off + 60 * ccomps * dcomps);

            auto g_y_0_xxx_xy = cbuffer.data(fd_geom_10_off + 61 * ccomps * dcomps);

            auto g_y_0_xxx_xz = cbuffer.data(fd_geom_10_off + 62 * ccomps * dcomps);

            auto g_y_0_xxy_xx = cbuffer.data(fd_geom_10_off + 66 * ccomps * dcomps);

            auto g_y_0_xxy_xy = cbuffer.data(fd_geom_10_off + 67 * ccomps * dcomps);

            auto g_y_0_xxy_xz = cbuffer.data(fd_geom_10_off + 68 * ccomps * dcomps);

            auto g_y_0_xxz_xx = cbuffer.data(fd_geom_10_off + 72 * ccomps * dcomps);

            auto g_y_0_xxz_xy = cbuffer.data(fd_geom_10_off + 73 * ccomps * dcomps);

            auto g_y_0_xxz_xz = cbuffer.data(fd_geom_10_off + 74 * ccomps * dcomps);

            auto g_y_0_xyy_xx = cbuffer.data(fd_geom_10_off + 78 * ccomps * dcomps);

            auto g_y_0_xyy_xy = cbuffer.data(fd_geom_10_off + 79 * ccomps * dcomps);

            auto g_y_0_xyy_xz = cbuffer.data(fd_geom_10_off + 80 * ccomps * dcomps);

            auto g_y_0_xyz_xx = cbuffer.data(fd_geom_10_off + 84 * ccomps * dcomps);

            auto g_y_0_xyz_xy = cbuffer.data(fd_geom_10_off + 85 * ccomps * dcomps);

            auto g_y_0_xyz_xz = cbuffer.data(fd_geom_10_off + 86 * ccomps * dcomps);

            auto g_y_0_xzz_xx = cbuffer.data(fd_geom_10_off + 90 * ccomps * dcomps);

            auto g_y_0_xzz_xy = cbuffer.data(fd_geom_10_off + 91 * ccomps * dcomps);

            auto g_y_0_xzz_xz = cbuffer.data(fd_geom_10_off + 92 * ccomps * dcomps);

            auto g_y_0_yyy_xx = cbuffer.data(fd_geom_10_off + 96 * ccomps * dcomps);

            auto g_y_0_yyy_xy = cbuffer.data(fd_geom_10_off + 97 * ccomps * dcomps);

            auto g_y_0_yyy_xz = cbuffer.data(fd_geom_10_off + 98 * ccomps * dcomps);

            auto g_y_0_yyy_yy = cbuffer.data(fd_geom_10_off + 99 * ccomps * dcomps);

            auto g_y_0_yyy_yz = cbuffer.data(fd_geom_10_off + 100 * ccomps * dcomps);

            auto g_y_0_yyy_zz = cbuffer.data(fd_geom_10_off + 101 * ccomps * dcomps);

            auto g_y_0_yyz_xx = cbuffer.data(fd_geom_10_off + 102 * ccomps * dcomps);

            auto g_y_0_yyz_xy = cbuffer.data(fd_geom_10_off + 103 * ccomps * dcomps);

            auto g_y_0_yyz_xz = cbuffer.data(fd_geom_10_off + 104 * ccomps * dcomps);

            auto g_y_0_yyz_yz = cbuffer.data(fd_geom_10_off + 106 * ccomps * dcomps);

            auto g_y_0_yyz_zz = cbuffer.data(fd_geom_10_off + 107 * ccomps * dcomps);

            auto g_y_0_yzz_xx = cbuffer.data(fd_geom_10_off + 108 * ccomps * dcomps);

            auto g_y_0_yzz_xy = cbuffer.data(fd_geom_10_off + 109 * ccomps * dcomps);

            auto g_y_0_yzz_xz = cbuffer.data(fd_geom_10_off + 110 * ccomps * dcomps);

            auto g_y_0_yzz_yz = cbuffer.data(fd_geom_10_off + 112 * ccomps * dcomps);

            auto g_y_0_yzz_zz = cbuffer.data(fd_geom_10_off + 113 * ccomps * dcomps);

            auto g_y_0_zzz_xx = cbuffer.data(fd_geom_10_off + 114 * ccomps * dcomps);

            auto g_y_0_zzz_xy = cbuffer.data(fd_geom_10_off + 115 * ccomps * dcomps);

            auto g_y_0_zzz_xz = cbuffer.data(fd_geom_10_off + 116 * ccomps * dcomps);

            auto g_y_0_zzz_yz = cbuffer.data(fd_geom_10_off + 118 * ccomps * dcomps);

            auto g_y_0_zzz_zz = cbuffer.data(fd_geom_10_off + 119 * ccomps * dcomps);

            auto g_z_0_xxx_xx = cbuffer.data(fd_geom_10_off + 120 * ccomps * dcomps);

            auto g_z_0_xxx_xy = cbuffer.data(fd_geom_10_off + 121 * ccomps * dcomps);

            auto g_z_0_xxx_xz = cbuffer.data(fd_geom_10_off + 122 * ccomps * dcomps);

            auto g_z_0_xxy_xx = cbuffer.data(fd_geom_10_off + 126 * ccomps * dcomps);

            auto g_z_0_xxy_xy = cbuffer.data(fd_geom_10_off + 127 * ccomps * dcomps);

            auto g_z_0_xxy_xz = cbuffer.data(fd_geom_10_off + 128 * ccomps * dcomps);

            auto g_z_0_xxz_xx = cbuffer.data(fd_geom_10_off + 132 * ccomps * dcomps);

            auto g_z_0_xxz_xy = cbuffer.data(fd_geom_10_off + 133 * ccomps * dcomps);

            auto g_z_0_xxz_xz = cbuffer.data(fd_geom_10_off + 134 * ccomps * dcomps);

            auto g_z_0_xyy_xx = cbuffer.data(fd_geom_10_off + 138 * ccomps * dcomps);

            auto g_z_0_xyy_xy = cbuffer.data(fd_geom_10_off + 139 * ccomps * dcomps);

            auto g_z_0_xyy_xz = cbuffer.data(fd_geom_10_off + 140 * ccomps * dcomps);

            auto g_z_0_xyz_xx = cbuffer.data(fd_geom_10_off + 144 * ccomps * dcomps);

            auto g_z_0_xyz_xy = cbuffer.data(fd_geom_10_off + 145 * ccomps * dcomps);

            auto g_z_0_xyz_xz = cbuffer.data(fd_geom_10_off + 146 * ccomps * dcomps);

            auto g_z_0_xzz_xx = cbuffer.data(fd_geom_10_off + 150 * ccomps * dcomps);

            auto g_z_0_xzz_xy = cbuffer.data(fd_geom_10_off + 151 * ccomps * dcomps);

            auto g_z_0_xzz_xz = cbuffer.data(fd_geom_10_off + 152 * ccomps * dcomps);

            auto g_z_0_yyy_xx = cbuffer.data(fd_geom_10_off + 156 * ccomps * dcomps);

            auto g_z_0_yyy_xy = cbuffer.data(fd_geom_10_off + 157 * ccomps * dcomps);

            auto g_z_0_yyy_xz = cbuffer.data(fd_geom_10_off + 158 * ccomps * dcomps);

            auto g_z_0_yyy_yy = cbuffer.data(fd_geom_10_off + 159 * ccomps * dcomps);

            auto g_z_0_yyy_yz = cbuffer.data(fd_geom_10_off + 160 * ccomps * dcomps);

            auto g_z_0_yyz_xx = cbuffer.data(fd_geom_10_off + 162 * ccomps * dcomps);

            auto g_z_0_yyz_xy = cbuffer.data(fd_geom_10_off + 163 * ccomps * dcomps);

            auto g_z_0_yyz_xz = cbuffer.data(fd_geom_10_off + 164 * ccomps * dcomps);

            auto g_z_0_yyz_yy = cbuffer.data(fd_geom_10_off + 165 * ccomps * dcomps);

            auto g_z_0_yyz_yz = cbuffer.data(fd_geom_10_off + 166 * ccomps * dcomps);

            auto g_z_0_yzz_xx = cbuffer.data(fd_geom_10_off + 168 * ccomps * dcomps);

            auto g_z_0_yzz_xy = cbuffer.data(fd_geom_10_off + 169 * ccomps * dcomps);

            auto g_z_0_yzz_xz = cbuffer.data(fd_geom_10_off + 170 * ccomps * dcomps);

            auto g_z_0_yzz_yy = cbuffer.data(fd_geom_10_off + 171 * ccomps * dcomps);

            auto g_z_0_yzz_yz = cbuffer.data(fd_geom_10_off + 172 * ccomps * dcomps);

            auto g_z_0_zzz_xx = cbuffer.data(fd_geom_10_off + 174 * ccomps * dcomps);

            auto g_z_0_zzz_xy = cbuffer.data(fd_geom_10_off + 175 * ccomps * dcomps);

            auto g_z_0_zzz_xz = cbuffer.data(fd_geom_10_off + 176 * ccomps * dcomps);

            auto g_z_0_zzz_yy = cbuffer.data(fd_geom_10_off + 177 * ccomps * dcomps);

            auto g_z_0_zzz_yz = cbuffer.data(fd_geom_10_off + 178 * ccomps * dcomps);

            auto g_z_0_zzz_zz = cbuffer.data(fd_geom_10_off + 179 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_gpxx

            const auto gp_geom_10_off = idx_geom_10_gpxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxx_x = cbuffer.data(gp_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxx_y = cbuffer.data(gp_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxx_z = cbuffer.data(gp_geom_10_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxx_x, g_x_0_xxx_xx, g_x_0_xxx_xy, g_x_0_xxx_xz, g_x_0_xxx_y, g_x_0_xxx_z, g_x_0_xxxx_x, g_x_0_xxxx_y, g_x_0_xxxx_z, g_xxx_x, g_xxx_y, g_xxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxx_x[k] = -g_xxx_x[k] - g_x_0_xxx_x[k] * ab_x + g_x_0_xxx_xx[k];

                g_x_0_xxxx_y[k] = -g_xxx_y[k] - g_x_0_xxx_y[k] * ab_x + g_x_0_xxx_xy[k];

                g_x_0_xxxx_z[k] = -g_xxx_z[k] - g_x_0_xxx_z[k] * ab_x + g_x_0_xxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxy_x = cbuffer.data(gp_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxy_y = cbuffer.data(gp_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxy_z = cbuffer.data(gp_geom_10_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxx_x, g_x_0_xxx_xy, g_x_0_xxx_y, g_x_0_xxx_yy, g_x_0_xxx_yz, g_x_0_xxx_z, g_x_0_xxxy_x, g_x_0_xxxy_y, g_x_0_xxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxy_x[k] = -g_x_0_xxx_x[k] * ab_y + g_x_0_xxx_xy[k];

                g_x_0_xxxy_y[k] = -g_x_0_xxx_y[k] * ab_y + g_x_0_xxx_yy[k];

                g_x_0_xxxy_z[k] = -g_x_0_xxx_z[k] * ab_y + g_x_0_xxx_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxz_x = cbuffer.data(gp_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxz_y = cbuffer.data(gp_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxz_z = cbuffer.data(gp_geom_10_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxx_x, g_x_0_xxx_xz, g_x_0_xxx_y, g_x_0_xxx_yz, g_x_0_xxx_z, g_x_0_xxx_zz, g_x_0_xxxz_x, g_x_0_xxxz_y, g_x_0_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxz_x[k] = -g_x_0_xxx_x[k] * ab_z + g_x_0_xxx_xz[k];

                g_x_0_xxxz_y[k] = -g_x_0_xxx_y[k] * ab_z + g_x_0_xxx_yz[k];

                g_x_0_xxxz_z[k] = -g_x_0_xxx_z[k] * ab_z + g_x_0_xxx_zz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyy_x = cbuffer.data(gp_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxyy_y = cbuffer.data(gp_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxyy_z = cbuffer.data(gp_geom_10_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxy_x, g_x_0_xxy_xy, g_x_0_xxy_y, g_x_0_xxy_yy, g_x_0_xxy_yz, g_x_0_xxy_z, g_x_0_xxyy_x, g_x_0_xxyy_y, g_x_0_xxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyy_x[k] = -g_x_0_xxy_x[k] * ab_y + g_x_0_xxy_xy[k];

                g_x_0_xxyy_y[k] = -g_x_0_xxy_y[k] * ab_y + g_x_0_xxy_yy[k];

                g_x_0_xxyy_z[k] = -g_x_0_xxy_z[k] * ab_y + g_x_0_xxy_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyz_x = cbuffer.data(gp_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxyz_y = cbuffer.data(gp_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxyz_z = cbuffer.data(gp_geom_10_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyz_x, g_x_0_xxyz_y, g_x_0_xxyz_z, g_x_0_xxz_x, g_x_0_xxz_xy, g_x_0_xxz_y, g_x_0_xxz_yy, g_x_0_xxz_yz, g_x_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyz_x[k] = -g_x_0_xxz_x[k] * ab_y + g_x_0_xxz_xy[k];

                g_x_0_xxyz_y[k] = -g_x_0_xxz_y[k] * ab_y + g_x_0_xxz_yy[k];

                g_x_0_xxyz_z[k] = -g_x_0_xxz_z[k] * ab_y + g_x_0_xxz_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzz_x = cbuffer.data(gp_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxzz_y = cbuffer.data(gp_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxzz_z = cbuffer.data(gp_geom_10_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxz_x, g_x_0_xxz_xz, g_x_0_xxz_y, g_x_0_xxz_yz, g_x_0_xxz_z, g_x_0_xxz_zz, g_x_0_xxzz_x, g_x_0_xxzz_y, g_x_0_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzz_x[k] = -g_x_0_xxz_x[k] * ab_z + g_x_0_xxz_xz[k];

                g_x_0_xxzz_y[k] = -g_x_0_xxz_y[k] * ab_z + g_x_0_xxz_yz[k];

                g_x_0_xxzz_z[k] = -g_x_0_xxz_z[k] * ab_z + g_x_0_xxz_zz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyy_x = cbuffer.data(gp_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xyyy_y = cbuffer.data(gp_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xyyy_z = cbuffer.data(gp_geom_10_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyy_x, g_x_0_xyy_xy, g_x_0_xyy_y, g_x_0_xyy_yy, g_x_0_xyy_yz, g_x_0_xyy_z, g_x_0_xyyy_x, g_x_0_xyyy_y, g_x_0_xyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyy_x[k] = -g_x_0_xyy_x[k] * ab_y + g_x_0_xyy_xy[k];

                g_x_0_xyyy_y[k] = -g_x_0_xyy_y[k] * ab_y + g_x_0_xyy_yy[k];

                g_x_0_xyyy_z[k] = -g_x_0_xyy_z[k] * ab_y + g_x_0_xyy_yz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyz_x = cbuffer.data(gp_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xyyz_y = cbuffer.data(gp_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xyyz_z = cbuffer.data(gp_geom_10_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyz_x, g_x_0_xyyz_y, g_x_0_xyyz_z, g_x_0_xyz_x, g_x_0_xyz_xy, g_x_0_xyz_y, g_x_0_xyz_yy, g_x_0_xyz_yz, g_x_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyz_x[k] = -g_x_0_xyz_x[k] * ab_y + g_x_0_xyz_xy[k];

                g_x_0_xyyz_y[k] = -g_x_0_xyz_y[k] * ab_y + g_x_0_xyz_yy[k];

                g_x_0_xyyz_z[k] = -g_x_0_xyz_z[k] * ab_y + g_x_0_xyz_yz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzz_x = cbuffer.data(gp_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xyzz_y = cbuffer.data(gp_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xyzz_z = cbuffer.data(gp_geom_10_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyzz_x, g_x_0_xyzz_y, g_x_0_xyzz_z, g_x_0_xzz_x, g_x_0_xzz_xy, g_x_0_xzz_y, g_x_0_xzz_yy, g_x_0_xzz_yz, g_x_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzz_x[k] = -g_x_0_xzz_x[k] * ab_y + g_x_0_xzz_xy[k];

                g_x_0_xyzz_y[k] = -g_x_0_xzz_y[k] * ab_y + g_x_0_xzz_yy[k];

                g_x_0_xyzz_z[k] = -g_x_0_xzz_z[k] * ab_y + g_x_0_xzz_yz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzz_x = cbuffer.data(gp_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xzzz_y = cbuffer.data(gp_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xzzz_z = cbuffer.data(gp_geom_10_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xzz_x, g_x_0_xzz_xz, g_x_0_xzz_y, g_x_0_xzz_yz, g_x_0_xzz_z, g_x_0_xzz_zz, g_x_0_xzzz_x, g_x_0_xzzz_y, g_x_0_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzz_x[k] = -g_x_0_xzz_x[k] * ab_z + g_x_0_xzz_xz[k];

                g_x_0_xzzz_y[k] = -g_x_0_xzz_y[k] * ab_z + g_x_0_xzz_yz[k];

                g_x_0_xzzz_z[k] = -g_x_0_xzz_z[k] * ab_z + g_x_0_xzz_zz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyy_x = cbuffer.data(gp_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_yyyy_y = cbuffer.data(gp_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_yyyy_z = cbuffer.data(gp_geom_10_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyy_x, g_x_0_yyy_xy, g_x_0_yyy_y, g_x_0_yyy_yy, g_x_0_yyy_yz, g_x_0_yyy_z, g_x_0_yyyy_x, g_x_0_yyyy_y, g_x_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyy_x[k] = -g_x_0_yyy_x[k] * ab_y + g_x_0_yyy_xy[k];

                g_x_0_yyyy_y[k] = -g_x_0_yyy_y[k] * ab_y + g_x_0_yyy_yy[k];

                g_x_0_yyyy_z[k] = -g_x_0_yyy_z[k] * ab_y + g_x_0_yyy_yz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyz_x = cbuffer.data(gp_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_yyyz_y = cbuffer.data(gp_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_yyyz_z = cbuffer.data(gp_geom_10_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyz_x, g_x_0_yyyz_y, g_x_0_yyyz_z, g_x_0_yyz_x, g_x_0_yyz_xy, g_x_0_yyz_y, g_x_0_yyz_yy, g_x_0_yyz_yz, g_x_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyz_x[k] = -g_x_0_yyz_x[k] * ab_y + g_x_0_yyz_xy[k];

                g_x_0_yyyz_y[k] = -g_x_0_yyz_y[k] * ab_y + g_x_0_yyz_yy[k];

                g_x_0_yyyz_z[k] = -g_x_0_yyz_z[k] * ab_y + g_x_0_yyz_yz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzz_x = cbuffer.data(gp_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_yyzz_y = cbuffer.data(gp_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_yyzz_z = cbuffer.data(gp_geom_10_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyzz_x, g_x_0_yyzz_y, g_x_0_yyzz_z, g_x_0_yzz_x, g_x_0_yzz_xy, g_x_0_yzz_y, g_x_0_yzz_yy, g_x_0_yzz_yz, g_x_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzz_x[k] = -g_x_0_yzz_x[k] * ab_y + g_x_0_yzz_xy[k];

                g_x_0_yyzz_y[k] = -g_x_0_yzz_y[k] * ab_y + g_x_0_yzz_yy[k];

                g_x_0_yyzz_z[k] = -g_x_0_yzz_z[k] * ab_y + g_x_0_yzz_yz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzz_x = cbuffer.data(gp_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_yzzz_y = cbuffer.data(gp_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_yzzz_z = cbuffer.data(gp_geom_10_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yzzz_x, g_x_0_yzzz_y, g_x_0_yzzz_z, g_x_0_zzz_x, g_x_0_zzz_xy, g_x_0_zzz_y, g_x_0_zzz_yy, g_x_0_zzz_yz, g_x_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzz_x[k] = -g_x_0_zzz_x[k] * ab_y + g_x_0_zzz_xy[k];

                g_x_0_yzzz_y[k] = -g_x_0_zzz_y[k] * ab_y + g_x_0_zzz_yy[k];

                g_x_0_yzzz_z[k] = -g_x_0_zzz_z[k] * ab_y + g_x_0_zzz_yz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzz_x = cbuffer.data(gp_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_zzzz_y = cbuffer.data(gp_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_zzzz_z = cbuffer.data(gp_geom_10_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzz_x, g_x_0_zzz_xz, g_x_0_zzz_y, g_x_0_zzz_yz, g_x_0_zzz_z, g_x_0_zzz_zz, g_x_0_zzzz_x, g_x_0_zzzz_y, g_x_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzz_x[k] = -g_x_0_zzz_x[k] * ab_z + g_x_0_zzz_xz[k];

                g_x_0_zzzz_y[k] = -g_x_0_zzz_y[k] * ab_z + g_x_0_zzz_yz[k];

                g_x_0_zzzz_z[k] = -g_x_0_zzz_z[k] * ab_z + g_x_0_zzz_zz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxx_x = cbuffer.data(gp_geom_10_off + 45 * ccomps * dcomps);

            auto g_y_0_xxxx_y = cbuffer.data(gp_geom_10_off + 46 * ccomps * dcomps);

            auto g_y_0_xxxx_z = cbuffer.data(gp_geom_10_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxx_x, g_y_0_xxx_xx, g_y_0_xxx_xy, g_y_0_xxx_xz, g_y_0_xxx_y, g_y_0_xxx_z, g_y_0_xxxx_x, g_y_0_xxxx_y, g_y_0_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxx_x[k] = -g_y_0_xxx_x[k] * ab_x + g_y_0_xxx_xx[k];

                g_y_0_xxxx_y[k] = -g_y_0_xxx_y[k] * ab_x + g_y_0_xxx_xy[k];

                g_y_0_xxxx_z[k] = -g_y_0_xxx_z[k] * ab_x + g_y_0_xxx_xz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxy_x = cbuffer.data(gp_geom_10_off + 48 * ccomps * dcomps);

            auto g_y_0_xxxy_y = cbuffer.data(gp_geom_10_off + 49 * ccomps * dcomps);

            auto g_y_0_xxxy_z = cbuffer.data(gp_geom_10_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxy_x, g_y_0_xxxy_y, g_y_0_xxxy_z, g_y_0_xxy_x, g_y_0_xxy_xx, g_y_0_xxy_xy, g_y_0_xxy_xz, g_y_0_xxy_y, g_y_0_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxy_x[k] = -g_y_0_xxy_x[k] * ab_x + g_y_0_xxy_xx[k];

                g_y_0_xxxy_y[k] = -g_y_0_xxy_y[k] * ab_x + g_y_0_xxy_xy[k];

                g_y_0_xxxy_z[k] = -g_y_0_xxy_z[k] * ab_x + g_y_0_xxy_xz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxz_x = cbuffer.data(gp_geom_10_off + 51 * ccomps * dcomps);

            auto g_y_0_xxxz_y = cbuffer.data(gp_geom_10_off + 52 * ccomps * dcomps);

            auto g_y_0_xxxz_z = cbuffer.data(gp_geom_10_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxz_x, g_y_0_xxxz_y, g_y_0_xxxz_z, g_y_0_xxz_x, g_y_0_xxz_xx, g_y_0_xxz_xy, g_y_0_xxz_xz, g_y_0_xxz_y, g_y_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxz_x[k] = -g_y_0_xxz_x[k] * ab_x + g_y_0_xxz_xx[k];

                g_y_0_xxxz_y[k] = -g_y_0_xxz_y[k] * ab_x + g_y_0_xxz_xy[k];

                g_y_0_xxxz_z[k] = -g_y_0_xxz_z[k] * ab_x + g_y_0_xxz_xz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyy_x = cbuffer.data(gp_geom_10_off + 54 * ccomps * dcomps);

            auto g_y_0_xxyy_y = cbuffer.data(gp_geom_10_off + 55 * ccomps * dcomps);

            auto g_y_0_xxyy_z = cbuffer.data(gp_geom_10_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyy_x, g_y_0_xxyy_y, g_y_0_xxyy_z, g_y_0_xyy_x, g_y_0_xyy_xx, g_y_0_xyy_xy, g_y_0_xyy_xz, g_y_0_xyy_y, g_y_0_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyy_x[k] = -g_y_0_xyy_x[k] * ab_x + g_y_0_xyy_xx[k];

                g_y_0_xxyy_y[k] = -g_y_0_xyy_y[k] * ab_x + g_y_0_xyy_xy[k];

                g_y_0_xxyy_z[k] = -g_y_0_xyy_z[k] * ab_x + g_y_0_xyy_xz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyz_x = cbuffer.data(gp_geom_10_off + 57 * ccomps * dcomps);

            auto g_y_0_xxyz_y = cbuffer.data(gp_geom_10_off + 58 * ccomps * dcomps);

            auto g_y_0_xxyz_z = cbuffer.data(gp_geom_10_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyz_x, g_y_0_xxyz_y, g_y_0_xxyz_z, g_y_0_xyz_x, g_y_0_xyz_xx, g_y_0_xyz_xy, g_y_0_xyz_xz, g_y_0_xyz_y, g_y_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyz_x[k] = -g_y_0_xyz_x[k] * ab_x + g_y_0_xyz_xx[k];

                g_y_0_xxyz_y[k] = -g_y_0_xyz_y[k] * ab_x + g_y_0_xyz_xy[k];

                g_y_0_xxyz_z[k] = -g_y_0_xyz_z[k] * ab_x + g_y_0_xyz_xz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzz_x = cbuffer.data(gp_geom_10_off + 60 * ccomps * dcomps);

            auto g_y_0_xxzz_y = cbuffer.data(gp_geom_10_off + 61 * ccomps * dcomps);

            auto g_y_0_xxzz_z = cbuffer.data(gp_geom_10_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxzz_x, g_y_0_xxzz_y, g_y_0_xxzz_z, g_y_0_xzz_x, g_y_0_xzz_xx, g_y_0_xzz_xy, g_y_0_xzz_xz, g_y_0_xzz_y, g_y_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzz_x[k] = -g_y_0_xzz_x[k] * ab_x + g_y_0_xzz_xx[k];

                g_y_0_xxzz_y[k] = -g_y_0_xzz_y[k] * ab_x + g_y_0_xzz_xy[k];

                g_y_0_xxzz_z[k] = -g_y_0_xzz_z[k] * ab_x + g_y_0_xzz_xz[k];
            }

            /// Set up 63-66 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyy_x = cbuffer.data(gp_geom_10_off + 63 * ccomps * dcomps);

            auto g_y_0_xyyy_y = cbuffer.data(gp_geom_10_off + 64 * ccomps * dcomps);

            auto g_y_0_xyyy_z = cbuffer.data(gp_geom_10_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyy_x, g_y_0_xyyy_y, g_y_0_xyyy_z, g_y_0_yyy_x, g_y_0_yyy_xx, g_y_0_yyy_xy, g_y_0_yyy_xz, g_y_0_yyy_y, g_y_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyy_x[k] = -g_y_0_yyy_x[k] * ab_x + g_y_0_yyy_xx[k];

                g_y_0_xyyy_y[k] = -g_y_0_yyy_y[k] * ab_x + g_y_0_yyy_xy[k];

                g_y_0_xyyy_z[k] = -g_y_0_yyy_z[k] * ab_x + g_y_0_yyy_xz[k];
            }

            /// Set up 66-69 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyz_x = cbuffer.data(gp_geom_10_off + 66 * ccomps * dcomps);

            auto g_y_0_xyyz_y = cbuffer.data(gp_geom_10_off + 67 * ccomps * dcomps);

            auto g_y_0_xyyz_z = cbuffer.data(gp_geom_10_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyz_x, g_y_0_xyyz_y, g_y_0_xyyz_z, g_y_0_yyz_x, g_y_0_yyz_xx, g_y_0_yyz_xy, g_y_0_yyz_xz, g_y_0_yyz_y, g_y_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyz_x[k] = -g_y_0_yyz_x[k] * ab_x + g_y_0_yyz_xx[k];

                g_y_0_xyyz_y[k] = -g_y_0_yyz_y[k] * ab_x + g_y_0_yyz_xy[k];

                g_y_0_xyyz_z[k] = -g_y_0_yyz_z[k] * ab_x + g_y_0_yyz_xz[k];
            }

            /// Set up 69-72 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzz_x = cbuffer.data(gp_geom_10_off + 69 * ccomps * dcomps);

            auto g_y_0_xyzz_y = cbuffer.data(gp_geom_10_off + 70 * ccomps * dcomps);

            auto g_y_0_xyzz_z = cbuffer.data(gp_geom_10_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyzz_x, g_y_0_xyzz_y, g_y_0_xyzz_z, g_y_0_yzz_x, g_y_0_yzz_xx, g_y_0_yzz_xy, g_y_0_yzz_xz, g_y_0_yzz_y, g_y_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzz_x[k] = -g_y_0_yzz_x[k] * ab_x + g_y_0_yzz_xx[k];

                g_y_0_xyzz_y[k] = -g_y_0_yzz_y[k] * ab_x + g_y_0_yzz_xy[k];

                g_y_0_xyzz_z[k] = -g_y_0_yzz_z[k] * ab_x + g_y_0_yzz_xz[k];
            }

            /// Set up 72-75 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzz_x = cbuffer.data(gp_geom_10_off + 72 * ccomps * dcomps);

            auto g_y_0_xzzz_y = cbuffer.data(gp_geom_10_off + 73 * ccomps * dcomps);

            auto g_y_0_xzzz_z = cbuffer.data(gp_geom_10_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xzzz_x, g_y_0_xzzz_y, g_y_0_xzzz_z, g_y_0_zzz_x, g_y_0_zzz_xx, g_y_0_zzz_xy, g_y_0_zzz_xz, g_y_0_zzz_y, g_y_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzz_x[k] = -g_y_0_zzz_x[k] * ab_x + g_y_0_zzz_xx[k];

                g_y_0_xzzz_y[k] = -g_y_0_zzz_y[k] * ab_x + g_y_0_zzz_xy[k];

                g_y_0_xzzz_z[k] = -g_y_0_zzz_z[k] * ab_x + g_y_0_zzz_xz[k];
            }

            /// Set up 75-78 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyy_x = cbuffer.data(gp_geom_10_off + 75 * ccomps * dcomps);

            auto g_y_0_yyyy_y = cbuffer.data(gp_geom_10_off + 76 * ccomps * dcomps);

            auto g_y_0_yyyy_z = cbuffer.data(gp_geom_10_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyy_x, g_y_0_yyy_xy, g_y_0_yyy_y, g_y_0_yyy_yy, g_y_0_yyy_yz, g_y_0_yyy_z, g_y_0_yyyy_x, g_y_0_yyyy_y, g_y_0_yyyy_z, g_yyy_x, g_yyy_y, g_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyy_x[k] = -g_yyy_x[k] - g_y_0_yyy_x[k] * ab_y + g_y_0_yyy_xy[k];

                g_y_0_yyyy_y[k] = -g_yyy_y[k] - g_y_0_yyy_y[k] * ab_y + g_y_0_yyy_yy[k];

                g_y_0_yyyy_z[k] = -g_yyy_z[k] - g_y_0_yyy_z[k] * ab_y + g_y_0_yyy_yz[k];
            }

            /// Set up 78-81 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyz_x = cbuffer.data(gp_geom_10_off + 78 * ccomps * dcomps);

            auto g_y_0_yyyz_y = cbuffer.data(gp_geom_10_off + 79 * ccomps * dcomps);

            auto g_y_0_yyyz_z = cbuffer.data(gp_geom_10_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyy_x, g_y_0_yyy_xz, g_y_0_yyy_y, g_y_0_yyy_yz, g_y_0_yyy_z, g_y_0_yyy_zz, g_y_0_yyyz_x, g_y_0_yyyz_y, g_y_0_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyz_x[k] = -g_y_0_yyy_x[k] * ab_z + g_y_0_yyy_xz[k];

                g_y_0_yyyz_y[k] = -g_y_0_yyy_y[k] * ab_z + g_y_0_yyy_yz[k];

                g_y_0_yyyz_z[k] = -g_y_0_yyy_z[k] * ab_z + g_y_0_yyy_zz[k];
            }

            /// Set up 81-84 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzz_x = cbuffer.data(gp_geom_10_off + 81 * ccomps * dcomps);

            auto g_y_0_yyzz_y = cbuffer.data(gp_geom_10_off + 82 * ccomps * dcomps);

            auto g_y_0_yyzz_z = cbuffer.data(gp_geom_10_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyz_x, g_y_0_yyz_xz, g_y_0_yyz_y, g_y_0_yyz_yz, g_y_0_yyz_z, g_y_0_yyz_zz, g_y_0_yyzz_x, g_y_0_yyzz_y, g_y_0_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzz_x[k] = -g_y_0_yyz_x[k] * ab_z + g_y_0_yyz_xz[k];

                g_y_0_yyzz_y[k] = -g_y_0_yyz_y[k] * ab_z + g_y_0_yyz_yz[k];

                g_y_0_yyzz_z[k] = -g_y_0_yyz_z[k] * ab_z + g_y_0_yyz_zz[k];
            }

            /// Set up 84-87 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzz_x = cbuffer.data(gp_geom_10_off + 84 * ccomps * dcomps);

            auto g_y_0_yzzz_y = cbuffer.data(gp_geom_10_off + 85 * ccomps * dcomps);

            auto g_y_0_yzzz_z = cbuffer.data(gp_geom_10_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yzz_x, g_y_0_yzz_xz, g_y_0_yzz_y, g_y_0_yzz_yz, g_y_0_yzz_z, g_y_0_yzz_zz, g_y_0_yzzz_x, g_y_0_yzzz_y, g_y_0_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzz_x[k] = -g_y_0_yzz_x[k] * ab_z + g_y_0_yzz_xz[k];

                g_y_0_yzzz_y[k] = -g_y_0_yzz_y[k] * ab_z + g_y_0_yzz_yz[k];

                g_y_0_yzzz_z[k] = -g_y_0_yzz_z[k] * ab_z + g_y_0_yzz_zz[k];
            }

            /// Set up 87-90 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzz_x = cbuffer.data(gp_geom_10_off + 87 * ccomps * dcomps);

            auto g_y_0_zzzz_y = cbuffer.data(gp_geom_10_off + 88 * ccomps * dcomps);

            auto g_y_0_zzzz_z = cbuffer.data(gp_geom_10_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzz_x, g_y_0_zzz_xz, g_y_0_zzz_y, g_y_0_zzz_yz, g_y_0_zzz_z, g_y_0_zzz_zz, g_y_0_zzzz_x, g_y_0_zzzz_y, g_y_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzz_x[k] = -g_y_0_zzz_x[k] * ab_z + g_y_0_zzz_xz[k];

                g_y_0_zzzz_y[k] = -g_y_0_zzz_y[k] * ab_z + g_y_0_zzz_yz[k];

                g_y_0_zzzz_z[k] = -g_y_0_zzz_z[k] * ab_z + g_y_0_zzz_zz[k];
            }

            /// Set up 90-93 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxx_x = cbuffer.data(gp_geom_10_off + 90 * ccomps * dcomps);

            auto g_z_0_xxxx_y = cbuffer.data(gp_geom_10_off + 91 * ccomps * dcomps);

            auto g_z_0_xxxx_z = cbuffer.data(gp_geom_10_off + 92 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxx_x, g_z_0_xxx_xx, g_z_0_xxx_xy, g_z_0_xxx_xz, g_z_0_xxx_y, g_z_0_xxx_z, g_z_0_xxxx_x, g_z_0_xxxx_y, g_z_0_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxx_x[k] = -g_z_0_xxx_x[k] * ab_x + g_z_0_xxx_xx[k];

                g_z_0_xxxx_y[k] = -g_z_0_xxx_y[k] * ab_x + g_z_0_xxx_xy[k];

                g_z_0_xxxx_z[k] = -g_z_0_xxx_z[k] * ab_x + g_z_0_xxx_xz[k];
            }

            /// Set up 93-96 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxy_x = cbuffer.data(gp_geom_10_off + 93 * ccomps * dcomps);

            auto g_z_0_xxxy_y = cbuffer.data(gp_geom_10_off + 94 * ccomps * dcomps);

            auto g_z_0_xxxy_z = cbuffer.data(gp_geom_10_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxy_x, g_z_0_xxxy_y, g_z_0_xxxy_z, g_z_0_xxy_x, g_z_0_xxy_xx, g_z_0_xxy_xy, g_z_0_xxy_xz, g_z_0_xxy_y, g_z_0_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxy_x[k] = -g_z_0_xxy_x[k] * ab_x + g_z_0_xxy_xx[k];

                g_z_0_xxxy_y[k] = -g_z_0_xxy_y[k] * ab_x + g_z_0_xxy_xy[k];

                g_z_0_xxxy_z[k] = -g_z_0_xxy_z[k] * ab_x + g_z_0_xxy_xz[k];
            }

            /// Set up 96-99 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxz_x = cbuffer.data(gp_geom_10_off + 96 * ccomps * dcomps);

            auto g_z_0_xxxz_y = cbuffer.data(gp_geom_10_off + 97 * ccomps * dcomps);

            auto g_z_0_xxxz_z = cbuffer.data(gp_geom_10_off + 98 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxz_x, g_z_0_xxxz_y, g_z_0_xxxz_z, g_z_0_xxz_x, g_z_0_xxz_xx, g_z_0_xxz_xy, g_z_0_xxz_xz, g_z_0_xxz_y, g_z_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxz_x[k] = -g_z_0_xxz_x[k] * ab_x + g_z_0_xxz_xx[k];

                g_z_0_xxxz_y[k] = -g_z_0_xxz_y[k] * ab_x + g_z_0_xxz_xy[k];

                g_z_0_xxxz_z[k] = -g_z_0_xxz_z[k] * ab_x + g_z_0_xxz_xz[k];
            }

            /// Set up 99-102 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyy_x = cbuffer.data(gp_geom_10_off + 99 * ccomps * dcomps);

            auto g_z_0_xxyy_y = cbuffer.data(gp_geom_10_off + 100 * ccomps * dcomps);

            auto g_z_0_xxyy_z = cbuffer.data(gp_geom_10_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyy_x, g_z_0_xxyy_y, g_z_0_xxyy_z, g_z_0_xyy_x, g_z_0_xyy_xx, g_z_0_xyy_xy, g_z_0_xyy_xz, g_z_0_xyy_y, g_z_0_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyy_x[k] = -g_z_0_xyy_x[k] * ab_x + g_z_0_xyy_xx[k];

                g_z_0_xxyy_y[k] = -g_z_0_xyy_y[k] * ab_x + g_z_0_xyy_xy[k];

                g_z_0_xxyy_z[k] = -g_z_0_xyy_z[k] * ab_x + g_z_0_xyy_xz[k];
            }

            /// Set up 102-105 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyz_x = cbuffer.data(gp_geom_10_off + 102 * ccomps * dcomps);

            auto g_z_0_xxyz_y = cbuffer.data(gp_geom_10_off + 103 * ccomps * dcomps);

            auto g_z_0_xxyz_z = cbuffer.data(gp_geom_10_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyz_x, g_z_0_xxyz_y, g_z_0_xxyz_z, g_z_0_xyz_x, g_z_0_xyz_xx, g_z_0_xyz_xy, g_z_0_xyz_xz, g_z_0_xyz_y, g_z_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyz_x[k] = -g_z_0_xyz_x[k] * ab_x + g_z_0_xyz_xx[k];

                g_z_0_xxyz_y[k] = -g_z_0_xyz_y[k] * ab_x + g_z_0_xyz_xy[k];

                g_z_0_xxyz_z[k] = -g_z_0_xyz_z[k] * ab_x + g_z_0_xyz_xz[k];
            }

            /// Set up 105-108 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzz_x = cbuffer.data(gp_geom_10_off + 105 * ccomps * dcomps);

            auto g_z_0_xxzz_y = cbuffer.data(gp_geom_10_off + 106 * ccomps * dcomps);

            auto g_z_0_xxzz_z = cbuffer.data(gp_geom_10_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxzz_x, g_z_0_xxzz_y, g_z_0_xxzz_z, g_z_0_xzz_x, g_z_0_xzz_xx, g_z_0_xzz_xy, g_z_0_xzz_xz, g_z_0_xzz_y, g_z_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzz_x[k] = -g_z_0_xzz_x[k] * ab_x + g_z_0_xzz_xx[k];

                g_z_0_xxzz_y[k] = -g_z_0_xzz_y[k] * ab_x + g_z_0_xzz_xy[k];

                g_z_0_xxzz_z[k] = -g_z_0_xzz_z[k] * ab_x + g_z_0_xzz_xz[k];
            }

            /// Set up 108-111 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyy_x = cbuffer.data(gp_geom_10_off + 108 * ccomps * dcomps);

            auto g_z_0_xyyy_y = cbuffer.data(gp_geom_10_off + 109 * ccomps * dcomps);

            auto g_z_0_xyyy_z = cbuffer.data(gp_geom_10_off + 110 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyy_x, g_z_0_xyyy_y, g_z_0_xyyy_z, g_z_0_yyy_x, g_z_0_yyy_xx, g_z_0_yyy_xy, g_z_0_yyy_xz, g_z_0_yyy_y, g_z_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyy_x[k] = -g_z_0_yyy_x[k] * ab_x + g_z_0_yyy_xx[k];

                g_z_0_xyyy_y[k] = -g_z_0_yyy_y[k] * ab_x + g_z_0_yyy_xy[k];

                g_z_0_xyyy_z[k] = -g_z_0_yyy_z[k] * ab_x + g_z_0_yyy_xz[k];
            }

            /// Set up 111-114 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyz_x = cbuffer.data(gp_geom_10_off + 111 * ccomps * dcomps);

            auto g_z_0_xyyz_y = cbuffer.data(gp_geom_10_off + 112 * ccomps * dcomps);

            auto g_z_0_xyyz_z = cbuffer.data(gp_geom_10_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyz_x, g_z_0_xyyz_y, g_z_0_xyyz_z, g_z_0_yyz_x, g_z_0_yyz_xx, g_z_0_yyz_xy, g_z_0_yyz_xz, g_z_0_yyz_y, g_z_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyz_x[k] = -g_z_0_yyz_x[k] * ab_x + g_z_0_yyz_xx[k];

                g_z_0_xyyz_y[k] = -g_z_0_yyz_y[k] * ab_x + g_z_0_yyz_xy[k];

                g_z_0_xyyz_z[k] = -g_z_0_yyz_z[k] * ab_x + g_z_0_yyz_xz[k];
            }

            /// Set up 114-117 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzz_x = cbuffer.data(gp_geom_10_off + 114 * ccomps * dcomps);

            auto g_z_0_xyzz_y = cbuffer.data(gp_geom_10_off + 115 * ccomps * dcomps);

            auto g_z_0_xyzz_z = cbuffer.data(gp_geom_10_off + 116 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyzz_x, g_z_0_xyzz_y, g_z_0_xyzz_z, g_z_0_yzz_x, g_z_0_yzz_xx, g_z_0_yzz_xy, g_z_0_yzz_xz, g_z_0_yzz_y, g_z_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzz_x[k] = -g_z_0_yzz_x[k] * ab_x + g_z_0_yzz_xx[k];

                g_z_0_xyzz_y[k] = -g_z_0_yzz_y[k] * ab_x + g_z_0_yzz_xy[k];

                g_z_0_xyzz_z[k] = -g_z_0_yzz_z[k] * ab_x + g_z_0_yzz_xz[k];
            }

            /// Set up 117-120 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzz_x = cbuffer.data(gp_geom_10_off + 117 * ccomps * dcomps);

            auto g_z_0_xzzz_y = cbuffer.data(gp_geom_10_off + 118 * ccomps * dcomps);

            auto g_z_0_xzzz_z = cbuffer.data(gp_geom_10_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xzzz_x, g_z_0_xzzz_y, g_z_0_xzzz_z, g_z_0_zzz_x, g_z_0_zzz_xx, g_z_0_zzz_xy, g_z_0_zzz_xz, g_z_0_zzz_y, g_z_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzz_x[k] = -g_z_0_zzz_x[k] * ab_x + g_z_0_zzz_xx[k];

                g_z_0_xzzz_y[k] = -g_z_0_zzz_y[k] * ab_x + g_z_0_zzz_xy[k];

                g_z_0_xzzz_z[k] = -g_z_0_zzz_z[k] * ab_x + g_z_0_zzz_xz[k];
            }

            /// Set up 120-123 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyy_x = cbuffer.data(gp_geom_10_off + 120 * ccomps * dcomps);

            auto g_z_0_yyyy_y = cbuffer.data(gp_geom_10_off + 121 * ccomps * dcomps);

            auto g_z_0_yyyy_z = cbuffer.data(gp_geom_10_off + 122 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyy_x, g_z_0_yyy_xy, g_z_0_yyy_y, g_z_0_yyy_yy, g_z_0_yyy_yz, g_z_0_yyy_z, g_z_0_yyyy_x, g_z_0_yyyy_y, g_z_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyy_x[k] = -g_z_0_yyy_x[k] * ab_y + g_z_0_yyy_xy[k];

                g_z_0_yyyy_y[k] = -g_z_0_yyy_y[k] * ab_y + g_z_0_yyy_yy[k];

                g_z_0_yyyy_z[k] = -g_z_0_yyy_z[k] * ab_y + g_z_0_yyy_yz[k];
            }

            /// Set up 123-126 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyz_x = cbuffer.data(gp_geom_10_off + 123 * ccomps * dcomps);

            auto g_z_0_yyyz_y = cbuffer.data(gp_geom_10_off + 124 * ccomps * dcomps);

            auto g_z_0_yyyz_z = cbuffer.data(gp_geom_10_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyz_x, g_z_0_yyyz_y, g_z_0_yyyz_z, g_z_0_yyz_x, g_z_0_yyz_xy, g_z_0_yyz_y, g_z_0_yyz_yy, g_z_0_yyz_yz, g_z_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyz_x[k] = -g_z_0_yyz_x[k] * ab_y + g_z_0_yyz_xy[k];

                g_z_0_yyyz_y[k] = -g_z_0_yyz_y[k] * ab_y + g_z_0_yyz_yy[k];

                g_z_0_yyyz_z[k] = -g_z_0_yyz_z[k] * ab_y + g_z_0_yyz_yz[k];
            }

            /// Set up 126-129 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzz_x = cbuffer.data(gp_geom_10_off + 126 * ccomps * dcomps);

            auto g_z_0_yyzz_y = cbuffer.data(gp_geom_10_off + 127 * ccomps * dcomps);

            auto g_z_0_yyzz_z = cbuffer.data(gp_geom_10_off + 128 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyzz_x, g_z_0_yyzz_y, g_z_0_yyzz_z, g_z_0_yzz_x, g_z_0_yzz_xy, g_z_0_yzz_y, g_z_0_yzz_yy, g_z_0_yzz_yz, g_z_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzz_x[k] = -g_z_0_yzz_x[k] * ab_y + g_z_0_yzz_xy[k];

                g_z_0_yyzz_y[k] = -g_z_0_yzz_y[k] * ab_y + g_z_0_yzz_yy[k];

                g_z_0_yyzz_z[k] = -g_z_0_yzz_z[k] * ab_y + g_z_0_yzz_yz[k];
            }

            /// Set up 129-132 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzz_x = cbuffer.data(gp_geom_10_off + 129 * ccomps * dcomps);

            auto g_z_0_yzzz_y = cbuffer.data(gp_geom_10_off + 130 * ccomps * dcomps);

            auto g_z_0_yzzz_z = cbuffer.data(gp_geom_10_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yzzz_x, g_z_0_yzzz_y, g_z_0_yzzz_z, g_z_0_zzz_x, g_z_0_zzz_xy, g_z_0_zzz_y, g_z_0_zzz_yy, g_z_0_zzz_yz, g_z_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzz_x[k] = -g_z_0_zzz_x[k] * ab_y + g_z_0_zzz_xy[k];

                g_z_0_yzzz_y[k] = -g_z_0_zzz_y[k] * ab_y + g_z_0_zzz_yy[k];

                g_z_0_yzzz_z[k] = -g_z_0_zzz_z[k] * ab_y + g_z_0_zzz_yz[k];
            }

            /// Set up 132-135 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzz_x = cbuffer.data(gp_geom_10_off + 132 * ccomps * dcomps);

            auto g_z_0_zzzz_y = cbuffer.data(gp_geom_10_off + 133 * ccomps * dcomps);

            auto g_z_0_zzzz_z = cbuffer.data(gp_geom_10_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzz_x, g_z_0_zzz_xz, g_z_0_zzz_y, g_z_0_zzz_yz, g_z_0_zzz_z, g_z_0_zzz_zz, g_z_0_zzzz_x, g_z_0_zzzz_y, g_z_0_zzzz_z, g_zzz_x, g_zzz_y, g_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzz_x[k] = -g_zzz_x[k] - g_z_0_zzz_x[k] * ab_z + g_z_0_zzz_xz[k];

                g_z_0_zzzz_y[k] = -g_zzz_y[k] - g_z_0_zzz_y[k] * ab_z + g_z_0_zzz_yz[k];

                g_z_0_zzzz_z[k] = -g_zzz_z[k] - g_z_0_zzz_z[k] * ab_z + g_z_0_zzz_zz[k];
            }
        }
    }
}

} // erirec namespace

