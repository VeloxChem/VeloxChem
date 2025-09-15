#include "ElectronRepulsionGeom0100ContrRecGPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_gpxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_gpxx,
                                            const size_t idx_fpxx,
                                            const size_t idx_geom_01_fpxx,
                                            const size_t idx_geom_01_fdxx,
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

            auto g_xxy_x = cbuffer.data(fp_off + 3 * ccomps * dcomps);

            auto g_xxy_y = cbuffer.data(fp_off + 4 * ccomps * dcomps);

            auto g_xxy_z = cbuffer.data(fp_off + 5 * ccomps * dcomps);

            auto g_xxz_x = cbuffer.data(fp_off + 6 * ccomps * dcomps);

            auto g_xxz_y = cbuffer.data(fp_off + 7 * ccomps * dcomps);

            auto g_xxz_z = cbuffer.data(fp_off + 8 * ccomps * dcomps);

            auto g_xyy_x = cbuffer.data(fp_off + 9 * ccomps * dcomps);

            auto g_xyy_y = cbuffer.data(fp_off + 10 * ccomps * dcomps);

            auto g_xyy_z = cbuffer.data(fp_off + 11 * ccomps * dcomps);

            auto g_xyz_x = cbuffer.data(fp_off + 12 * ccomps * dcomps);

            auto g_xyz_y = cbuffer.data(fp_off + 13 * ccomps * dcomps);

            auto g_xyz_z = cbuffer.data(fp_off + 14 * ccomps * dcomps);

            auto g_xzz_x = cbuffer.data(fp_off + 15 * ccomps * dcomps);

            auto g_xzz_y = cbuffer.data(fp_off + 16 * ccomps * dcomps);

            auto g_xzz_z = cbuffer.data(fp_off + 17 * ccomps * dcomps);

            auto g_yyy_x = cbuffer.data(fp_off + 18 * ccomps * dcomps);

            auto g_yyy_y = cbuffer.data(fp_off + 19 * ccomps * dcomps);

            auto g_yyy_z = cbuffer.data(fp_off + 20 * ccomps * dcomps);

            auto g_yyz_x = cbuffer.data(fp_off + 21 * ccomps * dcomps);

            auto g_yyz_y = cbuffer.data(fp_off + 22 * ccomps * dcomps);

            auto g_yyz_z = cbuffer.data(fp_off + 23 * ccomps * dcomps);

            auto g_yzz_x = cbuffer.data(fp_off + 24 * ccomps * dcomps);

            auto g_yzz_y = cbuffer.data(fp_off + 25 * ccomps * dcomps);

            auto g_yzz_z = cbuffer.data(fp_off + 26 * ccomps * dcomps);

            auto g_zzz_x = cbuffer.data(fp_off + 27 * ccomps * dcomps);

            auto g_zzz_y = cbuffer.data(fp_off + 28 * ccomps * dcomps);

            auto g_zzz_z = cbuffer.data(fp_off + 29 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FPSS

            const auto fp_geom_01_off = idx_geom_01_fpxx + i * dcomps + j;

            auto g_0_x_xxx_x = cbuffer.data(fp_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxx_y = cbuffer.data(fp_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxx_z = cbuffer.data(fp_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxy_x = cbuffer.data(fp_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxy_y = cbuffer.data(fp_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxy_z = cbuffer.data(fp_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxz_x = cbuffer.data(fp_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxz_y = cbuffer.data(fp_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxz_z = cbuffer.data(fp_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xyy_x = cbuffer.data(fp_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xyy_y = cbuffer.data(fp_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xyy_z = cbuffer.data(fp_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xyz_x = cbuffer.data(fp_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xyz_y = cbuffer.data(fp_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xyz_z = cbuffer.data(fp_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xzz_x = cbuffer.data(fp_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xzz_y = cbuffer.data(fp_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xzz_z = cbuffer.data(fp_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_yyy_x = cbuffer.data(fp_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_yyy_y = cbuffer.data(fp_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_yyy_z = cbuffer.data(fp_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_yyz_x = cbuffer.data(fp_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_yyz_y = cbuffer.data(fp_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_yyz_z = cbuffer.data(fp_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_yzz_x = cbuffer.data(fp_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_yzz_y = cbuffer.data(fp_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_yzz_z = cbuffer.data(fp_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_zzz_x = cbuffer.data(fp_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_zzz_y = cbuffer.data(fp_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_zzz_z = cbuffer.data(fp_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_y_xxx_x = cbuffer.data(fp_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_y_xxx_y = cbuffer.data(fp_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_y_xxx_z = cbuffer.data(fp_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_y_xxy_x = cbuffer.data(fp_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_y_xxy_y = cbuffer.data(fp_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_y_xxy_z = cbuffer.data(fp_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_y_xxz_x = cbuffer.data(fp_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_y_xxz_y = cbuffer.data(fp_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_y_xxz_z = cbuffer.data(fp_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_y_xyy_x = cbuffer.data(fp_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_y_xyy_y = cbuffer.data(fp_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_y_xyy_z = cbuffer.data(fp_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_y_xyz_x = cbuffer.data(fp_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_y_xyz_y = cbuffer.data(fp_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_y_xyz_z = cbuffer.data(fp_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_y_xzz_x = cbuffer.data(fp_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_y_xzz_y = cbuffer.data(fp_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_y_xzz_z = cbuffer.data(fp_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_y_yyy_x = cbuffer.data(fp_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_y_yyy_y = cbuffer.data(fp_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_y_yyy_z = cbuffer.data(fp_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_y_yyz_x = cbuffer.data(fp_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_y_yyz_y = cbuffer.data(fp_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_y_yyz_z = cbuffer.data(fp_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_y_yzz_x = cbuffer.data(fp_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_y_yzz_y = cbuffer.data(fp_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_y_yzz_z = cbuffer.data(fp_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_y_zzz_x = cbuffer.data(fp_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_y_zzz_y = cbuffer.data(fp_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_y_zzz_z = cbuffer.data(fp_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_z_xxx_x = cbuffer.data(fp_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_z_xxx_y = cbuffer.data(fp_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_z_xxx_z = cbuffer.data(fp_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_z_xxy_x = cbuffer.data(fp_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_z_xxy_y = cbuffer.data(fp_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_z_xxy_z = cbuffer.data(fp_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_z_xxz_x = cbuffer.data(fp_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_z_xxz_y = cbuffer.data(fp_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_z_xxz_z = cbuffer.data(fp_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_z_xyy_x = cbuffer.data(fp_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_z_xyy_y = cbuffer.data(fp_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_z_xyy_z = cbuffer.data(fp_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_z_xyz_x = cbuffer.data(fp_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_z_xyz_y = cbuffer.data(fp_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_z_xyz_z = cbuffer.data(fp_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_z_xzz_x = cbuffer.data(fp_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_z_xzz_y = cbuffer.data(fp_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_z_xzz_z = cbuffer.data(fp_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_z_yyy_x = cbuffer.data(fp_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_z_yyy_y = cbuffer.data(fp_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_z_yyy_z = cbuffer.data(fp_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_z_yyz_x = cbuffer.data(fp_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_z_yyz_y = cbuffer.data(fp_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_z_yyz_z = cbuffer.data(fp_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_z_yzz_x = cbuffer.data(fp_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_z_yzz_y = cbuffer.data(fp_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_z_yzz_z = cbuffer.data(fp_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_z_zzz_x = cbuffer.data(fp_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_z_zzz_y = cbuffer.data(fp_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_z_zzz_z = cbuffer.data(fp_geom_01_off + 89 * ccomps * dcomps);

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

            /// set up bra offset for contr_buffer_gpxx

            const auto gp_geom_01_off = idx_geom_01_gpxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxx_x = cbuffer.data(gp_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxx_y = cbuffer.data(gp_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxx_z = cbuffer.data(gp_geom_01_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxx_x, g_0_x_xxx_xx, g_0_x_xxx_xy, g_0_x_xxx_xz, g_0_x_xxx_y, g_0_x_xxx_z, g_0_x_xxxx_x, g_0_x_xxxx_y, g_0_x_xxxx_z, g_xxx_x, g_xxx_y, g_xxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxx_x[k] = g_xxx_x[k] - g_0_x_xxx_x[k] * ab_x + g_0_x_xxx_xx[k];

                g_0_x_xxxx_y[k] = g_xxx_y[k] - g_0_x_xxx_y[k] * ab_x + g_0_x_xxx_xy[k];

                g_0_x_xxxx_z[k] = g_xxx_z[k] - g_0_x_xxx_z[k] * ab_x + g_0_x_xxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxy_x = cbuffer.data(gp_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxy_y = cbuffer.data(gp_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxy_z = cbuffer.data(gp_geom_01_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxx_x, g_0_x_xxx_xy, g_0_x_xxx_y, g_0_x_xxx_yy, g_0_x_xxx_yz, g_0_x_xxx_z, g_0_x_xxxy_x, g_0_x_xxxy_y, g_0_x_xxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxy_x[k] = -g_0_x_xxx_x[k] * ab_y + g_0_x_xxx_xy[k];

                g_0_x_xxxy_y[k] = -g_0_x_xxx_y[k] * ab_y + g_0_x_xxx_yy[k];

                g_0_x_xxxy_z[k] = -g_0_x_xxx_z[k] * ab_y + g_0_x_xxx_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxz_x = cbuffer.data(gp_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxz_y = cbuffer.data(gp_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxz_z = cbuffer.data(gp_geom_01_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxx_x, g_0_x_xxx_xz, g_0_x_xxx_y, g_0_x_xxx_yz, g_0_x_xxx_z, g_0_x_xxx_zz, g_0_x_xxxz_x, g_0_x_xxxz_y, g_0_x_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxz_x[k] = -g_0_x_xxx_x[k] * ab_z + g_0_x_xxx_xz[k];

                g_0_x_xxxz_y[k] = -g_0_x_xxx_y[k] * ab_z + g_0_x_xxx_yz[k];

                g_0_x_xxxz_z[k] = -g_0_x_xxx_z[k] * ab_z + g_0_x_xxx_zz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyy_x = cbuffer.data(gp_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxyy_y = cbuffer.data(gp_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxyy_z = cbuffer.data(gp_geom_01_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxy_x, g_0_x_xxy_xy, g_0_x_xxy_y, g_0_x_xxy_yy, g_0_x_xxy_yz, g_0_x_xxy_z, g_0_x_xxyy_x, g_0_x_xxyy_y, g_0_x_xxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyy_x[k] = -g_0_x_xxy_x[k] * ab_y + g_0_x_xxy_xy[k];

                g_0_x_xxyy_y[k] = -g_0_x_xxy_y[k] * ab_y + g_0_x_xxy_yy[k];

                g_0_x_xxyy_z[k] = -g_0_x_xxy_z[k] * ab_y + g_0_x_xxy_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyz_x = cbuffer.data(gp_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxyz_y = cbuffer.data(gp_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxyz_z = cbuffer.data(gp_geom_01_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyz_x, g_0_x_xxyz_y, g_0_x_xxyz_z, g_0_x_xxz_x, g_0_x_xxz_xy, g_0_x_xxz_y, g_0_x_xxz_yy, g_0_x_xxz_yz, g_0_x_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyz_x[k] = -g_0_x_xxz_x[k] * ab_y + g_0_x_xxz_xy[k];

                g_0_x_xxyz_y[k] = -g_0_x_xxz_y[k] * ab_y + g_0_x_xxz_yy[k];

                g_0_x_xxyz_z[k] = -g_0_x_xxz_z[k] * ab_y + g_0_x_xxz_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxzz_x = cbuffer.data(gp_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxzz_y = cbuffer.data(gp_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxzz_z = cbuffer.data(gp_geom_01_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxz_x, g_0_x_xxz_xz, g_0_x_xxz_y, g_0_x_xxz_yz, g_0_x_xxz_z, g_0_x_xxz_zz, g_0_x_xxzz_x, g_0_x_xxzz_y, g_0_x_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxzz_x[k] = -g_0_x_xxz_x[k] * ab_z + g_0_x_xxz_xz[k];

                g_0_x_xxzz_y[k] = -g_0_x_xxz_y[k] * ab_z + g_0_x_xxz_yz[k];

                g_0_x_xxzz_z[k] = -g_0_x_xxz_z[k] * ab_z + g_0_x_xxz_zz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyy_x = cbuffer.data(gp_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xyyy_y = cbuffer.data(gp_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xyyy_z = cbuffer.data(gp_geom_01_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyy_x, g_0_x_xyy_xy, g_0_x_xyy_y, g_0_x_xyy_yy, g_0_x_xyy_yz, g_0_x_xyy_z, g_0_x_xyyy_x, g_0_x_xyyy_y, g_0_x_xyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyy_x[k] = -g_0_x_xyy_x[k] * ab_y + g_0_x_xyy_xy[k];

                g_0_x_xyyy_y[k] = -g_0_x_xyy_y[k] * ab_y + g_0_x_xyy_yy[k];

                g_0_x_xyyy_z[k] = -g_0_x_xyy_z[k] * ab_y + g_0_x_xyy_yz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyz_x = cbuffer.data(gp_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xyyz_y = cbuffer.data(gp_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xyyz_z = cbuffer.data(gp_geom_01_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyz_x, g_0_x_xyyz_y, g_0_x_xyyz_z, g_0_x_xyz_x, g_0_x_xyz_xy, g_0_x_xyz_y, g_0_x_xyz_yy, g_0_x_xyz_yz, g_0_x_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyz_x[k] = -g_0_x_xyz_x[k] * ab_y + g_0_x_xyz_xy[k];

                g_0_x_xyyz_y[k] = -g_0_x_xyz_y[k] * ab_y + g_0_x_xyz_yy[k];

                g_0_x_xyyz_z[k] = -g_0_x_xyz_z[k] * ab_y + g_0_x_xyz_yz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyzz_x = cbuffer.data(gp_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xyzz_y = cbuffer.data(gp_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xyzz_z = cbuffer.data(gp_geom_01_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyzz_x, g_0_x_xyzz_y, g_0_x_xyzz_z, g_0_x_xzz_x, g_0_x_xzz_xy, g_0_x_xzz_y, g_0_x_xzz_yy, g_0_x_xzz_yz, g_0_x_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyzz_x[k] = -g_0_x_xzz_x[k] * ab_y + g_0_x_xzz_xy[k];

                g_0_x_xyzz_y[k] = -g_0_x_xzz_y[k] * ab_y + g_0_x_xzz_yy[k];

                g_0_x_xyzz_z[k] = -g_0_x_xzz_z[k] * ab_y + g_0_x_xzz_yz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_0_x_xzzz_x = cbuffer.data(gp_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xzzz_y = cbuffer.data(gp_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xzzz_z = cbuffer.data(gp_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xzz_x, g_0_x_xzz_xz, g_0_x_xzz_y, g_0_x_xzz_yz, g_0_x_xzz_z, g_0_x_xzz_zz, g_0_x_xzzz_x, g_0_x_xzzz_y, g_0_x_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzzz_x[k] = -g_0_x_xzz_x[k] * ab_z + g_0_x_xzz_xz[k];

                g_0_x_xzzz_y[k] = -g_0_x_xzz_y[k] * ab_z + g_0_x_xzz_yz[k];

                g_0_x_xzzz_z[k] = -g_0_x_xzz_z[k] * ab_z + g_0_x_xzz_zz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyy_x = cbuffer.data(gp_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_yyyy_y = cbuffer.data(gp_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_yyyy_z = cbuffer.data(gp_geom_01_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyy_x, g_0_x_yyy_xy, g_0_x_yyy_y, g_0_x_yyy_yy, g_0_x_yyy_yz, g_0_x_yyy_z, g_0_x_yyyy_x, g_0_x_yyyy_y, g_0_x_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyy_x[k] = -g_0_x_yyy_x[k] * ab_y + g_0_x_yyy_xy[k];

                g_0_x_yyyy_y[k] = -g_0_x_yyy_y[k] * ab_y + g_0_x_yyy_yy[k];

                g_0_x_yyyy_z[k] = -g_0_x_yyy_z[k] * ab_y + g_0_x_yyy_yz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyz_x = cbuffer.data(gp_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_yyyz_y = cbuffer.data(gp_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_yyyz_z = cbuffer.data(gp_geom_01_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyz_x, g_0_x_yyyz_y, g_0_x_yyyz_z, g_0_x_yyz_x, g_0_x_yyz_xy, g_0_x_yyz_y, g_0_x_yyz_yy, g_0_x_yyz_yz, g_0_x_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyz_x[k] = -g_0_x_yyz_x[k] * ab_y + g_0_x_yyz_xy[k];

                g_0_x_yyyz_y[k] = -g_0_x_yyz_y[k] * ab_y + g_0_x_yyz_yy[k];

                g_0_x_yyyz_z[k] = -g_0_x_yyz_z[k] * ab_y + g_0_x_yyz_yz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyzz_x = cbuffer.data(gp_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_yyzz_y = cbuffer.data(gp_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_yyzz_z = cbuffer.data(gp_geom_01_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyzz_x, g_0_x_yyzz_y, g_0_x_yyzz_z, g_0_x_yzz_x, g_0_x_yzz_xy, g_0_x_yzz_y, g_0_x_yzz_yy, g_0_x_yzz_yz, g_0_x_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyzz_x[k] = -g_0_x_yzz_x[k] * ab_y + g_0_x_yzz_xy[k];

                g_0_x_yyzz_y[k] = -g_0_x_yzz_y[k] * ab_y + g_0_x_yzz_yy[k];

                g_0_x_yyzz_z[k] = -g_0_x_yzz_z[k] * ab_y + g_0_x_yzz_yz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_0_x_yzzz_x = cbuffer.data(gp_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_yzzz_y = cbuffer.data(gp_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_yzzz_z = cbuffer.data(gp_geom_01_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yzzz_x, g_0_x_yzzz_y, g_0_x_yzzz_z, g_0_x_zzz_x, g_0_x_zzz_xy, g_0_x_zzz_y, g_0_x_zzz_yy, g_0_x_zzz_yz, g_0_x_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzzz_x[k] = -g_0_x_zzz_x[k] * ab_y + g_0_x_zzz_xy[k];

                g_0_x_yzzz_y[k] = -g_0_x_zzz_y[k] * ab_y + g_0_x_zzz_yy[k];

                g_0_x_yzzz_z[k] = -g_0_x_zzz_z[k] * ab_y + g_0_x_zzz_yz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_0_x_zzzz_x = cbuffer.data(gp_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_zzzz_y = cbuffer.data(gp_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_zzzz_z = cbuffer.data(gp_geom_01_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzz_x, g_0_x_zzz_xz, g_0_x_zzz_y, g_0_x_zzz_yz, g_0_x_zzz_z, g_0_x_zzz_zz, g_0_x_zzzz_x, g_0_x_zzzz_y, g_0_x_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzzz_x[k] = -g_0_x_zzz_x[k] * ab_z + g_0_x_zzz_xz[k];

                g_0_x_zzzz_y[k] = -g_0_x_zzz_y[k] * ab_z + g_0_x_zzz_yz[k];

                g_0_x_zzzz_z[k] = -g_0_x_zzz_z[k] * ab_z + g_0_x_zzz_zz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxx_x = cbuffer.data(gp_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_y_xxxx_y = cbuffer.data(gp_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_y_xxxx_z = cbuffer.data(gp_geom_01_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxx_x, g_0_y_xxx_xx, g_0_y_xxx_xy, g_0_y_xxx_xz, g_0_y_xxx_y, g_0_y_xxx_z, g_0_y_xxxx_x, g_0_y_xxxx_y, g_0_y_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxx_x[k] = -g_0_y_xxx_x[k] * ab_x + g_0_y_xxx_xx[k];

                g_0_y_xxxx_y[k] = -g_0_y_xxx_y[k] * ab_x + g_0_y_xxx_xy[k];

                g_0_y_xxxx_z[k] = -g_0_y_xxx_z[k] * ab_x + g_0_y_xxx_xz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxy_x = cbuffer.data(gp_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_y_xxxy_y = cbuffer.data(gp_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_y_xxxy_z = cbuffer.data(gp_geom_01_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxy_x, g_0_y_xxxy_y, g_0_y_xxxy_z, g_0_y_xxy_x, g_0_y_xxy_xx, g_0_y_xxy_xy, g_0_y_xxy_xz, g_0_y_xxy_y, g_0_y_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxy_x[k] = -g_0_y_xxy_x[k] * ab_x + g_0_y_xxy_xx[k];

                g_0_y_xxxy_y[k] = -g_0_y_xxy_y[k] * ab_x + g_0_y_xxy_xy[k];

                g_0_y_xxxy_z[k] = -g_0_y_xxy_z[k] * ab_x + g_0_y_xxy_xz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxz_x = cbuffer.data(gp_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_y_xxxz_y = cbuffer.data(gp_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_y_xxxz_z = cbuffer.data(gp_geom_01_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxz_x, g_0_y_xxxz_y, g_0_y_xxxz_z, g_0_y_xxz_x, g_0_y_xxz_xx, g_0_y_xxz_xy, g_0_y_xxz_xz, g_0_y_xxz_y, g_0_y_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxz_x[k] = -g_0_y_xxz_x[k] * ab_x + g_0_y_xxz_xx[k];

                g_0_y_xxxz_y[k] = -g_0_y_xxz_y[k] * ab_x + g_0_y_xxz_xy[k];

                g_0_y_xxxz_z[k] = -g_0_y_xxz_z[k] * ab_x + g_0_y_xxz_xz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyy_x = cbuffer.data(gp_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_y_xxyy_y = cbuffer.data(gp_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_y_xxyy_z = cbuffer.data(gp_geom_01_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyy_x, g_0_y_xxyy_y, g_0_y_xxyy_z, g_0_y_xyy_x, g_0_y_xyy_xx, g_0_y_xyy_xy, g_0_y_xyy_xz, g_0_y_xyy_y, g_0_y_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyy_x[k] = -g_0_y_xyy_x[k] * ab_x + g_0_y_xyy_xx[k];

                g_0_y_xxyy_y[k] = -g_0_y_xyy_y[k] * ab_x + g_0_y_xyy_xy[k];

                g_0_y_xxyy_z[k] = -g_0_y_xyy_z[k] * ab_x + g_0_y_xyy_xz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyz_x = cbuffer.data(gp_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_y_xxyz_y = cbuffer.data(gp_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_y_xxyz_z = cbuffer.data(gp_geom_01_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyz_x, g_0_y_xxyz_y, g_0_y_xxyz_z, g_0_y_xyz_x, g_0_y_xyz_xx, g_0_y_xyz_xy, g_0_y_xyz_xz, g_0_y_xyz_y, g_0_y_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyz_x[k] = -g_0_y_xyz_x[k] * ab_x + g_0_y_xyz_xx[k];

                g_0_y_xxyz_y[k] = -g_0_y_xyz_y[k] * ab_x + g_0_y_xyz_xy[k];

                g_0_y_xxyz_z[k] = -g_0_y_xyz_z[k] * ab_x + g_0_y_xyz_xz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxzz_x = cbuffer.data(gp_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_y_xxzz_y = cbuffer.data(gp_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_y_xxzz_z = cbuffer.data(gp_geom_01_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxzz_x, g_0_y_xxzz_y, g_0_y_xxzz_z, g_0_y_xzz_x, g_0_y_xzz_xx, g_0_y_xzz_xy, g_0_y_xzz_xz, g_0_y_xzz_y, g_0_y_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxzz_x[k] = -g_0_y_xzz_x[k] * ab_x + g_0_y_xzz_xx[k];

                g_0_y_xxzz_y[k] = -g_0_y_xzz_y[k] * ab_x + g_0_y_xzz_xy[k];

                g_0_y_xxzz_z[k] = -g_0_y_xzz_z[k] * ab_x + g_0_y_xzz_xz[k];
            }

            /// Set up 63-66 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyy_x = cbuffer.data(gp_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_y_xyyy_y = cbuffer.data(gp_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_y_xyyy_z = cbuffer.data(gp_geom_01_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyy_x, g_0_y_xyyy_y, g_0_y_xyyy_z, g_0_y_yyy_x, g_0_y_yyy_xx, g_0_y_yyy_xy, g_0_y_yyy_xz, g_0_y_yyy_y, g_0_y_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyy_x[k] = -g_0_y_yyy_x[k] * ab_x + g_0_y_yyy_xx[k];

                g_0_y_xyyy_y[k] = -g_0_y_yyy_y[k] * ab_x + g_0_y_yyy_xy[k];

                g_0_y_xyyy_z[k] = -g_0_y_yyy_z[k] * ab_x + g_0_y_yyy_xz[k];
            }

            /// Set up 66-69 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyz_x = cbuffer.data(gp_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_y_xyyz_y = cbuffer.data(gp_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_y_xyyz_z = cbuffer.data(gp_geom_01_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyz_x, g_0_y_xyyz_y, g_0_y_xyyz_z, g_0_y_yyz_x, g_0_y_yyz_xx, g_0_y_yyz_xy, g_0_y_yyz_xz, g_0_y_yyz_y, g_0_y_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyz_x[k] = -g_0_y_yyz_x[k] * ab_x + g_0_y_yyz_xx[k];

                g_0_y_xyyz_y[k] = -g_0_y_yyz_y[k] * ab_x + g_0_y_yyz_xy[k];

                g_0_y_xyyz_z[k] = -g_0_y_yyz_z[k] * ab_x + g_0_y_yyz_xz[k];
            }

            /// Set up 69-72 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyzz_x = cbuffer.data(gp_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_y_xyzz_y = cbuffer.data(gp_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_y_xyzz_z = cbuffer.data(gp_geom_01_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyzz_x, g_0_y_xyzz_y, g_0_y_xyzz_z, g_0_y_yzz_x, g_0_y_yzz_xx, g_0_y_yzz_xy, g_0_y_yzz_xz, g_0_y_yzz_y, g_0_y_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyzz_x[k] = -g_0_y_yzz_x[k] * ab_x + g_0_y_yzz_xx[k];

                g_0_y_xyzz_y[k] = -g_0_y_yzz_y[k] * ab_x + g_0_y_yzz_xy[k];

                g_0_y_xyzz_z[k] = -g_0_y_yzz_z[k] * ab_x + g_0_y_yzz_xz[k];
            }

            /// Set up 72-75 components of targeted buffer : cbuffer.data(

            auto g_0_y_xzzz_x = cbuffer.data(gp_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_y_xzzz_y = cbuffer.data(gp_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_y_xzzz_z = cbuffer.data(gp_geom_01_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xzzz_x, g_0_y_xzzz_y, g_0_y_xzzz_z, g_0_y_zzz_x, g_0_y_zzz_xx, g_0_y_zzz_xy, g_0_y_zzz_xz, g_0_y_zzz_y, g_0_y_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzzz_x[k] = -g_0_y_zzz_x[k] * ab_x + g_0_y_zzz_xx[k];

                g_0_y_xzzz_y[k] = -g_0_y_zzz_y[k] * ab_x + g_0_y_zzz_xy[k];

                g_0_y_xzzz_z[k] = -g_0_y_zzz_z[k] * ab_x + g_0_y_zzz_xz[k];
            }

            /// Set up 75-78 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyy_x = cbuffer.data(gp_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_y_yyyy_y = cbuffer.data(gp_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_y_yyyy_z = cbuffer.data(gp_geom_01_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyy_x, g_0_y_yyy_xy, g_0_y_yyy_y, g_0_y_yyy_yy, g_0_y_yyy_yz, g_0_y_yyy_z, g_0_y_yyyy_x, g_0_y_yyyy_y, g_0_y_yyyy_z, g_yyy_x, g_yyy_y, g_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyy_x[k] = g_yyy_x[k] - g_0_y_yyy_x[k] * ab_y + g_0_y_yyy_xy[k];

                g_0_y_yyyy_y[k] = g_yyy_y[k] - g_0_y_yyy_y[k] * ab_y + g_0_y_yyy_yy[k];

                g_0_y_yyyy_z[k] = g_yyy_z[k] - g_0_y_yyy_z[k] * ab_y + g_0_y_yyy_yz[k];
            }

            /// Set up 78-81 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyz_x = cbuffer.data(gp_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_y_yyyz_y = cbuffer.data(gp_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_y_yyyz_z = cbuffer.data(gp_geom_01_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyy_x, g_0_y_yyy_xz, g_0_y_yyy_y, g_0_y_yyy_yz, g_0_y_yyy_z, g_0_y_yyy_zz, g_0_y_yyyz_x, g_0_y_yyyz_y, g_0_y_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyz_x[k] = -g_0_y_yyy_x[k] * ab_z + g_0_y_yyy_xz[k];

                g_0_y_yyyz_y[k] = -g_0_y_yyy_y[k] * ab_z + g_0_y_yyy_yz[k];

                g_0_y_yyyz_z[k] = -g_0_y_yyy_z[k] * ab_z + g_0_y_yyy_zz[k];
            }

            /// Set up 81-84 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyzz_x = cbuffer.data(gp_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_y_yyzz_y = cbuffer.data(gp_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_y_yyzz_z = cbuffer.data(gp_geom_01_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyz_x, g_0_y_yyz_xz, g_0_y_yyz_y, g_0_y_yyz_yz, g_0_y_yyz_z, g_0_y_yyz_zz, g_0_y_yyzz_x, g_0_y_yyzz_y, g_0_y_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyzz_x[k] = -g_0_y_yyz_x[k] * ab_z + g_0_y_yyz_xz[k];

                g_0_y_yyzz_y[k] = -g_0_y_yyz_y[k] * ab_z + g_0_y_yyz_yz[k];

                g_0_y_yyzz_z[k] = -g_0_y_yyz_z[k] * ab_z + g_0_y_yyz_zz[k];
            }

            /// Set up 84-87 components of targeted buffer : cbuffer.data(

            auto g_0_y_yzzz_x = cbuffer.data(gp_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_y_yzzz_y = cbuffer.data(gp_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_y_yzzz_z = cbuffer.data(gp_geom_01_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yzz_x, g_0_y_yzz_xz, g_0_y_yzz_y, g_0_y_yzz_yz, g_0_y_yzz_z, g_0_y_yzz_zz, g_0_y_yzzz_x, g_0_y_yzzz_y, g_0_y_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzzz_x[k] = -g_0_y_yzz_x[k] * ab_z + g_0_y_yzz_xz[k];

                g_0_y_yzzz_y[k] = -g_0_y_yzz_y[k] * ab_z + g_0_y_yzz_yz[k];

                g_0_y_yzzz_z[k] = -g_0_y_yzz_z[k] * ab_z + g_0_y_yzz_zz[k];
            }

            /// Set up 87-90 components of targeted buffer : cbuffer.data(

            auto g_0_y_zzzz_x = cbuffer.data(gp_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_y_zzzz_y = cbuffer.data(gp_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_y_zzzz_z = cbuffer.data(gp_geom_01_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzz_x, g_0_y_zzz_xz, g_0_y_zzz_y, g_0_y_zzz_yz, g_0_y_zzz_z, g_0_y_zzz_zz, g_0_y_zzzz_x, g_0_y_zzzz_y, g_0_y_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzzz_x[k] = -g_0_y_zzz_x[k] * ab_z + g_0_y_zzz_xz[k];

                g_0_y_zzzz_y[k] = -g_0_y_zzz_y[k] * ab_z + g_0_y_zzz_yz[k];

                g_0_y_zzzz_z[k] = -g_0_y_zzz_z[k] * ab_z + g_0_y_zzz_zz[k];
            }

            /// Set up 90-93 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxx_x = cbuffer.data(gp_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_z_xxxx_y = cbuffer.data(gp_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_z_xxxx_z = cbuffer.data(gp_geom_01_off + 92 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxx_x, g_0_z_xxx_xx, g_0_z_xxx_xy, g_0_z_xxx_xz, g_0_z_xxx_y, g_0_z_xxx_z, g_0_z_xxxx_x, g_0_z_xxxx_y, g_0_z_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxx_x[k] = -g_0_z_xxx_x[k] * ab_x + g_0_z_xxx_xx[k];

                g_0_z_xxxx_y[k] = -g_0_z_xxx_y[k] * ab_x + g_0_z_xxx_xy[k];

                g_0_z_xxxx_z[k] = -g_0_z_xxx_z[k] * ab_x + g_0_z_xxx_xz[k];
            }

            /// Set up 93-96 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxy_x = cbuffer.data(gp_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_z_xxxy_y = cbuffer.data(gp_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_z_xxxy_z = cbuffer.data(gp_geom_01_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxy_x, g_0_z_xxxy_y, g_0_z_xxxy_z, g_0_z_xxy_x, g_0_z_xxy_xx, g_0_z_xxy_xy, g_0_z_xxy_xz, g_0_z_xxy_y, g_0_z_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxy_x[k] = -g_0_z_xxy_x[k] * ab_x + g_0_z_xxy_xx[k];

                g_0_z_xxxy_y[k] = -g_0_z_xxy_y[k] * ab_x + g_0_z_xxy_xy[k];

                g_0_z_xxxy_z[k] = -g_0_z_xxy_z[k] * ab_x + g_0_z_xxy_xz[k];
            }

            /// Set up 96-99 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxz_x = cbuffer.data(gp_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_z_xxxz_y = cbuffer.data(gp_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_z_xxxz_z = cbuffer.data(gp_geom_01_off + 98 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxz_x, g_0_z_xxxz_y, g_0_z_xxxz_z, g_0_z_xxz_x, g_0_z_xxz_xx, g_0_z_xxz_xy, g_0_z_xxz_xz, g_0_z_xxz_y, g_0_z_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxz_x[k] = -g_0_z_xxz_x[k] * ab_x + g_0_z_xxz_xx[k];

                g_0_z_xxxz_y[k] = -g_0_z_xxz_y[k] * ab_x + g_0_z_xxz_xy[k];

                g_0_z_xxxz_z[k] = -g_0_z_xxz_z[k] * ab_x + g_0_z_xxz_xz[k];
            }

            /// Set up 99-102 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyy_x = cbuffer.data(gp_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_z_xxyy_y = cbuffer.data(gp_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_z_xxyy_z = cbuffer.data(gp_geom_01_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyy_x, g_0_z_xxyy_y, g_0_z_xxyy_z, g_0_z_xyy_x, g_0_z_xyy_xx, g_0_z_xyy_xy, g_0_z_xyy_xz, g_0_z_xyy_y, g_0_z_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyy_x[k] = -g_0_z_xyy_x[k] * ab_x + g_0_z_xyy_xx[k];

                g_0_z_xxyy_y[k] = -g_0_z_xyy_y[k] * ab_x + g_0_z_xyy_xy[k];

                g_0_z_xxyy_z[k] = -g_0_z_xyy_z[k] * ab_x + g_0_z_xyy_xz[k];
            }

            /// Set up 102-105 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyz_x = cbuffer.data(gp_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_z_xxyz_y = cbuffer.data(gp_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_z_xxyz_z = cbuffer.data(gp_geom_01_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyz_x, g_0_z_xxyz_y, g_0_z_xxyz_z, g_0_z_xyz_x, g_0_z_xyz_xx, g_0_z_xyz_xy, g_0_z_xyz_xz, g_0_z_xyz_y, g_0_z_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyz_x[k] = -g_0_z_xyz_x[k] * ab_x + g_0_z_xyz_xx[k];

                g_0_z_xxyz_y[k] = -g_0_z_xyz_y[k] * ab_x + g_0_z_xyz_xy[k];

                g_0_z_xxyz_z[k] = -g_0_z_xyz_z[k] * ab_x + g_0_z_xyz_xz[k];
            }

            /// Set up 105-108 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxzz_x = cbuffer.data(gp_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_z_xxzz_y = cbuffer.data(gp_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_z_xxzz_z = cbuffer.data(gp_geom_01_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxzz_x, g_0_z_xxzz_y, g_0_z_xxzz_z, g_0_z_xzz_x, g_0_z_xzz_xx, g_0_z_xzz_xy, g_0_z_xzz_xz, g_0_z_xzz_y, g_0_z_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxzz_x[k] = -g_0_z_xzz_x[k] * ab_x + g_0_z_xzz_xx[k];

                g_0_z_xxzz_y[k] = -g_0_z_xzz_y[k] * ab_x + g_0_z_xzz_xy[k];

                g_0_z_xxzz_z[k] = -g_0_z_xzz_z[k] * ab_x + g_0_z_xzz_xz[k];
            }

            /// Set up 108-111 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyy_x = cbuffer.data(gp_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_z_xyyy_y = cbuffer.data(gp_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_z_xyyy_z = cbuffer.data(gp_geom_01_off + 110 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyy_x, g_0_z_xyyy_y, g_0_z_xyyy_z, g_0_z_yyy_x, g_0_z_yyy_xx, g_0_z_yyy_xy, g_0_z_yyy_xz, g_0_z_yyy_y, g_0_z_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyy_x[k] = -g_0_z_yyy_x[k] * ab_x + g_0_z_yyy_xx[k];

                g_0_z_xyyy_y[k] = -g_0_z_yyy_y[k] * ab_x + g_0_z_yyy_xy[k];

                g_0_z_xyyy_z[k] = -g_0_z_yyy_z[k] * ab_x + g_0_z_yyy_xz[k];
            }

            /// Set up 111-114 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyz_x = cbuffer.data(gp_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_z_xyyz_y = cbuffer.data(gp_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_z_xyyz_z = cbuffer.data(gp_geom_01_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyz_x, g_0_z_xyyz_y, g_0_z_xyyz_z, g_0_z_yyz_x, g_0_z_yyz_xx, g_0_z_yyz_xy, g_0_z_yyz_xz, g_0_z_yyz_y, g_0_z_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyz_x[k] = -g_0_z_yyz_x[k] * ab_x + g_0_z_yyz_xx[k];

                g_0_z_xyyz_y[k] = -g_0_z_yyz_y[k] * ab_x + g_0_z_yyz_xy[k];

                g_0_z_xyyz_z[k] = -g_0_z_yyz_z[k] * ab_x + g_0_z_yyz_xz[k];
            }

            /// Set up 114-117 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyzz_x = cbuffer.data(gp_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_z_xyzz_y = cbuffer.data(gp_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_z_xyzz_z = cbuffer.data(gp_geom_01_off + 116 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyzz_x, g_0_z_xyzz_y, g_0_z_xyzz_z, g_0_z_yzz_x, g_0_z_yzz_xx, g_0_z_yzz_xy, g_0_z_yzz_xz, g_0_z_yzz_y, g_0_z_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyzz_x[k] = -g_0_z_yzz_x[k] * ab_x + g_0_z_yzz_xx[k];

                g_0_z_xyzz_y[k] = -g_0_z_yzz_y[k] * ab_x + g_0_z_yzz_xy[k];

                g_0_z_xyzz_z[k] = -g_0_z_yzz_z[k] * ab_x + g_0_z_yzz_xz[k];
            }

            /// Set up 117-120 components of targeted buffer : cbuffer.data(

            auto g_0_z_xzzz_x = cbuffer.data(gp_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_z_xzzz_y = cbuffer.data(gp_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_z_xzzz_z = cbuffer.data(gp_geom_01_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzzz_x, g_0_z_xzzz_y, g_0_z_xzzz_z, g_0_z_zzz_x, g_0_z_zzz_xx, g_0_z_zzz_xy, g_0_z_zzz_xz, g_0_z_zzz_y, g_0_z_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzzz_x[k] = -g_0_z_zzz_x[k] * ab_x + g_0_z_zzz_xx[k];

                g_0_z_xzzz_y[k] = -g_0_z_zzz_y[k] * ab_x + g_0_z_zzz_xy[k];

                g_0_z_xzzz_z[k] = -g_0_z_zzz_z[k] * ab_x + g_0_z_zzz_xz[k];
            }

            /// Set up 120-123 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyy_x = cbuffer.data(gp_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_z_yyyy_y = cbuffer.data(gp_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_z_yyyy_z = cbuffer.data(gp_geom_01_off + 122 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyy_x, g_0_z_yyy_xy, g_0_z_yyy_y, g_0_z_yyy_yy, g_0_z_yyy_yz, g_0_z_yyy_z, g_0_z_yyyy_x, g_0_z_yyyy_y, g_0_z_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyy_x[k] = -g_0_z_yyy_x[k] * ab_y + g_0_z_yyy_xy[k];

                g_0_z_yyyy_y[k] = -g_0_z_yyy_y[k] * ab_y + g_0_z_yyy_yy[k];

                g_0_z_yyyy_z[k] = -g_0_z_yyy_z[k] * ab_y + g_0_z_yyy_yz[k];
            }

            /// Set up 123-126 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyz_x = cbuffer.data(gp_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_z_yyyz_y = cbuffer.data(gp_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_z_yyyz_z = cbuffer.data(gp_geom_01_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyz_x, g_0_z_yyyz_y, g_0_z_yyyz_z, g_0_z_yyz_x, g_0_z_yyz_xy, g_0_z_yyz_y, g_0_z_yyz_yy, g_0_z_yyz_yz, g_0_z_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyz_x[k] = -g_0_z_yyz_x[k] * ab_y + g_0_z_yyz_xy[k];

                g_0_z_yyyz_y[k] = -g_0_z_yyz_y[k] * ab_y + g_0_z_yyz_yy[k];

                g_0_z_yyyz_z[k] = -g_0_z_yyz_z[k] * ab_y + g_0_z_yyz_yz[k];
            }

            /// Set up 126-129 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyzz_x = cbuffer.data(gp_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_z_yyzz_y = cbuffer.data(gp_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_z_yyzz_z = cbuffer.data(gp_geom_01_off + 128 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyzz_x, g_0_z_yyzz_y, g_0_z_yyzz_z, g_0_z_yzz_x, g_0_z_yzz_xy, g_0_z_yzz_y, g_0_z_yzz_yy, g_0_z_yzz_yz, g_0_z_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyzz_x[k] = -g_0_z_yzz_x[k] * ab_y + g_0_z_yzz_xy[k];

                g_0_z_yyzz_y[k] = -g_0_z_yzz_y[k] * ab_y + g_0_z_yzz_yy[k];

                g_0_z_yyzz_z[k] = -g_0_z_yzz_z[k] * ab_y + g_0_z_yzz_yz[k];
            }

            /// Set up 129-132 components of targeted buffer : cbuffer.data(

            auto g_0_z_yzzz_x = cbuffer.data(gp_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_z_yzzz_y = cbuffer.data(gp_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_z_yzzz_z = cbuffer.data(gp_geom_01_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzzz_x, g_0_z_yzzz_y, g_0_z_yzzz_z, g_0_z_zzz_x, g_0_z_zzz_xy, g_0_z_zzz_y, g_0_z_zzz_yy, g_0_z_zzz_yz, g_0_z_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzzz_x[k] = -g_0_z_zzz_x[k] * ab_y + g_0_z_zzz_xy[k];

                g_0_z_yzzz_y[k] = -g_0_z_zzz_y[k] * ab_y + g_0_z_zzz_yy[k];

                g_0_z_yzzz_z[k] = -g_0_z_zzz_z[k] * ab_y + g_0_z_zzz_yz[k];
            }

            /// Set up 132-135 components of targeted buffer : cbuffer.data(

            auto g_0_z_zzzz_x = cbuffer.data(gp_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_z_zzzz_y = cbuffer.data(gp_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_z_zzzz_z = cbuffer.data(gp_geom_01_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzz_x, g_0_z_zzz_xz, g_0_z_zzz_y, g_0_z_zzz_yz, g_0_z_zzz_z, g_0_z_zzz_zz, g_0_z_zzzz_x, g_0_z_zzzz_y, g_0_z_zzzz_z, g_zzz_x, g_zzz_y, g_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzzz_x[k] = g_zzz_x[k] - g_0_z_zzz_x[k] * ab_z + g_0_z_zzz_xz[k];

                g_0_z_zzzz_y[k] = g_zzz_y[k] - g_0_z_zzz_y[k] * ab_z + g_0_z_zzz_yz[k];

                g_0_z_zzzz_z[k] = g_zzz_z[k] - g_0_z_zzz_z[k] * ab_z + g_0_z_zzz_zz[k];
            }
        }
    }
}

} // erirec namespace

