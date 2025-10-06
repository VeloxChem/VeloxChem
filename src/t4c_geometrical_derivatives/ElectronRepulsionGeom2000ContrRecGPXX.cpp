#include "ElectronRepulsionGeom2000ContrRecGPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_gpxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_gpxx,
                                            const size_t idx_geom_10_fpxx,
                                            const size_t idx_geom_20_fpxx,
                                            const size_t idx_geom_20_fdxx,
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

            /// Set up components of auxilary buffer : FPSS

            const auto fp_geom_20_off = idx_geom_20_fpxx + i * dcomps + j;

            auto g_xx_0_xxx_x = cbuffer.data(fp_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxx_y = cbuffer.data(fp_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxx_z = cbuffer.data(fp_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxy_x = cbuffer.data(fp_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxy_y = cbuffer.data(fp_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxy_z = cbuffer.data(fp_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xxz_x = cbuffer.data(fp_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxz_y = cbuffer.data(fp_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxz_z = cbuffer.data(fp_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xyy_x = cbuffer.data(fp_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xyy_y = cbuffer.data(fp_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xyy_z = cbuffer.data(fp_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xyz_x = cbuffer.data(fp_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xyz_y = cbuffer.data(fp_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xyz_z = cbuffer.data(fp_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xzz_x = cbuffer.data(fp_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xzz_y = cbuffer.data(fp_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xzz_z = cbuffer.data(fp_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_yyy_x = cbuffer.data(fp_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_yyy_y = cbuffer.data(fp_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_yyy_z = cbuffer.data(fp_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_yyz_x = cbuffer.data(fp_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_yyz_y = cbuffer.data(fp_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_yyz_z = cbuffer.data(fp_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_yzz_x = cbuffer.data(fp_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_yzz_y = cbuffer.data(fp_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_yzz_z = cbuffer.data(fp_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_zzz_x = cbuffer.data(fp_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_zzz_y = cbuffer.data(fp_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_zzz_z = cbuffer.data(fp_geom_20_off + 29 * ccomps * dcomps);

            auto g_xy_0_xxx_x = cbuffer.data(fp_geom_20_off + 30 * ccomps * dcomps);

            auto g_xy_0_xxx_y = cbuffer.data(fp_geom_20_off + 31 * ccomps * dcomps);

            auto g_xy_0_xxx_z = cbuffer.data(fp_geom_20_off + 32 * ccomps * dcomps);

            auto g_xy_0_xxy_x = cbuffer.data(fp_geom_20_off + 33 * ccomps * dcomps);

            auto g_xy_0_xxy_y = cbuffer.data(fp_geom_20_off + 34 * ccomps * dcomps);

            auto g_xy_0_xxy_z = cbuffer.data(fp_geom_20_off + 35 * ccomps * dcomps);

            auto g_xy_0_xxz_x = cbuffer.data(fp_geom_20_off + 36 * ccomps * dcomps);

            auto g_xy_0_xxz_y = cbuffer.data(fp_geom_20_off + 37 * ccomps * dcomps);

            auto g_xy_0_xxz_z = cbuffer.data(fp_geom_20_off + 38 * ccomps * dcomps);

            auto g_xy_0_xyy_x = cbuffer.data(fp_geom_20_off + 39 * ccomps * dcomps);

            auto g_xy_0_xyy_y = cbuffer.data(fp_geom_20_off + 40 * ccomps * dcomps);

            auto g_xy_0_xyy_z = cbuffer.data(fp_geom_20_off + 41 * ccomps * dcomps);

            auto g_xy_0_xyz_x = cbuffer.data(fp_geom_20_off + 42 * ccomps * dcomps);

            auto g_xy_0_xyz_y = cbuffer.data(fp_geom_20_off + 43 * ccomps * dcomps);

            auto g_xy_0_xyz_z = cbuffer.data(fp_geom_20_off + 44 * ccomps * dcomps);

            auto g_xy_0_xzz_x = cbuffer.data(fp_geom_20_off + 45 * ccomps * dcomps);

            auto g_xy_0_xzz_y = cbuffer.data(fp_geom_20_off + 46 * ccomps * dcomps);

            auto g_xy_0_xzz_z = cbuffer.data(fp_geom_20_off + 47 * ccomps * dcomps);

            auto g_xy_0_yyy_x = cbuffer.data(fp_geom_20_off + 48 * ccomps * dcomps);

            auto g_xy_0_yyy_y = cbuffer.data(fp_geom_20_off + 49 * ccomps * dcomps);

            auto g_xy_0_yyy_z = cbuffer.data(fp_geom_20_off + 50 * ccomps * dcomps);

            auto g_xy_0_yyz_x = cbuffer.data(fp_geom_20_off + 51 * ccomps * dcomps);

            auto g_xy_0_yyz_y = cbuffer.data(fp_geom_20_off + 52 * ccomps * dcomps);

            auto g_xy_0_yyz_z = cbuffer.data(fp_geom_20_off + 53 * ccomps * dcomps);

            auto g_xy_0_yzz_x = cbuffer.data(fp_geom_20_off + 54 * ccomps * dcomps);

            auto g_xy_0_yzz_y = cbuffer.data(fp_geom_20_off + 55 * ccomps * dcomps);

            auto g_xy_0_yzz_z = cbuffer.data(fp_geom_20_off + 56 * ccomps * dcomps);

            auto g_xy_0_zzz_x = cbuffer.data(fp_geom_20_off + 57 * ccomps * dcomps);

            auto g_xy_0_zzz_y = cbuffer.data(fp_geom_20_off + 58 * ccomps * dcomps);

            auto g_xy_0_zzz_z = cbuffer.data(fp_geom_20_off + 59 * ccomps * dcomps);

            auto g_xz_0_xxx_x = cbuffer.data(fp_geom_20_off + 60 * ccomps * dcomps);

            auto g_xz_0_xxx_y = cbuffer.data(fp_geom_20_off + 61 * ccomps * dcomps);

            auto g_xz_0_xxx_z = cbuffer.data(fp_geom_20_off + 62 * ccomps * dcomps);

            auto g_xz_0_xxy_x = cbuffer.data(fp_geom_20_off + 63 * ccomps * dcomps);

            auto g_xz_0_xxy_y = cbuffer.data(fp_geom_20_off + 64 * ccomps * dcomps);

            auto g_xz_0_xxy_z = cbuffer.data(fp_geom_20_off + 65 * ccomps * dcomps);

            auto g_xz_0_xxz_x = cbuffer.data(fp_geom_20_off + 66 * ccomps * dcomps);

            auto g_xz_0_xxz_y = cbuffer.data(fp_geom_20_off + 67 * ccomps * dcomps);

            auto g_xz_0_xxz_z = cbuffer.data(fp_geom_20_off + 68 * ccomps * dcomps);

            auto g_xz_0_xyy_x = cbuffer.data(fp_geom_20_off + 69 * ccomps * dcomps);

            auto g_xz_0_xyy_y = cbuffer.data(fp_geom_20_off + 70 * ccomps * dcomps);

            auto g_xz_0_xyy_z = cbuffer.data(fp_geom_20_off + 71 * ccomps * dcomps);

            auto g_xz_0_xyz_x = cbuffer.data(fp_geom_20_off + 72 * ccomps * dcomps);

            auto g_xz_0_xyz_y = cbuffer.data(fp_geom_20_off + 73 * ccomps * dcomps);

            auto g_xz_0_xyz_z = cbuffer.data(fp_geom_20_off + 74 * ccomps * dcomps);

            auto g_xz_0_xzz_x = cbuffer.data(fp_geom_20_off + 75 * ccomps * dcomps);

            auto g_xz_0_xzz_y = cbuffer.data(fp_geom_20_off + 76 * ccomps * dcomps);

            auto g_xz_0_xzz_z = cbuffer.data(fp_geom_20_off + 77 * ccomps * dcomps);

            auto g_xz_0_yyy_x = cbuffer.data(fp_geom_20_off + 78 * ccomps * dcomps);

            auto g_xz_0_yyy_y = cbuffer.data(fp_geom_20_off + 79 * ccomps * dcomps);

            auto g_xz_0_yyy_z = cbuffer.data(fp_geom_20_off + 80 * ccomps * dcomps);

            auto g_xz_0_yyz_x = cbuffer.data(fp_geom_20_off + 81 * ccomps * dcomps);

            auto g_xz_0_yyz_y = cbuffer.data(fp_geom_20_off + 82 * ccomps * dcomps);

            auto g_xz_0_yyz_z = cbuffer.data(fp_geom_20_off + 83 * ccomps * dcomps);

            auto g_xz_0_yzz_x = cbuffer.data(fp_geom_20_off + 84 * ccomps * dcomps);

            auto g_xz_0_yzz_y = cbuffer.data(fp_geom_20_off + 85 * ccomps * dcomps);

            auto g_xz_0_yzz_z = cbuffer.data(fp_geom_20_off + 86 * ccomps * dcomps);

            auto g_xz_0_zzz_x = cbuffer.data(fp_geom_20_off + 87 * ccomps * dcomps);

            auto g_xz_0_zzz_y = cbuffer.data(fp_geom_20_off + 88 * ccomps * dcomps);

            auto g_xz_0_zzz_z = cbuffer.data(fp_geom_20_off + 89 * ccomps * dcomps);

            auto g_yy_0_xxx_x = cbuffer.data(fp_geom_20_off + 90 * ccomps * dcomps);

            auto g_yy_0_xxx_y = cbuffer.data(fp_geom_20_off + 91 * ccomps * dcomps);

            auto g_yy_0_xxx_z = cbuffer.data(fp_geom_20_off + 92 * ccomps * dcomps);

            auto g_yy_0_xxy_x = cbuffer.data(fp_geom_20_off + 93 * ccomps * dcomps);

            auto g_yy_0_xxy_y = cbuffer.data(fp_geom_20_off + 94 * ccomps * dcomps);

            auto g_yy_0_xxy_z = cbuffer.data(fp_geom_20_off + 95 * ccomps * dcomps);

            auto g_yy_0_xxz_x = cbuffer.data(fp_geom_20_off + 96 * ccomps * dcomps);

            auto g_yy_0_xxz_y = cbuffer.data(fp_geom_20_off + 97 * ccomps * dcomps);

            auto g_yy_0_xxz_z = cbuffer.data(fp_geom_20_off + 98 * ccomps * dcomps);

            auto g_yy_0_xyy_x = cbuffer.data(fp_geom_20_off + 99 * ccomps * dcomps);

            auto g_yy_0_xyy_y = cbuffer.data(fp_geom_20_off + 100 * ccomps * dcomps);

            auto g_yy_0_xyy_z = cbuffer.data(fp_geom_20_off + 101 * ccomps * dcomps);

            auto g_yy_0_xyz_x = cbuffer.data(fp_geom_20_off + 102 * ccomps * dcomps);

            auto g_yy_0_xyz_y = cbuffer.data(fp_geom_20_off + 103 * ccomps * dcomps);

            auto g_yy_0_xyz_z = cbuffer.data(fp_geom_20_off + 104 * ccomps * dcomps);

            auto g_yy_0_xzz_x = cbuffer.data(fp_geom_20_off + 105 * ccomps * dcomps);

            auto g_yy_0_xzz_y = cbuffer.data(fp_geom_20_off + 106 * ccomps * dcomps);

            auto g_yy_0_xzz_z = cbuffer.data(fp_geom_20_off + 107 * ccomps * dcomps);

            auto g_yy_0_yyy_x = cbuffer.data(fp_geom_20_off + 108 * ccomps * dcomps);

            auto g_yy_0_yyy_y = cbuffer.data(fp_geom_20_off + 109 * ccomps * dcomps);

            auto g_yy_0_yyy_z = cbuffer.data(fp_geom_20_off + 110 * ccomps * dcomps);

            auto g_yy_0_yyz_x = cbuffer.data(fp_geom_20_off + 111 * ccomps * dcomps);

            auto g_yy_0_yyz_y = cbuffer.data(fp_geom_20_off + 112 * ccomps * dcomps);

            auto g_yy_0_yyz_z = cbuffer.data(fp_geom_20_off + 113 * ccomps * dcomps);

            auto g_yy_0_yzz_x = cbuffer.data(fp_geom_20_off + 114 * ccomps * dcomps);

            auto g_yy_0_yzz_y = cbuffer.data(fp_geom_20_off + 115 * ccomps * dcomps);

            auto g_yy_0_yzz_z = cbuffer.data(fp_geom_20_off + 116 * ccomps * dcomps);

            auto g_yy_0_zzz_x = cbuffer.data(fp_geom_20_off + 117 * ccomps * dcomps);

            auto g_yy_0_zzz_y = cbuffer.data(fp_geom_20_off + 118 * ccomps * dcomps);

            auto g_yy_0_zzz_z = cbuffer.data(fp_geom_20_off + 119 * ccomps * dcomps);

            auto g_yz_0_xxx_x = cbuffer.data(fp_geom_20_off + 120 * ccomps * dcomps);

            auto g_yz_0_xxx_y = cbuffer.data(fp_geom_20_off + 121 * ccomps * dcomps);

            auto g_yz_0_xxx_z = cbuffer.data(fp_geom_20_off + 122 * ccomps * dcomps);

            auto g_yz_0_xxy_x = cbuffer.data(fp_geom_20_off + 123 * ccomps * dcomps);

            auto g_yz_0_xxy_y = cbuffer.data(fp_geom_20_off + 124 * ccomps * dcomps);

            auto g_yz_0_xxy_z = cbuffer.data(fp_geom_20_off + 125 * ccomps * dcomps);

            auto g_yz_0_xxz_x = cbuffer.data(fp_geom_20_off + 126 * ccomps * dcomps);

            auto g_yz_0_xxz_y = cbuffer.data(fp_geom_20_off + 127 * ccomps * dcomps);

            auto g_yz_0_xxz_z = cbuffer.data(fp_geom_20_off + 128 * ccomps * dcomps);

            auto g_yz_0_xyy_x = cbuffer.data(fp_geom_20_off + 129 * ccomps * dcomps);

            auto g_yz_0_xyy_y = cbuffer.data(fp_geom_20_off + 130 * ccomps * dcomps);

            auto g_yz_0_xyy_z = cbuffer.data(fp_geom_20_off + 131 * ccomps * dcomps);

            auto g_yz_0_xyz_x = cbuffer.data(fp_geom_20_off + 132 * ccomps * dcomps);

            auto g_yz_0_xyz_y = cbuffer.data(fp_geom_20_off + 133 * ccomps * dcomps);

            auto g_yz_0_xyz_z = cbuffer.data(fp_geom_20_off + 134 * ccomps * dcomps);

            auto g_yz_0_xzz_x = cbuffer.data(fp_geom_20_off + 135 * ccomps * dcomps);

            auto g_yz_0_xzz_y = cbuffer.data(fp_geom_20_off + 136 * ccomps * dcomps);

            auto g_yz_0_xzz_z = cbuffer.data(fp_geom_20_off + 137 * ccomps * dcomps);

            auto g_yz_0_yyy_x = cbuffer.data(fp_geom_20_off + 138 * ccomps * dcomps);

            auto g_yz_0_yyy_y = cbuffer.data(fp_geom_20_off + 139 * ccomps * dcomps);

            auto g_yz_0_yyy_z = cbuffer.data(fp_geom_20_off + 140 * ccomps * dcomps);

            auto g_yz_0_yyz_x = cbuffer.data(fp_geom_20_off + 141 * ccomps * dcomps);

            auto g_yz_0_yyz_y = cbuffer.data(fp_geom_20_off + 142 * ccomps * dcomps);

            auto g_yz_0_yyz_z = cbuffer.data(fp_geom_20_off + 143 * ccomps * dcomps);

            auto g_yz_0_yzz_x = cbuffer.data(fp_geom_20_off + 144 * ccomps * dcomps);

            auto g_yz_0_yzz_y = cbuffer.data(fp_geom_20_off + 145 * ccomps * dcomps);

            auto g_yz_0_yzz_z = cbuffer.data(fp_geom_20_off + 146 * ccomps * dcomps);

            auto g_yz_0_zzz_x = cbuffer.data(fp_geom_20_off + 147 * ccomps * dcomps);

            auto g_yz_0_zzz_y = cbuffer.data(fp_geom_20_off + 148 * ccomps * dcomps);

            auto g_yz_0_zzz_z = cbuffer.data(fp_geom_20_off + 149 * ccomps * dcomps);

            auto g_zz_0_xxx_x = cbuffer.data(fp_geom_20_off + 150 * ccomps * dcomps);

            auto g_zz_0_xxx_y = cbuffer.data(fp_geom_20_off + 151 * ccomps * dcomps);

            auto g_zz_0_xxx_z = cbuffer.data(fp_geom_20_off + 152 * ccomps * dcomps);

            auto g_zz_0_xxy_x = cbuffer.data(fp_geom_20_off + 153 * ccomps * dcomps);

            auto g_zz_0_xxy_y = cbuffer.data(fp_geom_20_off + 154 * ccomps * dcomps);

            auto g_zz_0_xxy_z = cbuffer.data(fp_geom_20_off + 155 * ccomps * dcomps);

            auto g_zz_0_xxz_x = cbuffer.data(fp_geom_20_off + 156 * ccomps * dcomps);

            auto g_zz_0_xxz_y = cbuffer.data(fp_geom_20_off + 157 * ccomps * dcomps);

            auto g_zz_0_xxz_z = cbuffer.data(fp_geom_20_off + 158 * ccomps * dcomps);

            auto g_zz_0_xyy_x = cbuffer.data(fp_geom_20_off + 159 * ccomps * dcomps);

            auto g_zz_0_xyy_y = cbuffer.data(fp_geom_20_off + 160 * ccomps * dcomps);

            auto g_zz_0_xyy_z = cbuffer.data(fp_geom_20_off + 161 * ccomps * dcomps);

            auto g_zz_0_xyz_x = cbuffer.data(fp_geom_20_off + 162 * ccomps * dcomps);

            auto g_zz_0_xyz_y = cbuffer.data(fp_geom_20_off + 163 * ccomps * dcomps);

            auto g_zz_0_xyz_z = cbuffer.data(fp_geom_20_off + 164 * ccomps * dcomps);

            auto g_zz_0_xzz_x = cbuffer.data(fp_geom_20_off + 165 * ccomps * dcomps);

            auto g_zz_0_xzz_y = cbuffer.data(fp_geom_20_off + 166 * ccomps * dcomps);

            auto g_zz_0_xzz_z = cbuffer.data(fp_geom_20_off + 167 * ccomps * dcomps);

            auto g_zz_0_yyy_x = cbuffer.data(fp_geom_20_off + 168 * ccomps * dcomps);

            auto g_zz_0_yyy_y = cbuffer.data(fp_geom_20_off + 169 * ccomps * dcomps);

            auto g_zz_0_yyy_z = cbuffer.data(fp_geom_20_off + 170 * ccomps * dcomps);

            auto g_zz_0_yyz_x = cbuffer.data(fp_geom_20_off + 171 * ccomps * dcomps);

            auto g_zz_0_yyz_y = cbuffer.data(fp_geom_20_off + 172 * ccomps * dcomps);

            auto g_zz_0_yyz_z = cbuffer.data(fp_geom_20_off + 173 * ccomps * dcomps);

            auto g_zz_0_yzz_x = cbuffer.data(fp_geom_20_off + 174 * ccomps * dcomps);

            auto g_zz_0_yzz_y = cbuffer.data(fp_geom_20_off + 175 * ccomps * dcomps);

            auto g_zz_0_yzz_z = cbuffer.data(fp_geom_20_off + 176 * ccomps * dcomps);

            auto g_zz_0_zzz_x = cbuffer.data(fp_geom_20_off + 177 * ccomps * dcomps);

            auto g_zz_0_zzz_y = cbuffer.data(fp_geom_20_off + 178 * ccomps * dcomps);

            auto g_zz_0_zzz_z = cbuffer.data(fp_geom_20_off + 179 * ccomps * dcomps);

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

            /// set up bra offset for contr_buffer_gpxx

            const auto gp_geom_20_off = idx_geom_20_gpxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxx_x = cbuffer.data(gp_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxxx_y = cbuffer.data(gp_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxxx_z = cbuffer.data(gp_geom_20_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxx_x, g_x_0_xxx_y, g_x_0_xxx_z, g_xx_0_xxx_x, g_xx_0_xxx_xx, g_xx_0_xxx_xy, g_xx_0_xxx_xz, g_xx_0_xxx_y, g_xx_0_xxx_z, g_xx_0_xxxx_x, g_xx_0_xxxx_y, g_xx_0_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxx_x[k] = -2.0 * g_x_0_xxx_x[k] - g_xx_0_xxx_x[k] * ab_x + g_xx_0_xxx_xx[k];

                g_xx_0_xxxx_y[k] = -2.0 * g_x_0_xxx_y[k] - g_xx_0_xxx_y[k] * ab_x + g_xx_0_xxx_xy[k];

                g_xx_0_xxxx_z[k] = -2.0 * g_x_0_xxx_z[k] - g_xx_0_xxx_z[k] * ab_x + g_xx_0_xxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxy_x = cbuffer.data(gp_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxxy_y = cbuffer.data(gp_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxxy_z = cbuffer.data(gp_geom_20_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxx_x, g_xx_0_xxx_xy, g_xx_0_xxx_y, g_xx_0_xxx_yy, g_xx_0_xxx_yz, g_xx_0_xxx_z, g_xx_0_xxxy_x, g_xx_0_xxxy_y, g_xx_0_xxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxy_x[k] = -g_xx_0_xxx_x[k] * ab_y + g_xx_0_xxx_xy[k];

                g_xx_0_xxxy_y[k] = -g_xx_0_xxx_y[k] * ab_y + g_xx_0_xxx_yy[k];

                g_xx_0_xxxy_z[k] = -g_xx_0_xxx_z[k] * ab_y + g_xx_0_xxx_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxz_x = cbuffer.data(gp_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxxz_y = cbuffer.data(gp_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxxz_z = cbuffer.data(gp_geom_20_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxx_x, g_xx_0_xxx_xz, g_xx_0_xxx_y, g_xx_0_xxx_yz, g_xx_0_xxx_z, g_xx_0_xxx_zz, g_xx_0_xxxz_x, g_xx_0_xxxz_y, g_xx_0_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxz_x[k] = -g_xx_0_xxx_x[k] * ab_z + g_xx_0_xxx_xz[k];

                g_xx_0_xxxz_y[k] = -g_xx_0_xxx_y[k] * ab_z + g_xx_0_xxx_yz[k];

                g_xx_0_xxxz_z[k] = -g_xx_0_xxx_z[k] * ab_z + g_xx_0_xxx_zz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyy_x = cbuffer.data(gp_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxyy_y = cbuffer.data(gp_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxyy_z = cbuffer.data(gp_geom_20_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxy_x, g_xx_0_xxy_xy, g_xx_0_xxy_y, g_xx_0_xxy_yy, g_xx_0_xxy_yz, g_xx_0_xxy_z, g_xx_0_xxyy_x, g_xx_0_xxyy_y, g_xx_0_xxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyy_x[k] = -g_xx_0_xxy_x[k] * ab_y + g_xx_0_xxy_xy[k];

                g_xx_0_xxyy_y[k] = -g_xx_0_xxy_y[k] * ab_y + g_xx_0_xxy_yy[k];

                g_xx_0_xxyy_z[k] = -g_xx_0_xxy_z[k] * ab_y + g_xx_0_xxy_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyz_x = cbuffer.data(gp_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxyz_y = cbuffer.data(gp_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxyz_z = cbuffer.data(gp_geom_20_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyz_x, g_xx_0_xxyz_y, g_xx_0_xxyz_z, g_xx_0_xxz_x, g_xx_0_xxz_xy, g_xx_0_xxz_y, g_xx_0_xxz_yy, g_xx_0_xxz_yz, g_xx_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyz_x[k] = -g_xx_0_xxz_x[k] * ab_y + g_xx_0_xxz_xy[k];

                g_xx_0_xxyz_y[k] = -g_xx_0_xxz_y[k] * ab_y + g_xx_0_xxz_yy[k];

                g_xx_0_xxyz_z[k] = -g_xx_0_xxz_z[k] * ab_y + g_xx_0_xxz_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxzz_x = cbuffer.data(gp_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxzz_y = cbuffer.data(gp_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxzz_z = cbuffer.data(gp_geom_20_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxz_x, g_xx_0_xxz_xz, g_xx_0_xxz_y, g_xx_0_xxz_yz, g_xx_0_xxz_z, g_xx_0_xxz_zz, g_xx_0_xxzz_x, g_xx_0_xxzz_y, g_xx_0_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxzz_x[k] = -g_xx_0_xxz_x[k] * ab_z + g_xx_0_xxz_xz[k];

                g_xx_0_xxzz_y[k] = -g_xx_0_xxz_y[k] * ab_z + g_xx_0_xxz_yz[k];

                g_xx_0_xxzz_z[k] = -g_xx_0_xxz_z[k] * ab_z + g_xx_0_xxz_zz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyy_x = cbuffer.data(gp_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xyyy_y = cbuffer.data(gp_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xyyy_z = cbuffer.data(gp_geom_20_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyy_x, g_xx_0_xyy_xy, g_xx_0_xyy_y, g_xx_0_xyy_yy, g_xx_0_xyy_yz, g_xx_0_xyy_z, g_xx_0_xyyy_x, g_xx_0_xyyy_y, g_xx_0_xyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyy_x[k] = -g_xx_0_xyy_x[k] * ab_y + g_xx_0_xyy_xy[k];

                g_xx_0_xyyy_y[k] = -g_xx_0_xyy_y[k] * ab_y + g_xx_0_xyy_yy[k];

                g_xx_0_xyyy_z[k] = -g_xx_0_xyy_z[k] * ab_y + g_xx_0_xyy_yz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyz_x = cbuffer.data(gp_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xyyz_y = cbuffer.data(gp_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xyyz_z = cbuffer.data(gp_geom_20_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyz_x, g_xx_0_xyyz_y, g_xx_0_xyyz_z, g_xx_0_xyz_x, g_xx_0_xyz_xy, g_xx_0_xyz_y, g_xx_0_xyz_yy, g_xx_0_xyz_yz, g_xx_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyz_x[k] = -g_xx_0_xyz_x[k] * ab_y + g_xx_0_xyz_xy[k];

                g_xx_0_xyyz_y[k] = -g_xx_0_xyz_y[k] * ab_y + g_xx_0_xyz_yy[k];

                g_xx_0_xyyz_z[k] = -g_xx_0_xyz_z[k] * ab_y + g_xx_0_xyz_yz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyzz_x = cbuffer.data(gp_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xyzz_y = cbuffer.data(gp_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xyzz_z = cbuffer.data(gp_geom_20_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyzz_x, g_xx_0_xyzz_y, g_xx_0_xyzz_z, g_xx_0_xzz_x, g_xx_0_xzz_xy, g_xx_0_xzz_y, g_xx_0_xzz_yy, g_xx_0_xzz_yz, g_xx_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyzz_x[k] = -g_xx_0_xzz_x[k] * ab_y + g_xx_0_xzz_xy[k];

                g_xx_0_xyzz_y[k] = -g_xx_0_xzz_y[k] * ab_y + g_xx_0_xzz_yy[k];

                g_xx_0_xyzz_z[k] = -g_xx_0_xzz_z[k] * ab_y + g_xx_0_xzz_yz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xzzz_x = cbuffer.data(gp_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xzzz_y = cbuffer.data(gp_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xzzz_z = cbuffer.data(gp_geom_20_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xzz_x, g_xx_0_xzz_xz, g_xx_0_xzz_y, g_xx_0_xzz_yz, g_xx_0_xzz_z, g_xx_0_xzz_zz, g_xx_0_xzzz_x, g_xx_0_xzzz_y, g_xx_0_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xzzz_x[k] = -g_xx_0_xzz_x[k] * ab_z + g_xx_0_xzz_xz[k];

                g_xx_0_xzzz_y[k] = -g_xx_0_xzz_y[k] * ab_z + g_xx_0_xzz_yz[k];

                g_xx_0_xzzz_z[k] = -g_xx_0_xzz_z[k] * ab_z + g_xx_0_xzz_zz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyy_x = cbuffer.data(gp_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_yyyy_y = cbuffer.data(gp_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_yyyy_z = cbuffer.data(gp_geom_20_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyy_x, g_xx_0_yyy_xy, g_xx_0_yyy_y, g_xx_0_yyy_yy, g_xx_0_yyy_yz, g_xx_0_yyy_z, g_xx_0_yyyy_x, g_xx_0_yyyy_y, g_xx_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyy_x[k] = -g_xx_0_yyy_x[k] * ab_y + g_xx_0_yyy_xy[k];

                g_xx_0_yyyy_y[k] = -g_xx_0_yyy_y[k] * ab_y + g_xx_0_yyy_yy[k];

                g_xx_0_yyyy_z[k] = -g_xx_0_yyy_z[k] * ab_y + g_xx_0_yyy_yz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyz_x = cbuffer.data(gp_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_yyyz_y = cbuffer.data(gp_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_yyyz_z = cbuffer.data(gp_geom_20_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyz_x, g_xx_0_yyyz_y, g_xx_0_yyyz_z, g_xx_0_yyz_x, g_xx_0_yyz_xy, g_xx_0_yyz_y, g_xx_0_yyz_yy, g_xx_0_yyz_yz, g_xx_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyz_x[k] = -g_xx_0_yyz_x[k] * ab_y + g_xx_0_yyz_xy[k];

                g_xx_0_yyyz_y[k] = -g_xx_0_yyz_y[k] * ab_y + g_xx_0_yyz_yy[k];

                g_xx_0_yyyz_z[k] = -g_xx_0_yyz_z[k] * ab_y + g_xx_0_yyz_yz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyzz_x = cbuffer.data(gp_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_yyzz_y = cbuffer.data(gp_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_yyzz_z = cbuffer.data(gp_geom_20_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyzz_x, g_xx_0_yyzz_y, g_xx_0_yyzz_z, g_xx_0_yzz_x, g_xx_0_yzz_xy, g_xx_0_yzz_y, g_xx_0_yzz_yy, g_xx_0_yzz_yz, g_xx_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyzz_x[k] = -g_xx_0_yzz_x[k] * ab_y + g_xx_0_yzz_xy[k];

                g_xx_0_yyzz_y[k] = -g_xx_0_yzz_y[k] * ab_y + g_xx_0_yzz_yy[k];

                g_xx_0_yyzz_z[k] = -g_xx_0_yzz_z[k] * ab_y + g_xx_0_yzz_yz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yzzz_x = cbuffer.data(gp_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_yzzz_y = cbuffer.data(gp_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_yzzz_z = cbuffer.data(gp_geom_20_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yzzz_x, g_xx_0_yzzz_y, g_xx_0_yzzz_z, g_xx_0_zzz_x, g_xx_0_zzz_xy, g_xx_0_zzz_y, g_xx_0_zzz_yy, g_xx_0_zzz_yz, g_xx_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yzzz_x[k] = -g_xx_0_zzz_x[k] * ab_y + g_xx_0_zzz_xy[k];

                g_xx_0_yzzz_y[k] = -g_xx_0_zzz_y[k] * ab_y + g_xx_0_zzz_yy[k];

                g_xx_0_yzzz_z[k] = -g_xx_0_zzz_z[k] * ab_y + g_xx_0_zzz_yz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_xx_0_zzzz_x = cbuffer.data(gp_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_zzzz_y = cbuffer.data(gp_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_zzzz_z = cbuffer.data(gp_geom_20_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_zzz_x, g_xx_0_zzz_xz, g_xx_0_zzz_y, g_xx_0_zzz_yz, g_xx_0_zzz_z, g_xx_0_zzz_zz, g_xx_0_zzzz_x, g_xx_0_zzzz_y, g_xx_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_zzzz_x[k] = -g_xx_0_zzz_x[k] * ab_z + g_xx_0_zzz_xz[k];

                g_xx_0_zzzz_y[k] = -g_xx_0_zzz_y[k] * ab_z + g_xx_0_zzz_yz[k];

                g_xx_0_zzzz_z[k] = -g_xx_0_zzz_z[k] * ab_z + g_xx_0_zzz_zz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxx_x = cbuffer.data(gp_geom_20_off + 45 * ccomps * dcomps);

            auto g_xy_0_xxxx_y = cbuffer.data(gp_geom_20_off + 46 * ccomps * dcomps);

            auto g_xy_0_xxxx_z = cbuffer.data(gp_geom_20_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxx_x, g_xy_0_xxx_xx, g_xy_0_xxx_xy, g_xy_0_xxx_xz, g_xy_0_xxx_y, g_xy_0_xxx_z, g_xy_0_xxxx_x, g_xy_0_xxxx_y, g_xy_0_xxxx_z, g_y_0_xxx_x, g_y_0_xxx_y, g_y_0_xxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxx_x[k] = -g_y_0_xxx_x[k] - g_xy_0_xxx_x[k] * ab_x + g_xy_0_xxx_xx[k];

                g_xy_0_xxxx_y[k] = -g_y_0_xxx_y[k] - g_xy_0_xxx_y[k] * ab_x + g_xy_0_xxx_xy[k];

                g_xy_0_xxxx_z[k] = -g_y_0_xxx_z[k] - g_xy_0_xxx_z[k] * ab_x + g_xy_0_xxx_xz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxy_x = cbuffer.data(gp_geom_20_off + 48 * ccomps * dcomps);

            auto g_xy_0_xxxy_y = cbuffer.data(gp_geom_20_off + 49 * ccomps * dcomps);

            auto g_xy_0_xxxy_z = cbuffer.data(gp_geom_20_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxy_x, g_xy_0_xxxy_y, g_xy_0_xxxy_z, g_xy_0_xxy_x, g_xy_0_xxy_xx, g_xy_0_xxy_xy, g_xy_0_xxy_xz, g_xy_0_xxy_y, g_xy_0_xxy_z, g_y_0_xxy_x, g_y_0_xxy_y, g_y_0_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxy_x[k] = -g_y_0_xxy_x[k] - g_xy_0_xxy_x[k] * ab_x + g_xy_0_xxy_xx[k];

                g_xy_0_xxxy_y[k] = -g_y_0_xxy_y[k] - g_xy_0_xxy_y[k] * ab_x + g_xy_0_xxy_xy[k];

                g_xy_0_xxxy_z[k] = -g_y_0_xxy_z[k] - g_xy_0_xxy_z[k] * ab_x + g_xy_0_xxy_xz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxz_x = cbuffer.data(gp_geom_20_off + 51 * ccomps * dcomps);

            auto g_xy_0_xxxz_y = cbuffer.data(gp_geom_20_off + 52 * ccomps * dcomps);

            auto g_xy_0_xxxz_z = cbuffer.data(gp_geom_20_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxx_x, g_xy_0_xxx_xz, g_xy_0_xxx_y, g_xy_0_xxx_yz, g_xy_0_xxx_z, g_xy_0_xxx_zz, g_xy_0_xxxz_x, g_xy_0_xxxz_y, g_xy_0_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxz_x[k] = -g_xy_0_xxx_x[k] * ab_z + g_xy_0_xxx_xz[k];

                g_xy_0_xxxz_y[k] = -g_xy_0_xxx_y[k] * ab_z + g_xy_0_xxx_yz[k];

                g_xy_0_xxxz_z[k] = -g_xy_0_xxx_z[k] * ab_z + g_xy_0_xxx_zz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyy_x = cbuffer.data(gp_geom_20_off + 54 * ccomps * dcomps);

            auto g_xy_0_xxyy_y = cbuffer.data(gp_geom_20_off + 55 * ccomps * dcomps);

            auto g_xy_0_xxyy_z = cbuffer.data(gp_geom_20_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyy_x, g_xy_0_xxyy_y, g_xy_0_xxyy_z, g_xy_0_xyy_x, g_xy_0_xyy_xx, g_xy_0_xyy_xy, g_xy_0_xyy_xz, g_xy_0_xyy_y, g_xy_0_xyy_z, g_y_0_xyy_x, g_y_0_xyy_y, g_y_0_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyy_x[k] = -g_y_0_xyy_x[k] - g_xy_0_xyy_x[k] * ab_x + g_xy_0_xyy_xx[k];

                g_xy_0_xxyy_y[k] = -g_y_0_xyy_y[k] - g_xy_0_xyy_y[k] * ab_x + g_xy_0_xyy_xy[k];

                g_xy_0_xxyy_z[k] = -g_y_0_xyy_z[k] - g_xy_0_xyy_z[k] * ab_x + g_xy_0_xyy_xz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyz_x = cbuffer.data(gp_geom_20_off + 57 * ccomps * dcomps);

            auto g_xy_0_xxyz_y = cbuffer.data(gp_geom_20_off + 58 * ccomps * dcomps);

            auto g_xy_0_xxyz_z = cbuffer.data(gp_geom_20_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxy_x, g_xy_0_xxy_xz, g_xy_0_xxy_y, g_xy_0_xxy_yz, g_xy_0_xxy_z, g_xy_0_xxy_zz, g_xy_0_xxyz_x, g_xy_0_xxyz_y, g_xy_0_xxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyz_x[k] = -g_xy_0_xxy_x[k] * ab_z + g_xy_0_xxy_xz[k];

                g_xy_0_xxyz_y[k] = -g_xy_0_xxy_y[k] * ab_z + g_xy_0_xxy_yz[k];

                g_xy_0_xxyz_z[k] = -g_xy_0_xxy_z[k] * ab_z + g_xy_0_xxy_zz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxzz_x = cbuffer.data(gp_geom_20_off + 60 * ccomps * dcomps);

            auto g_xy_0_xxzz_y = cbuffer.data(gp_geom_20_off + 61 * ccomps * dcomps);

            auto g_xy_0_xxzz_z = cbuffer.data(gp_geom_20_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxz_x, g_xy_0_xxz_xz, g_xy_0_xxz_y, g_xy_0_xxz_yz, g_xy_0_xxz_z, g_xy_0_xxz_zz, g_xy_0_xxzz_x, g_xy_0_xxzz_y, g_xy_0_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxzz_x[k] = -g_xy_0_xxz_x[k] * ab_z + g_xy_0_xxz_xz[k];

                g_xy_0_xxzz_y[k] = -g_xy_0_xxz_y[k] * ab_z + g_xy_0_xxz_yz[k];

                g_xy_0_xxzz_z[k] = -g_xy_0_xxz_z[k] * ab_z + g_xy_0_xxz_zz[k];
            }

            /// Set up 63-66 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyy_x = cbuffer.data(gp_geom_20_off + 63 * ccomps * dcomps);

            auto g_xy_0_xyyy_y = cbuffer.data(gp_geom_20_off + 64 * ccomps * dcomps);

            auto g_xy_0_xyyy_z = cbuffer.data(gp_geom_20_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyy_x, g_xy_0_xyyy_y, g_xy_0_xyyy_z, g_xy_0_yyy_x, g_xy_0_yyy_xx, g_xy_0_yyy_xy, g_xy_0_yyy_xz, g_xy_0_yyy_y, g_xy_0_yyy_z, g_y_0_yyy_x, g_y_0_yyy_y, g_y_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyy_x[k] = -g_y_0_yyy_x[k] - g_xy_0_yyy_x[k] * ab_x + g_xy_0_yyy_xx[k];

                g_xy_0_xyyy_y[k] = -g_y_0_yyy_y[k] - g_xy_0_yyy_y[k] * ab_x + g_xy_0_yyy_xy[k];

                g_xy_0_xyyy_z[k] = -g_y_0_yyy_z[k] - g_xy_0_yyy_z[k] * ab_x + g_xy_0_yyy_xz[k];
            }

            /// Set up 66-69 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyz_x = cbuffer.data(gp_geom_20_off + 66 * ccomps * dcomps);

            auto g_xy_0_xyyz_y = cbuffer.data(gp_geom_20_off + 67 * ccomps * dcomps);

            auto g_xy_0_xyyz_z = cbuffer.data(gp_geom_20_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyy_x, g_xy_0_xyy_xz, g_xy_0_xyy_y, g_xy_0_xyy_yz, g_xy_0_xyy_z, g_xy_0_xyy_zz, g_xy_0_xyyz_x, g_xy_0_xyyz_y, g_xy_0_xyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyz_x[k] = -g_xy_0_xyy_x[k] * ab_z + g_xy_0_xyy_xz[k];

                g_xy_0_xyyz_y[k] = -g_xy_0_xyy_y[k] * ab_z + g_xy_0_xyy_yz[k];

                g_xy_0_xyyz_z[k] = -g_xy_0_xyy_z[k] * ab_z + g_xy_0_xyy_zz[k];
            }

            /// Set up 69-72 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyzz_x = cbuffer.data(gp_geom_20_off + 69 * ccomps * dcomps);

            auto g_xy_0_xyzz_y = cbuffer.data(gp_geom_20_off + 70 * ccomps * dcomps);

            auto g_xy_0_xyzz_z = cbuffer.data(gp_geom_20_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyz_x, g_xy_0_xyz_xz, g_xy_0_xyz_y, g_xy_0_xyz_yz, g_xy_0_xyz_z, g_xy_0_xyz_zz, g_xy_0_xyzz_x, g_xy_0_xyzz_y, g_xy_0_xyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyzz_x[k] = -g_xy_0_xyz_x[k] * ab_z + g_xy_0_xyz_xz[k];

                g_xy_0_xyzz_y[k] = -g_xy_0_xyz_y[k] * ab_z + g_xy_0_xyz_yz[k];

                g_xy_0_xyzz_z[k] = -g_xy_0_xyz_z[k] * ab_z + g_xy_0_xyz_zz[k];
            }

            /// Set up 72-75 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xzzz_x = cbuffer.data(gp_geom_20_off + 72 * ccomps * dcomps);

            auto g_xy_0_xzzz_y = cbuffer.data(gp_geom_20_off + 73 * ccomps * dcomps);

            auto g_xy_0_xzzz_z = cbuffer.data(gp_geom_20_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xzz_x, g_xy_0_xzz_xz, g_xy_0_xzz_y, g_xy_0_xzz_yz, g_xy_0_xzz_z, g_xy_0_xzz_zz, g_xy_0_xzzz_x, g_xy_0_xzzz_y, g_xy_0_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xzzz_x[k] = -g_xy_0_xzz_x[k] * ab_z + g_xy_0_xzz_xz[k];

                g_xy_0_xzzz_y[k] = -g_xy_0_xzz_y[k] * ab_z + g_xy_0_xzz_yz[k];

                g_xy_0_xzzz_z[k] = -g_xy_0_xzz_z[k] * ab_z + g_xy_0_xzz_zz[k];
            }

            /// Set up 75-78 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyy_x = cbuffer.data(gp_geom_20_off + 75 * ccomps * dcomps);

            auto g_xy_0_yyyy_y = cbuffer.data(gp_geom_20_off + 76 * ccomps * dcomps);

            auto g_xy_0_yyyy_z = cbuffer.data(gp_geom_20_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyy_x, g_x_0_yyy_y, g_x_0_yyy_z, g_xy_0_yyy_x, g_xy_0_yyy_xy, g_xy_0_yyy_y, g_xy_0_yyy_yy, g_xy_0_yyy_yz, g_xy_0_yyy_z, g_xy_0_yyyy_x, g_xy_0_yyyy_y, g_xy_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyy_x[k] = -g_x_0_yyy_x[k] - g_xy_0_yyy_x[k] * ab_y + g_xy_0_yyy_xy[k];

                g_xy_0_yyyy_y[k] = -g_x_0_yyy_y[k] - g_xy_0_yyy_y[k] * ab_y + g_xy_0_yyy_yy[k];

                g_xy_0_yyyy_z[k] = -g_x_0_yyy_z[k] - g_xy_0_yyy_z[k] * ab_y + g_xy_0_yyy_yz[k];
            }

            /// Set up 78-81 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyz_x = cbuffer.data(gp_geom_20_off + 78 * ccomps * dcomps);

            auto g_xy_0_yyyz_y = cbuffer.data(gp_geom_20_off + 79 * ccomps * dcomps);

            auto g_xy_0_yyyz_z = cbuffer.data(gp_geom_20_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyy_x, g_xy_0_yyy_xz, g_xy_0_yyy_y, g_xy_0_yyy_yz, g_xy_0_yyy_z, g_xy_0_yyy_zz, g_xy_0_yyyz_x, g_xy_0_yyyz_y, g_xy_0_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyz_x[k] = -g_xy_0_yyy_x[k] * ab_z + g_xy_0_yyy_xz[k];

                g_xy_0_yyyz_y[k] = -g_xy_0_yyy_y[k] * ab_z + g_xy_0_yyy_yz[k];

                g_xy_0_yyyz_z[k] = -g_xy_0_yyy_z[k] * ab_z + g_xy_0_yyy_zz[k];
            }

            /// Set up 81-84 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyzz_x = cbuffer.data(gp_geom_20_off + 81 * ccomps * dcomps);

            auto g_xy_0_yyzz_y = cbuffer.data(gp_geom_20_off + 82 * ccomps * dcomps);

            auto g_xy_0_yyzz_z = cbuffer.data(gp_geom_20_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyz_x, g_xy_0_yyz_xz, g_xy_0_yyz_y, g_xy_0_yyz_yz, g_xy_0_yyz_z, g_xy_0_yyz_zz, g_xy_0_yyzz_x, g_xy_0_yyzz_y, g_xy_0_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyzz_x[k] = -g_xy_0_yyz_x[k] * ab_z + g_xy_0_yyz_xz[k];

                g_xy_0_yyzz_y[k] = -g_xy_0_yyz_y[k] * ab_z + g_xy_0_yyz_yz[k];

                g_xy_0_yyzz_z[k] = -g_xy_0_yyz_z[k] * ab_z + g_xy_0_yyz_zz[k];
            }

            /// Set up 84-87 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yzzz_x = cbuffer.data(gp_geom_20_off + 84 * ccomps * dcomps);

            auto g_xy_0_yzzz_y = cbuffer.data(gp_geom_20_off + 85 * ccomps * dcomps);

            auto g_xy_0_yzzz_z = cbuffer.data(gp_geom_20_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yzz_x, g_xy_0_yzz_xz, g_xy_0_yzz_y, g_xy_0_yzz_yz, g_xy_0_yzz_z, g_xy_0_yzz_zz, g_xy_0_yzzz_x, g_xy_0_yzzz_y, g_xy_0_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yzzz_x[k] = -g_xy_0_yzz_x[k] * ab_z + g_xy_0_yzz_xz[k];

                g_xy_0_yzzz_y[k] = -g_xy_0_yzz_y[k] * ab_z + g_xy_0_yzz_yz[k];

                g_xy_0_yzzz_z[k] = -g_xy_0_yzz_z[k] * ab_z + g_xy_0_yzz_zz[k];
            }

            /// Set up 87-90 components of targeted buffer : cbuffer.data(

            auto g_xy_0_zzzz_x = cbuffer.data(gp_geom_20_off + 87 * ccomps * dcomps);

            auto g_xy_0_zzzz_y = cbuffer.data(gp_geom_20_off + 88 * ccomps * dcomps);

            auto g_xy_0_zzzz_z = cbuffer.data(gp_geom_20_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_zzz_x, g_xy_0_zzz_xz, g_xy_0_zzz_y, g_xy_0_zzz_yz, g_xy_0_zzz_z, g_xy_0_zzz_zz, g_xy_0_zzzz_x, g_xy_0_zzzz_y, g_xy_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_zzzz_x[k] = -g_xy_0_zzz_x[k] * ab_z + g_xy_0_zzz_xz[k];

                g_xy_0_zzzz_y[k] = -g_xy_0_zzz_y[k] * ab_z + g_xy_0_zzz_yz[k];

                g_xy_0_zzzz_z[k] = -g_xy_0_zzz_z[k] * ab_z + g_xy_0_zzz_zz[k];
            }

            /// Set up 90-93 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxx_x = cbuffer.data(gp_geom_20_off + 90 * ccomps * dcomps);

            auto g_xz_0_xxxx_y = cbuffer.data(gp_geom_20_off + 91 * ccomps * dcomps);

            auto g_xz_0_xxxx_z = cbuffer.data(gp_geom_20_off + 92 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxx_x, g_xz_0_xxx_xx, g_xz_0_xxx_xy, g_xz_0_xxx_xz, g_xz_0_xxx_y, g_xz_0_xxx_z, g_xz_0_xxxx_x, g_xz_0_xxxx_y, g_xz_0_xxxx_z, g_z_0_xxx_x, g_z_0_xxx_y, g_z_0_xxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxx_x[k] = -g_z_0_xxx_x[k] - g_xz_0_xxx_x[k] * ab_x + g_xz_0_xxx_xx[k];

                g_xz_0_xxxx_y[k] = -g_z_0_xxx_y[k] - g_xz_0_xxx_y[k] * ab_x + g_xz_0_xxx_xy[k];

                g_xz_0_xxxx_z[k] = -g_z_0_xxx_z[k] - g_xz_0_xxx_z[k] * ab_x + g_xz_0_xxx_xz[k];
            }

            /// Set up 93-96 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxy_x = cbuffer.data(gp_geom_20_off + 93 * ccomps * dcomps);

            auto g_xz_0_xxxy_y = cbuffer.data(gp_geom_20_off + 94 * ccomps * dcomps);

            auto g_xz_0_xxxy_z = cbuffer.data(gp_geom_20_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxx_x, g_xz_0_xxx_xy, g_xz_0_xxx_y, g_xz_0_xxx_yy, g_xz_0_xxx_yz, g_xz_0_xxx_z, g_xz_0_xxxy_x, g_xz_0_xxxy_y, g_xz_0_xxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxy_x[k] = -g_xz_0_xxx_x[k] * ab_y + g_xz_0_xxx_xy[k];

                g_xz_0_xxxy_y[k] = -g_xz_0_xxx_y[k] * ab_y + g_xz_0_xxx_yy[k];

                g_xz_0_xxxy_z[k] = -g_xz_0_xxx_z[k] * ab_y + g_xz_0_xxx_yz[k];
            }

            /// Set up 96-99 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxz_x = cbuffer.data(gp_geom_20_off + 96 * ccomps * dcomps);

            auto g_xz_0_xxxz_y = cbuffer.data(gp_geom_20_off + 97 * ccomps * dcomps);

            auto g_xz_0_xxxz_z = cbuffer.data(gp_geom_20_off + 98 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxz_x, g_xz_0_xxxz_y, g_xz_0_xxxz_z, g_xz_0_xxz_x, g_xz_0_xxz_xx, g_xz_0_xxz_xy, g_xz_0_xxz_xz, g_xz_0_xxz_y, g_xz_0_xxz_z, g_z_0_xxz_x, g_z_0_xxz_y, g_z_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxz_x[k] = -g_z_0_xxz_x[k] - g_xz_0_xxz_x[k] * ab_x + g_xz_0_xxz_xx[k];

                g_xz_0_xxxz_y[k] = -g_z_0_xxz_y[k] - g_xz_0_xxz_y[k] * ab_x + g_xz_0_xxz_xy[k];

                g_xz_0_xxxz_z[k] = -g_z_0_xxz_z[k] - g_xz_0_xxz_z[k] * ab_x + g_xz_0_xxz_xz[k];
            }

            /// Set up 99-102 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyy_x = cbuffer.data(gp_geom_20_off + 99 * ccomps * dcomps);

            auto g_xz_0_xxyy_y = cbuffer.data(gp_geom_20_off + 100 * ccomps * dcomps);

            auto g_xz_0_xxyy_z = cbuffer.data(gp_geom_20_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxy_x, g_xz_0_xxy_xy, g_xz_0_xxy_y, g_xz_0_xxy_yy, g_xz_0_xxy_yz, g_xz_0_xxy_z, g_xz_0_xxyy_x, g_xz_0_xxyy_y, g_xz_0_xxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyy_x[k] = -g_xz_0_xxy_x[k] * ab_y + g_xz_0_xxy_xy[k];

                g_xz_0_xxyy_y[k] = -g_xz_0_xxy_y[k] * ab_y + g_xz_0_xxy_yy[k];

                g_xz_0_xxyy_z[k] = -g_xz_0_xxy_z[k] * ab_y + g_xz_0_xxy_yz[k];
            }

            /// Set up 102-105 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyz_x = cbuffer.data(gp_geom_20_off + 102 * ccomps * dcomps);

            auto g_xz_0_xxyz_y = cbuffer.data(gp_geom_20_off + 103 * ccomps * dcomps);

            auto g_xz_0_xxyz_z = cbuffer.data(gp_geom_20_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyz_x, g_xz_0_xxyz_y, g_xz_0_xxyz_z, g_xz_0_xxz_x, g_xz_0_xxz_xy, g_xz_0_xxz_y, g_xz_0_xxz_yy, g_xz_0_xxz_yz, g_xz_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyz_x[k] = -g_xz_0_xxz_x[k] * ab_y + g_xz_0_xxz_xy[k];

                g_xz_0_xxyz_y[k] = -g_xz_0_xxz_y[k] * ab_y + g_xz_0_xxz_yy[k];

                g_xz_0_xxyz_z[k] = -g_xz_0_xxz_z[k] * ab_y + g_xz_0_xxz_yz[k];
            }

            /// Set up 105-108 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxzz_x = cbuffer.data(gp_geom_20_off + 105 * ccomps * dcomps);

            auto g_xz_0_xxzz_y = cbuffer.data(gp_geom_20_off + 106 * ccomps * dcomps);

            auto g_xz_0_xxzz_z = cbuffer.data(gp_geom_20_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxzz_x, g_xz_0_xxzz_y, g_xz_0_xxzz_z, g_xz_0_xzz_x, g_xz_0_xzz_xx, g_xz_0_xzz_xy, g_xz_0_xzz_xz, g_xz_0_xzz_y, g_xz_0_xzz_z, g_z_0_xzz_x, g_z_0_xzz_y, g_z_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxzz_x[k] = -g_z_0_xzz_x[k] - g_xz_0_xzz_x[k] * ab_x + g_xz_0_xzz_xx[k];

                g_xz_0_xxzz_y[k] = -g_z_0_xzz_y[k] - g_xz_0_xzz_y[k] * ab_x + g_xz_0_xzz_xy[k];

                g_xz_0_xxzz_z[k] = -g_z_0_xzz_z[k] - g_xz_0_xzz_z[k] * ab_x + g_xz_0_xzz_xz[k];
            }

            /// Set up 108-111 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyy_x = cbuffer.data(gp_geom_20_off + 108 * ccomps * dcomps);

            auto g_xz_0_xyyy_y = cbuffer.data(gp_geom_20_off + 109 * ccomps * dcomps);

            auto g_xz_0_xyyy_z = cbuffer.data(gp_geom_20_off + 110 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyy_x, g_xz_0_xyy_xy, g_xz_0_xyy_y, g_xz_0_xyy_yy, g_xz_0_xyy_yz, g_xz_0_xyy_z, g_xz_0_xyyy_x, g_xz_0_xyyy_y, g_xz_0_xyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyy_x[k] = -g_xz_0_xyy_x[k] * ab_y + g_xz_0_xyy_xy[k];

                g_xz_0_xyyy_y[k] = -g_xz_0_xyy_y[k] * ab_y + g_xz_0_xyy_yy[k];

                g_xz_0_xyyy_z[k] = -g_xz_0_xyy_z[k] * ab_y + g_xz_0_xyy_yz[k];
            }

            /// Set up 111-114 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyz_x = cbuffer.data(gp_geom_20_off + 111 * ccomps * dcomps);

            auto g_xz_0_xyyz_y = cbuffer.data(gp_geom_20_off + 112 * ccomps * dcomps);

            auto g_xz_0_xyyz_z = cbuffer.data(gp_geom_20_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyz_x, g_xz_0_xyyz_y, g_xz_0_xyyz_z, g_xz_0_xyz_x, g_xz_0_xyz_xy, g_xz_0_xyz_y, g_xz_0_xyz_yy, g_xz_0_xyz_yz, g_xz_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyz_x[k] = -g_xz_0_xyz_x[k] * ab_y + g_xz_0_xyz_xy[k];

                g_xz_0_xyyz_y[k] = -g_xz_0_xyz_y[k] * ab_y + g_xz_0_xyz_yy[k];

                g_xz_0_xyyz_z[k] = -g_xz_0_xyz_z[k] * ab_y + g_xz_0_xyz_yz[k];
            }

            /// Set up 114-117 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyzz_x = cbuffer.data(gp_geom_20_off + 114 * ccomps * dcomps);

            auto g_xz_0_xyzz_y = cbuffer.data(gp_geom_20_off + 115 * ccomps * dcomps);

            auto g_xz_0_xyzz_z = cbuffer.data(gp_geom_20_off + 116 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyzz_x, g_xz_0_xyzz_y, g_xz_0_xyzz_z, g_xz_0_xzz_x, g_xz_0_xzz_xy, g_xz_0_xzz_y, g_xz_0_xzz_yy, g_xz_0_xzz_yz, g_xz_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyzz_x[k] = -g_xz_0_xzz_x[k] * ab_y + g_xz_0_xzz_xy[k];

                g_xz_0_xyzz_y[k] = -g_xz_0_xzz_y[k] * ab_y + g_xz_0_xzz_yy[k];

                g_xz_0_xyzz_z[k] = -g_xz_0_xzz_z[k] * ab_y + g_xz_0_xzz_yz[k];
            }

            /// Set up 117-120 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xzzz_x = cbuffer.data(gp_geom_20_off + 117 * ccomps * dcomps);

            auto g_xz_0_xzzz_y = cbuffer.data(gp_geom_20_off + 118 * ccomps * dcomps);

            auto g_xz_0_xzzz_z = cbuffer.data(gp_geom_20_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xzzz_x, g_xz_0_xzzz_y, g_xz_0_xzzz_z, g_xz_0_zzz_x, g_xz_0_zzz_xx, g_xz_0_zzz_xy, g_xz_0_zzz_xz, g_xz_0_zzz_y, g_xz_0_zzz_z, g_z_0_zzz_x, g_z_0_zzz_y, g_z_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xzzz_x[k] = -g_z_0_zzz_x[k] - g_xz_0_zzz_x[k] * ab_x + g_xz_0_zzz_xx[k];

                g_xz_0_xzzz_y[k] = -g_z_0_zzz_y[k] - g_xz_0_zzz_y[k] * ab_x + g_xz_0_zzz_xy[k];

                g_xz_0_xzzz_z[k] = -g_z_0_zzz_z[k] - g_xz_0_zzz_z[k] * ab_x + g_xz_0_zzz_xz[k];
            }

            /// Set up 120-123 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyy_x = cbuffer.data(gp_geom_20_off + 120 * ccomps * dcomps);

            auto g_xz_0_yyyy_y = cbuffer.data(gp_geom_20_off + 121 * ccomps * dcomps);

            auto g_xz_0_yyyy_z = cbuffer.data(gp_geom_20_off + 122 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyy_x, g_xz_0_yyy_xy, g_xz_0_yyy_y, g_xz_0_yyy_yy, g_xz_0_yyy_yz, g_xz_0_yyy_z, g_xz_0_yyyy_x, g_xz_0_yyyy_y, g_xz_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyy_x[k] = -g_xz_0_yyy_x[k] * ab_y + g_xz_0_yyy_xy[k];

                g_xz_0_yyyy_y[k] = -g_xz_0_yyy_y[k] * ab_y + g_xz_0_yyy_yy[k];

                g_xz_0_yyyy_z[k] = -g_xz_0_yyy_z[k] * ab_y + g_xz_0_yyy_yz[k];
            }

            /// Set up 123-126 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyz_x = cbuffer.data(gp_geom_20_off + 123 * ccomps * dcomps);

            auto g_xz_0_yyyz_y = cbuffer.data(gp_geom_20_off + 124 * ccomps * dcomps);

            auto g_xz_0_yyyz_z = cbuffer.data(gp_geom_20_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyz_x, g_xz_0_yyyz_y, g_xz_0_yyyz_z, g_xz_0_yyz_x, g_xz_0_yyz_xy, g_xz_0_yyz_y, g_xz_0_yyz_yy, g_xz_0_yyz_yz, g_xz_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyz_x[k] = -g_xz_0_yyz_x[k] * ab_y + g_xz_0_yyz_xy[k];

                g_xz_0_yyyz_y[k] = -g_xz_0_yyz_y[k] * ab_y + g_xz_0_yyz_yy[k];

                g_xz_0_yyyz_z[k] = -g_xz_0_yyz_z[k] * ab_y + g_xz_0_yyz_yz[k];
            }

            /// Set up 126-129 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyzz_x = cbuffer.data(gp_geom_20_off + 126 * ccomps * dcomps);

            auto g_xz_0_yyzz_y = cbuffer.data(gp_geom_20_off + 127 * ccomps * dcomps);

            auto g_xz_0_yyzz_z = cbuffer.data(gp_geom_20_off + 128 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyzz_x, g_xz_0_yyzz_y, g_xz_0_yyzz_z, g_xz_0_yzz_x, g_xz_0_yzz_xy, g_xz_0_yzz_y, g_xz_0_yzz_yy, g_xz_0_yzz_yz, g_xz_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyzz_x[k] = -g_xz_0_yzz_x[k] * ab_y + g_xz_0_yzz_xy[k];

                g_xz_0_yyzz_y[k] = -g_xz_0_yzz_y[k] * ab_y + g_xz_0_yzz_yy[k];

                g_xz_0_yyzz_z[k] = -g_xz_0_yzz_z[k] * ab_y + g_xz_0_yzz_yz[k];
            }

            /// Set up 129-132 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yzzz_x = cbuffer.data(gp_geom_20_off + 129 * ccomps * dcomps);

            auto g_xz_0_yzzz_y = cbuffer.data(gp_geom_20_off + 130 * ccomps * dcomps);

            auto g_xz_0_yzzz_z = cbuffer.data(gp_geom_20_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yzzz_x, g_xz_0_yzzz_y, g_xz_0_yzzz_z, g_xz_0_zzz_x, g_xz_0_zzz_xy, g_xz_0_zzz_y, g_xz_0_zzz_yy, g_xz_0_zzz_yz, g_xz_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yzzz_x[k] = -g_xz_0_zzz_x[k] * ab_y + g_xz_0_zzz_xy[k];

                g_xz_0_yzzz_y[k] = -g_xz_0_zzz_y[k] * ab_y + g_xz_0_zzz_yy[k];

                g_xz_0_yzzz_z[k] = -g_xz_0_zzz_z[k] * ab_y + g_xz_0_zzz_yz[k];
            }

            /// Set up 132-135 components of targeted buffer : cbuffer.data(

            auto g_xz_0_zzzz_x = cbuffer.data(gp_geom_20_off + 132 * ccomps * dcomps);

            auto g_xz_0_zzzz_y = cbuffer.data(gp_geom_20_off + 133 * ccomps * dcomps);

            auto g_xz_0_zzzz_z = cbuffer.data(gp_geom_20_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzz_x, g_x_0_zzz_y, g_x_0_zzz_z, g_xz_0_zzz_x, g_xz_0_zzz_xz, g_xz_0_zzz_y, g_xz_0_zzz_yz, g_xz_0_zzz_z, g_xz_0_zzz_zz, g_xz_0_zzzz_x, g_xz_0_zzzz_y, g_xz_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_zzzz_x[k] = -g_x_0_zzz_x[k] - g_xz_0_zzz_x[k] * ab_z + g_xz_0_zzz_xz[k];

                g_xz_0_zzzz_y[k] = -g_x_0_zzz_y[k] - g_xz_0_zzz_y[k] * ab_z + g_xz_0_zzz_yz[k];

                g_xz_0_zzzz_z[k] = -g_x_0_zzz_z[k] - g_xz_0_zzz_z[k] * ab_z + g_xz_0_zzz_zz[k];
            }

            /// Set up 135-138 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxx_x = cbuffer.data(gp_geom_20_off + 135 * ccomps * dcomps);

            auto g_yy_0_xxxx_y = cbuffer.data(gp_geom_20_off + 136 * ccomps * dcomps);

            auto g_yy_0_xxxx_z = cbuffer.data(gp_geom_20_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxx_x, g_yy_0_xxx_xx, g_yy_0_xxx_xy, g_yy_0_xxx_xz, g_yy_0_xxx_y, g_yy_0_xxx_z, g_yy_0_xxxx_x, g_yy_0_xxxx_y, g_yy_0_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxx_x[k] = -g_yy_0_xxx_x[k] * ab_x + g_yy_0_xxx_xx[k];

                g_yy_0_xxxx_y[k] = -g_yy_0_xxx_y[k] * ab_x + g_yy_0_xxx_xy[k];

                g_yy_0_xxxx_z[k] = -g_yy_0_xxx_z[k] * ab_x + g_yy_0_xxx_xz[k];
            }

            /// Set up 138-141 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxy_x = cbuffer.data(gp_geom_20_off + 138 * ccomps * dcomps);

            auto g_yy_0_xxxy_y = cbuffer.data(gp_geom_20_off + 139 * ccomps * dcomps);

            auto g_yy_0_xxxy_z = cbuffer.data(gp_geom_20_off + 140 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxy_x, g_yy_0_xxxy_y, g_yy_0_xxxy_z, g_yy_0_xxy_x, g_yy_0_xxy_xx, g_yy_0_xxy_xy, g_yy_0_xxy_xz, g_yy_0_xxy_y, g_yy_0_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxy_x[k] = -g_yy_0_xxy_x[k] * ab_x + g_yy_0_xxy_xx[k];

                g_yy_0_xxxy_y[k] = -g_yy_0_xxy_y[k] * ab_x + g_yy_0_xxy_xy[k];

                g_yy_0_xxxy_z[k] = -g_yy_0_xxy_z[k] * ab_x + g_yy_0_xxy_xz[k];
            }

            /// Set up 141-144 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxz_x = cbuffer.data(gp_geom_20_off + 141 * ccomps * dcomps);

            auto g_yy_0_xxxz_y = cbuffer.data(gp_geom_20_off + 142 * ccomps * dcomps);

            auto g_yy_0_xxxz_z = cbuffer.data(gp_geom_20_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxz_x, g_yy_0_xxxz_y, g_yy_0_xxxz_z, g_yy_0_xxz_x, g_yy_0_xxz_xx, g_yy_0_xxz_xy, g_yy_0_xxz_xz, g_yy_0_xxz_y, g_yy_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxz_x[k] = -g_yy_0_xxz_x[k] * ab_x + g_yy_0_xxz_xx[k];

                g_yy_0_xxxz_y[k] = -g_yy_0_xxz_y[k] * ab_x + g_yy_0_xxz_xy[k];

                g_yy_0_xxxz_z[k] = -g_yy_0_xxz_z[k] * ab_x + g_yy_0_xxz_xz[k];
            }

            /// Set up 144-147 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyy_x = cbuffer.data(gp_geom_20_off + 144 * ccomps * dcomps);

            auto g_yy_0_xxyy_y = cbuffer.data(gp_geom_20_off + 145 * ccomps * dcomps);

            auto g_yy_0_xxyy_z = cbuffer.data(gp_geom_20_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyy_x, g_yy_0_xxyy_y, g_yy_0_xxyy_z, g_yy_0_xyy_x, g_yy_0_xyy_xx, g_yy_0_xyy_xy, g_yy_0_xyy_xz, g_yy_0_xyy_y, g_yy_0_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyy_x[k] = -g_yy_0_xyy_x[k] * ab_x + g_yy_0_xyy_xx[k];

                g_yy_0_xxyy_y[k] = -g_yy_0_xyy_y[k] * ab_x + g_yy_0_xyy_xy[k];

                g_yy_0_xxyy_z[k] = -g_yy_0_xyy_z[k] * ab_x + g_yy_0_xyy_xz[k];
            }

            /// Set up 147-150 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyz_x = cbuffer.data(gp_geom_20_off + 147 * ccomps * dcomps);

            auto g_yy_0_xxyz_y = cbuffer.data(gp_geom_20_off + 148 * ccomps * dcomps);

            auto g_yy_0_xxyz_z = cbuffer.data(gp_geom_20_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyz_x, g_yy_0_xxyz_y, g_yy_0_xxyz_z, g_yy_0_xyz_x, g_yy_0_xyz_xx, g_yy_0_xyz_xy, g_yy_0_xyz_xz, g_yy_0_xyz_y, g_yy_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyz_x[k] = -g_yy_0_xyz_x[k] * ab_x + g_yy_0_xyz_xx[k];

                g_yy_0_xxyz_y[k] = -g_yy_0_xyz_y[k] * ab_x + g_yy_0_xyz_xy[k];

                g_yy_0_xxyz_z[k] = -g_yy_0_xyz_z[k] * ab_x + g_yy_0_xyz_xz[k];
            }

            /// Set up 150-153 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxzz_x = cbuffer.data(gp_geom_20_off + 150 * ccomps * dcomps);

            auto g_yy_0_xxzz_y = cbuffer.data(gp_geom_20_off + 151 * ccomps * dcomps);

            auto g_yy_0_xxzz_z = cbuffer.data(gp_geom_20_off + 152 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxzz_x, g_yy_0_xxzz_y, g_yy_0_xxzz_z, g_yy_0_xzz_x, g_yy_0_xzz_xx, g_yy_0_xzz_xy, g_yy_0_xzz_xz, g_yy_0_xzz_y, g_yy_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxzz_x[k] = -g_yy_0_xzz_x[k] * ab_x + g_yy_0_xzz_xx[k];

                g_yy_0_xxzz_y[k] = -g_yy_0_xzz_y[k] * ab_x + g_yy_0_xzz_xy[k];

                g_yy_0_xxzz_z[k] = -g_yy_0_xzz_z[k] * ab_x + g_yy_0_xzz_xz[k];
            }

            /// Set up 153-156 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyy_x = cbuffer.data(gp_geom_20_off + 153 * ccomps * dcomps);

            auto g_yy_0_xyyy_y = cbuffer.data(gp_geom_20_off + 154 * ccomps * dcomps);

            auto g_yy_0_xyyy_z = cbuffer.data(gp_geom_20_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyy_x, g_yy_0_xyyy_y, g_yy_0_xyyy_z, g_yy_0_yyy_x, g_yy_0_yyy_xx, g_yy_0_yyy_xy, g_yy_0_yyy_xz, g_yy_0_yyy_y, g_yy_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyy_x[k] = -g_yy_0_yyy_x[k] * ab_x + g_yy_0_yyy_xx[k];

                g_yy_0_xyyy_y[k] = -g_yy_0_yyy_y[k] * ab_x + g_yy_0_yyy_xy[k];

                g_yy_0_xyyy_z[k] = -g_yy_0_yyy_z[k] * ab_x + g_yy_0_yyy_xz[k];
            }

            /// Set up 156-159 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyz_x = cbuffer.data(gp_geom_20_off + 156 * ccomps * dcomps);

            auto g_yy_0_xyyz_y = cbuffer.data(gp_geom_20_off + 157 * ccomps * dcomps);

            auto g_yy_0_xyyz_z = cbuffer.data(gp_geom_20_off + 158 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyz_x, g_yy_0_xyyz_y, g_yy_0_xyyz_z, g_yy_0_yyz_x, g_yy_0_yyz_xx, g_yy_0_yyz_xy, g_yy_0_yyz_xz, g_yy_0_yyz_y, g_yy_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyz_x[k] = -g_yy_0_yyz_x[k] * ab_x + g_yy_0_yyz_xx[k];

                g_yy_0_xyyz_y[k] = -g_yy_0_yyz_y[k] * ab_x + g_yy_0_yyz_xy[k];

                g_yy_0_xyyz_z[k] = -g_yy_0_yyz_z[k] * ab_x + g_yy_0_yyz_xz[k];
            }

            /// Set up 159-162 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyzz_x = cbuffer.data(gp_geom_20_off + 159 * ccomps * dcomps);

            auto g_yy_0_xyzz_y = cbuffer.data(gp_geom_20_off + 160 * ccomps * dcomps);

            auto g_yy_0_xyzz_z = cbuffer.data(gp_geom_20_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyzz_x, g_yy_0_xyzz_y, g_yy_0_xyzz_z, g_yy_0_yzz_x, g_yy_0_yzz_xx, g_yy_0_yzz_xy, g_yy_0_yzz_xz, g_yy_0_yzz_y, g_yy_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyzz_x[k] = -g_yy_0_yzz_x[k] * ab_x + g_yy_0_yzz_xx[k];

                g_yy_0_xyzz_y[k] = -g_yy_0_yzz_y[k] * ab_x + g_yy_0_yzz_xy[k];

                g_yy_0_xyzz_z[k] = -g_yy_0_yzz_z[k] * ab_x + g_yy_0_yzz_xz[k];
            }

            /// Set up 162-165 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xzzz_x = cbuffer.data(gp_geom_20_off + 162 * ccomps * dcomps);

            auto g_yy_0_xzzz_y = cbuffer.data(gp_geom_20_off + 163 * ccomps * dcomps);

            auto g_yy_0_xzzz_z = cbuffer.data(gp_geom_20_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xzzz_x, g_yy_0_xzzz_y, g_yy_0_xzzz_z, g_yy_0_zzz_x, g_yy_0_zzz_xx, g_yy_0_zzz_xy, g_yy_0_zzz_xz, g_yy_0_zzz_y, g_yy_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xzzz_x[k] = -g_yy_0_zzz_x[k] * ab_x + g_yy_0_zzz_xx[k];

                g_yy_0_xzzz_y[k] = -g_yy_0_zzz_y[k] * ab_x + g_yy_0_zzz_xy[k];

                g_yy_0_xzzz_z[k] = -g_yy_0_zzz_z[k] * ab_x + g_yy_0_zzz_xz[k];
            }

            /// Set up 165-168 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyy_x = cbuffer.data(gp_geom_20_off + 165 * ccomps * dcomps);

            auto g_yy_0_yyyy_y = cbuffer.data(gp_geom_20_off + 166 * ccomps * dcomps);

            auto g_yy_0_yyyy_z = cbuffer.data(gp_geom_20_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyy_x, g_y_0_yyy_y, g_y_0_yyy_z, g_yy_0_yyy_x, g_yy_0_yyy_xy, g_yy_0_yyy_y, g_yy_0_yyy_yy, g_yy_0_yyy_yz, g_yy_0_yyy_z, g_yy_0_yyyy_x, g_yy_0_yyyy_y, g_yy_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyy_x[k] = -2.0 * g_y_0_yyy_x[k] - g_yy_0_yyy_x[k] * ab_y + g_yy_0_yyy_xy[k];

                g_yy_0_yyyy_y[k] = -2.0 * g_y_0_yyy_y[k] - g_yy_0_yyy_y[k] * ab_y + g_yy_0_yyy_yy[k];

                g_yy_0_yyyy_z[k] = -2.0 * g_y_0_yyy_z[k] - g_yy_0_yyy_z[k] * ab_y + g_yy_0_yyy_yz[k];
            }

            /// Set up 168-171 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyz_x = cbuffer.data(gp_geom_20_off + 168 * ccomps * dcomps);

            auto g_yy_0_yyyz_y = cbuffer.data(gp_geom_20_off + 169 * ccomps * dcomps);

            auto g_yy_0_yyyz_z = cbuffer.data(gp_geom_20_off + 170 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyy_x, g_yy_0_yyy_xz, g_yy_0_yyy_y, g_yy_0_yyy_yz, g_yy_0_yyy_z, g_yy_0_yyy_zz, g_yy_0_yyyz_x, g_yy_0_yyyz_y, g_yy_0_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyz_x[k] = -g_yy_0_yyy_x[k] * ab_z + g_yy_0_yyy_xz[k];

                g_yy_0_yyyz_y[k] = -g_yy_0_yyy_y[k] * ab_z + g_yy_0_yyy_yz[k];

                g_yy_0_yyyz_z[k] = -g_yy_0_yyy_z[k] * ab_z + g_yy_0_yyy_zz[k];
            }

            /// Set up 171-174 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyzz_x = cbuffer.data(gp_geom_20_off + 171 * ccomps * dcomps);

            auto g_yy_0_yyzz_y = cbuffer.data(gp_geom_20_off + 172 * ccomps * dcomps);

            auto g_yy_0_yyzz_z = cbuffer.data(gp_geom_20_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyz_x, g_yy_0_yyz_xz, g_yy_0_yyz_y, g_yy_0_yyz_yz, g_yy_0_yyz_z, g_yy_0_yyz_zz, g_yy_0_yyzz_x, g_yy_0_yyzz_y, g_yy_0_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyzz_x[k] = -g_yy_0_yyz_x[k] * ab_z + g_yy_0_yyz_xz[k];

                g_yy_0_yyzz_y[k] = -g_yy_0_yyz_y[k] * ab_z + g_yy_0_yyz_yz[k];

                g_yy_0_yyzz_z[k] = -g_yy_0_yyz_z[k] * ab_z + g_yy_0_yyz_zz[k];
            }

            /// Set up 174-177 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yzzz_x = cbuffer.data(gp_geom_20_off + 174 * ccomps * dcomps);

            auto g_yy_0_yzzz_y = cbuffer.data(gp_geom_20_off + 175 * ccomps * dcomps);

            auto g_yy_0_yzzz_z = cbuffer.data(gp_geom_20_off + 176 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yzz_x, g_yy_0_yzz_xz, g_yy_0_yzz_y, g_yy_0_yzz_yz, g_yy_0_yzz_z, g_yy_0_yzz_zz, g_yy_0_yzzz_x, g_yy_0_yzzz_y, g_yy_0_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yzzz_x[k] = -g_yy_0_yzz_x[k] * ab_z + g_yy_0_yzz_xz[k];

                g_yy_0_yzzz_y[k] = -g_yy_0_yzz_y[k] * ab_z + g_yy_0_yzz_yz[k];

                g_yy_0_yzzz_z[k] = -g_yy_0_yzz_z[k] * ab_z + g_yy_0_yzz_zz[k];
            }

            /// Set up 177-180 components of targeted buffer : cbuffer.data(

            auto g_yy_0_zzzz_x = cbuffer.data(gp_geom_20_off + 177 * ccomps * dcomps);

            auto g_yy_0_zzzz_y = cbuffer.data(gp_geom_20_off + 178 * ccomps * dcomps);

            auto g_yy_0_zzzz_z = cbuffer.data(gp_geom_20_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_zzz_x, g_yy_0_zzz_xz, g_yy_0_zzz_y, g_yy_0_zzz_yz, g_yy_0_zzz_z, g_yy_0_zzz_zz, g_yy_0_zzzz_x, g_yy_0_zzzz_y, g_yy_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_zzzz_x[k] = -g_yy_0_zzz_x[k] * ab_z + g_yy_0_zzz_xz[k];

                g_yy_0_zzzz_y[k] = -g_yy_0_zzz_y[k] * ab_z + g_yy_0_zzz_yz[k];

                g_yy_0_zzzz_z[k] = -g_yy_0_zzz_z[k] * ab_z + g_yy_0_zzz_zz[k];
            }

            /// Set up 180-183 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxx_x = cbuffer.data(gp_geom_20_off + 180 * ccomps * dcomps);

            auto g_yz_0_xxxx_y = cbuffer.data(gp_geom_20_off + 181 * ccomps * dcomps);

            auto g_yz_0_xxxx_z = cbuffer.data(gp_geom_20_off + 182 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxx_x, g_yz_0_xxx_xx, g_yz_0_xxx_xy, g_yz_0_xxx_xz, g_yz_0_xxx_y, g_yz_0_xxx_z, g_yz_0_xxxx_x, g_yz_0_xxxx_y, g_yz_0_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxx_x[k] = -g_yz_0_xxx_x[k] * ab_x + g_yz_0_xxx_xx[k];

                g_yz_0_xxxx_y[k] = -g_yz_0_xxx_y[k] * ab_x + g_yz_0_xxx_xy[k];

                g_yz_0_xxxx_z[k] = -g_yz_0_xxx_z[k] * ab_x + g_yz_0_xxx_xz[k];
            }

            /// Set up 183-186 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxy_x = cbuffer.data(gp_geom_20_off + 183 * ccomps * dcomps);

            auto g_yz_0_xxxy_y = cbuffer.data(gp_geom_20_off + 184 * ccomps * dcomps);

            auto g_yz_0_xxxy_z = cbuffer.data(gp_geom_20_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxy_x, g_yz_0_xxxy_y, g_yz_0_xxxy_z, g_yz_0_xxy_x, g_yz_0_xxy_xx, g_yz_0_xxy_xy, g_yz_0_xxy_xz, g_yz_0_xxy_y, g_yz_0_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxy_x[k] = -g_yz_0_xxy_x[k] * ab_x + g_yz_0_xxy_xx[k];

                g_yz_0_xxxy_y[k] = -g_yz_0_xxy_y[k] * ab_x + g_yz_0_xxy_xy[k];

                g_yz_0_xxxy_z[k] = -g_yz_0_xxy_z[k] * ab_x + g_yz_0_xxy_xz[k];
            }

            /// Set up 186-189 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxz_x = cbuffer.data(gp_geom_20_off + 186 * ccomps * dcomps);

            auto g_yz_0_xxxz_y = cbuffer.data(gp_geom_20_off + 187 * ccomps * dcomps);

            auto g_yz_0_xxxz_z = cbuffer.data(gp_geom_20_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxz_x, g_yz_0_xxxz_y, g_yz_0_xxxz_z, g_yz_0_xxz_x, g_yz_0_xxz_xx, g_yz_0_xxz_xy, g_yz_0_xxz_xz, g_yz_0_xxz_y, g_yz_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxz_x[k] = -g_yz_0_xxz_x[k] * ab_x + g_yz_0_xxz_xx[k];

                g_yz_0_xxxz_y[k] = -g_yz_0_xxz_y[k] * ab_x + g_yz_0_xxz_xy[k];

                g_yz_0_xxxz_z[k] = -g_yz_0_xxz_z[k] * ab_x + g_yz_0_xxz_xz[k];
            }

            /// Set up 189-192 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyy_x = cbuffer.data(gp_geom_20_off + 189 * ccomps * dcomps);

            auto g_yz_0_xxyy_y = cbuffer.data(gp_geom_20_off + 190 * ccomps * dcomps);

            auto g_yz_0_xxyy_z = cbuffer.data(gp_geom_20_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyy_x, g_yz_0_xxyy_y, g_yz_0_xxyy_z, g_yz_0_xyy_x, g_yz_0_xyy_xx, g_yz_0_xyy_xy, g_yz_0_xyy_xz, g_yz_0_xyy_y, g_yz_0_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyy_x[k] = -g_yz_0_xyy_x[k] * ab_x + g_yz_0_xyy_xx[k];

                g_yz_0_xxyy_y[k] = -g_yz_0_xyy_y[k] * ab_x + g_yz_0_xyy_xy[k];

                g_yz_0_xxyy_z[k] = -g_yz_0_xyy_z[k] * ab_x + g_yz_0_xyy_xz[k];
            }

            /// Set up 192-195 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyz_x = cbuffer.data(gp_geom_20_off + 192 * ccomps * dcomps);

            auto g_yz_0_xxyz_y = cbuffer.data(gp_geom_20_off + 193 * ccomps * dcomps);

            auto g_yz_0_xxyz_z = cbuffer.data(gp_geom_20_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyz_x, g_yz_0_xxyz_y, g_yz_0_xxyz_z, g_yz_0_xyz_x, g_yz_0_xyz_xx, g_yz_0_xyz_xy, g_yz_0_xyz_xz, g_yz_0_xyz_y, g_yz_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyz_x[k] = -g_yz_0_xyz_x[k] * ab_x + g_yz_0_xyz_xx[k];

                g_yz_0_xxyz_y[k] = -g_yz_0_xyz_y[k] * ab_x + g_yz_0_xyz_xy[k];

                g_yz_0_xxyz_z[k] = -g_yz_0_xyz_z[k] * ab_x + g_yz_0_xyz_xz[k];
            }

            /// Set up 195-198 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxzz_x = cbuffer.data(gp_geom_20_off + 195 * ccomps * dcomps);

            auto g_yz_0_xxzz_y = cbuffer.data(gp_geom_20_off + 196 * ccomps * dcomps);

            auto g_yz_0_xxzz_z = cbuffer.data(gp_geom_20_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxzz_x, g_yz_0_xxzz_y, g_yz_0_xxzz_z, g_yz_0_xzz_x, g_yz_0_xzz_xx, g_yz_0_xzz_xy, g_yz_0_xzz_xz, g_yz_0_xzz_y, g_yz_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxzz_x[k] = -g_yz_0_xzz_x[k] * ab_x + g_yz_0_xzz_xx[k];

                g_yz_0_xxzz_y[k] = -g_yz_0_xzz_y[k] * ab_x + g_yz_0_xzz_xy[k];

                g_yz_0_xxzz_z[k] = -g_yz_0_xzz_z[k] * ab_x + g_yz_0_xzz_xz[k];
            }

            /// Set up 198-201 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyy_x = cbuffer.data(gp_geom_20_off + 198 * ccomps * dcomps);

            auto g_yz_0_xyyy_y = cbuffer.data(gp_geom_20_off + 199 * ccomps * dcomps);

            auto g_yz_0_xyyy_z = cbuffer.data(gp_geom_20_off + 200 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyy_x, g_yz_0_xyyy_y, g_yz_0_xyyy_z, g_yz_0_yyy_x, g_yz_0_yyy_xx, g_yz_0_yyy_xy, g_yz_0_yyy_xz, g_yz_0_yyy_y, g_yz_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyy_x[k] = -g_yz_0_yyy_x[k] * ab_x + g_yz_0_yyy_xx[k];

                g_yz_0_xyyy_y[k] = -g_yz_0_yyy_y[k] * ab_x + g_yz_0_yyy_xy[k];

                g_yz_0_xyyy_z[k] = -g_yz_0_yyy_z[k] * ab_x + g_yz_0_yyy_xz[k];
            }

            /// Set up 201-204 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyz_x = cbuffer.data(gp_geom_20_off + 201 * ccomps * dcomps);

            auto g_yz_0_xyyz_y = cbuffer.data(gp_geom_20_off + 202 * ccomps * dcomps);

            auto g_yz_0_xyyz_z = cbuffer.data(gp_geom_20_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyz_x, g_yz_0_xyyz_y, g_yz_0_xyyz_z, g_yz_0_yyz_x, g_yz_0_yyz_xx, g_yz_0_yyz_xy, g_yz_0_yyz_xz, g_yz_0_yyz_y, g_yz_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyz_x[k] = -g_yz_0_yyz_x[k] * ab_x + g_yz_0_yyz_xx[k];

                g_yz_0_xyyz_y[k] = -g_yz_0_yyz_y[k] * ab_x + g_yz_0_yyz_xy[k];

                g_yz_0_xyyz_z[k] = -g_yz_0_yyz_z[k] * ab_x + g_yz_0_yyz_xz[k];
            }

            /// Set up 204-207 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyzz_x = cbuffer.data(gp_geom_20_off + 204 * ccomps * dcomps);

            auto g_yz_0_xyzz_y = cbuffer.data(gp_geom_20_off + 205 * ccomps * dcomps);

            auto g_yz_0_xyzz_z = cbuffer.data(gp_geom_20_off + 206 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyzz_x, g_yz_0_xyzz_y, g_yz_0_xyzz_z, g_yz_0_yzz_x, g_yz_0_yzz_xx, g_yz_0_yzz_xy, g_yz_0_yzz_xz, g_yz_0_yzz_y, g_yz_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyzz_x[k] = -g_yz_0_yzz_x[k] * ab_x + g_yz_0_yzz_xx[k];

                g_yz_0_xyzz_y[k] = -g_yz_0_yzz_y[k] * ab_x + g_yz_0_yzz_xy[k];

                g_yz_0_xyzz_z[k] = -g_yz_0_yzz_z[k] * ab_x + g_yz_0_yzz_xz[k];
            }

            /// Set up 207-210 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xzzz_x = cbuffer.data(gp_geom_20_off + 207 * ccomps * dcomps);

            auto g_yz_0_xzzz_y = cbuffer.data(gp_geom_20_off + 208 * ccomps * dcomps);

            auto g_yz_0_xzzz_z = cbuffer.data(gp_geom_20_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xzzz_x, g_yz_0_xzzz_y, g_yz_0_xzzz_z, g_yz_0_zzz_x, g_yz_0_zzz_xx, g_yz_0_zzz_xy, g_yz_0_zzz_xz, g_yz_0_zzz_y, g_yz_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xzzz_x[k] = -g_yz_0_zzz_x[k] * ab_x + g_yz_0_zzz_xx[k];

                g_yz_0_xzzz_y[k] = -g_yz_0_zzz_y[k] * ab_x + g_yz_0_zzz_xy[k];

                g_yz_0_xzzz_z[k] = -g_yz_0_zzz_z[k] * ab_x + g_yz_0_zzz_xz[k];
            }

            /// Set up 210-213 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyy_x = cbuffer.data(gp_geom_20_off + 210 * ccomps * dcomps);

            auto g_yz_0_yyyy_y = cbuffer.data(gp_geom_20_off + 211 * ccomps * dcomps);

            auto g_yz_0_yyyy_z = cbuffer.data(gp_geom_20_off + 212 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyy_x, g_yz_0_yyy_xy, g_yz_0_yyy_y, g_yz_0_yyy_yy, g_yz_0_yyy_yz, g_yz_0_yyy_z, g_yz_0_yyyy_x, g_yz_0_yyyy_y, g_yz_0_yyyy_z, g_z_0_yyy_x, g_z_0_yyy_y, g_z_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyy_x[k] = -g_z_0_yyy_x[k] - g_yz_0_yyy_x[k] * ab_y + g_yz_0_yyy_xy[k];

                g_yz_0_yyyy_y[k] = -g_z_0_yyy_y[k] - g_yz_0_yyy_y[k] * ab_y + g_yz_0_yyy_yy[k];

                g_yz_0_yyyy_z[k] = -g_z_0_yyy_z[k] - g_yz_0_yyy_z[k] * ab_y + g_yz_0_yyy_yz[k];
            }

            /// Set up 213-216 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyz_x = cbuffer.data(gp_geom_20_off + 213 * ccomps * dcomps);

            auto g_yz_0_yyyz_y = cbuffer.data(gp_geom_20_off + 214 * ccomps * dcomps);

            auto g_yz_0_yyyz_z = cbuffer.data(gp_geom_20_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyz_x, g_yz_0_yyyz_y, g_yz_0_yyyz_z, g_yz_0_yyz_x, g_yz_0_yyz_xy, g_yz_0_yyz_y, g_yz_0_yyz_yy, g_yz_0_yyz_yz, g_yz_0_yyz_z, g_z_0_yyz_x, g_z_0_yyz_y, g_z_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyz_x[k] = -g_z_0_yyz_x[k] - g_yz_0_yyz_x[k] * ab_y + g_yz_0_yyz_xy[k];

                g_yz_0_yyyz_y[k] = -g_z_0_yyz_y[k] - g_yz_0_yyz_y[k] * ab_y + g_yz_0_yyz_yy[k];

                g_yz_0_yyyz_z[k] = -g_z_0_yyz_z[k] - g_yz_0_yyz_z[k] * ab_y + g_yz_0_yyz_yz[k];
            }

            /// Set up 216-219 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyzz_x = cbuffer.data(gp_geom_20_off + 216 * ccomps * dcomps);

            auto g_yz_0_yyzz_y = cbuffer.data(gp_geom_20_off + 217 * ccomps * dcomps);

            auto g_yz_0_yyzz_z = cbuffer.data(gp_geom_20_off + 218 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyzz_x, g_yz_0_yyzz_y, g_yz_0_yyzz_z, g_yz_0_yzz_x, g_yz_0_yzz_xy, g_yz_0_yzz_y, g_yz_0_yzz_yy, g_yz_0_yzz_yz, g_yz_0_yzz_z, g_z_0_yzz_x, g_z_0_yzz_y, g_z_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyzz_x[k] = -g_z_0_yzz_x[k] - g_yz_0_yzz_x[k] * ab_y + g_yz_0_yzz_xy[k];

                g_yz_0_yyzz_y[k] = -g_z_0_yzz_y[k] - g_yz_0_yzz_y[k] * ab_y + g_yz_0_yzz_yy[k];

                g_yz_0_yyzz_z[k] = -g_z_0_yzz_z[k] - g_yz_0_yzz_z[k] * ab_y + g_yz_0_yzz_yz[k];
            }

            /// Set up 219-222 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yzzz_x = cbuffer.data(gp_geom_20_off + 219 * ccomps * dcomps);

            auto g_yz_0_yzzz_y = cbuffer.data(gp_geom_20_off + 220 * ccomps * dcomps);

            auto g_yz_0_yzzz_z = cbuffer.data(gp_geom_20_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yzzz_x, g_yz_0_yzzz_y, g_yz_0_yzzz_z, g_yz_0_zzz_x, g_yz_0_zzz_xy, g_yz_0_zzz_y, g_yz_0_zzz_yy, g_yz_0_zzz_yz, g_yz_0_zzz_z, g_z_0_zzz_x, g_z_0_zzz_y, g_z_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yzzz_x[k] = -g_z_0_zzz_x[k] - g_yz_0_zzz_x[k] * ab_y + g_yz_0_zzz_xy[k];

                g_yz_0_yzzz_y[k] = -g_z_0_zzz_y[k] - g_yz_0_zzz_y[k] * ab_y + g_yz_0_zzz_yy[k];

                g_yz_0_yzzz_z[k] = -g_z_0_zzz_z[k] - g_yz_0_zzz_z[k] * ab_y + g_yz_0_zzz_yz[k];
            }

            /// Set up 222-225 components of targeted buffer : cbuffer.data(

            auto g_yz_0_zzzz_x = cbuffer.data(gp_geom_20_off + 222 * ccomps * dcomps);

            auto g_yz_0_zzzz_y = cbuffer.data(gp_geom_20_off + 223 * ccomps * dcomps);

            auto g_yz_0_zzzz_z = cbuffer.data(gp_geom_20_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzz_x, g_y_0_zzz_y, g_y_0_zzz_z, g_yz_0_zzz_x, g_yz_0_zzz_xz, g_yz_0_zzz_y, g_yz_0_zzz_yz, g_yz_0_zzz_z, g_yz_0_zzz_zz, g_yz_0_zzzz_x, g_yz_0_zzzz_y, g_yz_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_zzzz_x[k] = -g_y_0_zzz_x[k] - g_yz_0_zzz_x[k] * ab_z + g_yz_0_zzz_xz[k];

                g_yz_0_zzzz_y[k] = -g_y_0_zzz_y[k] - g_yz_0_zzz_y[k] * ab_z + g_yz_0_zzz_yz[k];

                g_yz_0_zzzz_z[k] = -g_y_0_zzz_z[k] - g_yz_0_zzz_z[k] * ab_z + g_yz_0_zzz_zz[k];
            }

            /// Set up 225-228 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxx_x = cbuffer.data(gp_geom_20_off + 225 * ccomps * dcomps);

            auto g_zz_0_xxxx_y = cbuffer.data(gp_geom_20_off + 226 * ccomps * dcomps);

            auto g_zz_0_xxxx_z = cbuffer.data(gp_geom_20_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxx_x, g_zz_0_xxx_xx, g_zz_0_xxx_xy, g_zz_0_xxx_xz, g_zz_0_xxx_y, g_zz_0_xxx_z, g_zz_0_xxxx_x, g_zz_0_xxxx_y, g_zz_0_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxx_x[k] = -g_zz_0_xxx_x[k] * ab_x + g_zz_0_xxx_xx[k];

                g_zz_0_xxxx_y[k] = -g_zz_0_xxx_y[k] * ab_x + g_zz_0_xxx_xy[k];

                g_zz_0_xxxx_z[k] = -g_zz_0_xxx_z[k] * ab_x + g_zz_0_xxx_xz[k];
            }

            /// Set up 228-231 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxy_x = cbuffer.data(gp_geom_20_off + 228 * ccomps * dcomps);

            auto g_zz_0_xxxy_y = cbuffer.data(gp_geom_20_off + 229 * ccomps * dcomps);

            auto g_zz_0_xxxy_z = cbuffer.data(gp_geom_20_off + 230 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxy_x, g_zz_0_xxxy_y, g_zz_0_xxxy_z, g_zz_0_xxy_x, g_zz_0_xxy_xx, g_zz_0_xxy_xy, g_zz_0_xxy_xz, g_zz_0_xxy_y, g_zz_0_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxy_x[k] = -g_zz_0_xxy_x[k] * ab_x + g_zz_0_xxy_xx[k];

                g_zz_0_xxxy_y[k] = -g_zz_0_xxy_y[k] * ab_x + g_zz_0_xxy_xy[k];

                g_zz_0_xxxy_z[k] = -g_zz_0_xxy_z[k] * ab_x + g_zz_0_xxy_xz[k];
            }

            /// Set up 231-234 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxz_x = cbuffer.data(gp_geom_20_off + 231 * ccomps * dcomps);

            auto g_zz_0_xxxz_y = cbuffer.data(gp_geom_20_off + 232 * ccomps * dcomps);

            auto g_zz_0_xxxz_z = cbuffer.data(gp_geom_20_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxz_x, g_zz_0_xxxz_y, g_zz_0_xxxz_z, g_zz_0_xxz_x, g_zz_0_xxz_xx, g_zz_0_xxz_xy, g_zz_0_xxz_xz, g_zz_0_xxz_y, g_zz_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxz_x[k] = -g_zz_0_xxz_x[k] * ab_x + g_zz_0_xxz_xx[k];

                g_zz_0_xxxz_y[k] = -g_zz_0_xxz_y[k] * ab_x + g_zz_0_xxz_xy[k];

                g_zz_0_xxxz_z[k] = -g_zz_0_xxz_z[k] * ab_x + g_zz_0_xxz_xz[k];
            }

            /// Set up 234-237 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyy_x = cbuffer.data(gp_geom_20_off + 234 * ccomps * dcomps);

            auto g_zz_0_xxyy_y = cbuffer.data(gp_geom_20_off + 235 * ccomps * dcomps);

            auto g_zz_0_xxyy_z = cbuffer.data(gp_geom_20_off + 236 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyy_x, g_zz_0_xxyy_y, g_zz_0_xxyy_z, g_zz_0_xyy_x, g_zz_0_xyy_xx, g_zz_0_xyy_xy, g_zz_0_xyy_xz, g_zz_0_xyy_y, g_zz_0_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyy_x[k] = -g_zz_0_xyy_x[k] * ab_x + g_zz_0_xyy_xx[k];

                g_zz_0_xxyy_y[k] = -g_zz_0_xyy_y[k] * ab_x + g_zz_0_xyy_xy[k];

                g_zz_0_xxyy_z[k] = -g_zz_0_xyy_z[k] * ab_x + g_zz_0_xyy_xz[k];
            }

            /// Set up 237-240 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyz_x = cbuffer.data(gp_geom_20_off + 237 * ccomps * dcomps);

            auto g_zz_0_xxyz_y = cbuffer.data(gp_geom_20_off + 238 * ccomps * dcomps);

            auto g_zz_0_xxyz_z = cbuffer.data(gp_geom_20_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyz_x, g_zz_0_xxyz_y, g_zz_0_xxyz_z, g_zz_0_xyz_x, g_zz_0_xyz_xx, g_zz_0_xyz_xy, g_zz_0_xyz_xz, g_zz_0_xyz_y, g_zz_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyz_x[k] = -g_zz_0_xyz_x[k] * ab_x + g_zz_0_xyz_xx[k];

                g_zz_0_xxyz_y[k] = -g_zz_0_xyz_y[k] * ab_x + g_zz_0_xyz_xy[k];

                g_zz_0_xxyz_z[k] = -g_zz_0_xyz_z[k] * ab_x + g_zz_0_xyz_xz[k];
            }

            /// Set up 240-243 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxzz_x = cbuffer.data(gp_geom_20_off + 240 * ccomps * dcomps);

            auto g_zz_0_xxzz_y = cbuffer.data(gp_geom_20_off + 241 * ccomps * dcomps);

            auto g_zz_0_xxzz_z = cbuffer.data(gp_geom_20_off + 242 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxzz_x, g_zz_0_xxzz_y, g_zz_0_xxzz_z, g_zz_0_xzz_x, g_zz_0_xzz_xx, g_zz_0_xzz_xy, g_zz_0_xzz_xz, g_zz_0_xzz_y, g_zz_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxzz_x[k] = -g_zz_0_xzz_x[k] * ab_x + g_zz_0_xzz_xx[k];

                g_zz_0_xxzz_y[k] = -g_zz_0_xzz_y[k] * ab_x + g_zz_0_xzz_xy[k];

                g_zz_0_xxzz_z[k] = -g_zz_0_xzz_z[k] * ab_x + g_zz_0_xzz_xz[k];
            }

            /// Set up 243-246 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyy_x = cbuffer.data(gp_geom_20_off + 243 * ccomps * dcomps);

            auto g_zz_0_xyyy_y = cbuffer.data(gp_geom_20_off + 244 * ccomps * dcomps);

            auto g_zz_0_xyyy_z = cbuffer.data(gp_geom_20_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyy_x, g_zz_0_xyyy_y, g_zz_0_xyyy_z, g_zz_0_yyy_x, g_zz_0_yyy_xx, g_zz_0_yyy_xy, g_zz_0_yyy_xz, g_zz_0_yyy_y, g_zz_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyy_x[k] = -g_zz_0_yyy_x[k] * ab_x + g_zz_0_yyy_xx[k];

                g_zz_0_xyyy_y[k] = -g_zz_0_yyy_y[k] * ab_x + g_zz_0_yyy_xy[k];

                g_zz_0_xyyy_z[k] = -g_zz_0_yyy_z[k] * ab_x + g_zz_0_yyy_xz[k];
            }

            /// Set up 246-249 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyz_x = cbuffer.data(gp_geom_20_off + 246 * ccomps * dcomps);

            auto g_zz_0_xyyz_y = cbuffer.data(gp_geom_20_off + 247 * ccomps * dcomps);

            auto g_zz_0_xyyz_z = cbuffer.data(gp_geom_20_off + 248 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyz_x, g_zz_0_xyyz_y, g_zz_0_xyyz_z, g_zz_0_yyz_x, g_zz_0_yyz_xx, g_zz_0_yyz_xy, g_zz_0_yyz_xz, g_zz_0_yyz_y, g_zz_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyz_x[k] = -g_zz_0_yyz_x[k] * ab_x + g_zz_0_yyz_xx[k];

                g_zz_0_xyyz_y[k] = -g_zz_0_yyz_y[k] * ab_x + g_zz_0_yyz_xy[k];

                g_zz_0_xyyz_z[k] = -g_zz_0_yyz_z[k] * ab_x + g_zz_0_yyz_xz[k];
            }

            /// Set up 249-252 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyzz_x = cbuffer.data(gp_geom_20_off + 249 * ccomps * dcomps);

            auto g_zz_0_xyzz_y = cbuffer.data(gp_geom_20_off + 250 * ccomps * dcomps);

            auto g_zz_0_xyzz_z = cbuffer.data(gp_geom_20_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyzz_x, g_zz_0_xyzz_y, g_zz_0_xyzz_z, g_zz_0_yzz_x, g_zz_0_yzz_xx, g_zz_0_yzz_xy, g_zz_0_yzz_xz, g_zz_0_yzz_y, g_zz_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyzz_x[k] = -g_zz_0_yzz_x[k] * ab_x + g_zz_0_yzz_xx[k];

                g_zz_0_xyzz_y[k] = -g_zz_0_yzz_y[k] * ab_x + g_zz_0_yzz_xy[k];

                g_zz_0_xyzz_z[k] = -g_zz_0_yzz_z[k] * ab_x + g_zz_0_yzz_xz[k];
            }

            /// Set up 252-255 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xzzz_x = cbuffer.data(gp_geom_20_off + 252 * ccomps * dcomps);

            auto g_zz_0_xzzz_y = cbuffer.data(gp_geom_20_off + 253 * ccomps * dcomps);

            auto g_zz_0_xzzz_z = cbuffer.data(gp_geom_20_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xzzz_x, g_zz_0_xzzz_y, g_zz_0_xzzz_z, g_zz_0_zzz_x, g_zz_0_zzz_xx, g_zz_0_zzz_xy, g_zz_0_zzz_xz, g_zz_0_zzz_y, g_zz_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xzzz_x[k] = -g_zz_0_zzz_x[k] * ab_x + g_zz_0_zzz_xx[k];

                g_zz_0_xzzz_y[k] = -g_zz_0_zzz_y[k] * ab_x + g_zz_0_zzz_xy[k];

                g_zz_0_xzzz_z[k] = -g_zz_0_zzz_z[k] * ab_x + g_zz_0_zzz_xz[k];
            }

            /// Set up 255-258 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyy_x = cbuffer.data(gp_geom_20_off + 255 * ccomps * dcomps);

            auto g_zz_0_yyyy_y = cbuffer.data(gp_geom_20_off + 256 * ccomps * dcomps);

            auto g_zz_0_yyyy_z = cbuffer.data(gp_geom_20_off + 257 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyy_x, g_zz_0_yyy_xy, g_zz_0_yyy_y, g_zz_0_yyy_yy, g_zz_0_yyy_yz, g_zz_0_yyy_z, g_zz_0_yyyy_x, g_zz_0_yyyy_y, g_zz_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyy_x[k] = -g_zz_0_yyy_x[k] * ab_y + g_zz_0_yyy_xy[k];

                g_zz_0_yyyy_y[k] = -g_zz_0_yyy_y[k] * ab_y + g_zz_0_yyy_yy[k];

                g_zz_0_yyyy_z[k] = -g_zz_0_yyy_z[k] * ab_y + g_zz_0_yyy_yz[k];
            }

            /// Set up 258-261 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyz_x = cbuffer.data(gp_geom_20_off + 258 * ccomps * dcomps);

            auto g_zz_0_yyyz_y = cbuffer.data(gp_geom_20_off + 259 * ccomps * dcomps);

            auto g_zz_0_yyyz_z = cbuffer.data(gp_geom_20_off + 260 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyz_x, g_zz_0_yyyz_y, g_zz_0_yyyz_z, g_zz_0_yyz_x, g_zz_0_yyz_xy, g_zz_0_yyz_y, g_zz_0_yyz_yy, g_zz_0_yyz_yz, g_zz_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyz_x[k] = -g_zz_0_yyz_x[k] * ab_y + g_zz_0_yyz_xy[k];

                g_zz_0_yyyz_y[k] = -g_zz_0_yyz_y[k] * ab_y + g_zz_0_yyz_yy[k];

                g_zz_0_yyyz_z[k] = -g_zz_0_yyz_z[k] * ab_y + g_zz_0_yyz_yz[k];
            }

            /// Set up 261-264 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyzz_x = cbuffer.data(gp_geom_20_off + 261 * ccomps * dcomps);

            auto g_zz_0_yyzz_y = cbuffer.data(gp_geom_20_off + 262 * ccomps * dcomps);

            auto g_zz_0_yyzz_z = cbuffer.data(gp_geom_20_off + 263 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyzz_x, g_zz_0_yyzz_y, g_zz_0_yyzz_z, g_zz_0_yzz_x, g_zz_0_yzz_xy, g_zz_0_yzz_y, g_zz_0_yzz_yy, g_zz_0_yzz_yz, g_zz_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyzz_x[k] = -g_zz_0_yzz_x[k] * ab_y + g_zz_0_yzz_xy[k];

                g_zz_0_yyzz_y[k] = -g_zz_0_yzz_y[k] * ab_y + g_zz_0_yzz_yy[k];

                g_zz_0_yyzz_z[k] = -g_zz_0_yzz_z[k] * ab_y + g_zz_0_yzz_yz[k];
            }

            /// Set up 264-267 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yzzz_x = cbuffer.data(gp_geom_20_off + 264 * ccomps * dcomps);

            auto g_zz_0_yzzz_y = cbuffer.data(gp_geom_20_off + 265 * ccomps * dcomps);

            auto g_zz_0_yzzz_z = cbuffer.data(gp_geom_20_off + 266 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yzzz_x, g_zz_0_yzzz_y, g_zz_0_yzzz_z, g_zz_0_zzz_x, g_zz_0_zzz_xy, g_zz_0_zzz_y, g_zz_0_zzz_yy, g_zz_0_zzz_yz, g_zz_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yzzz_x[k] = -g_zz_0_zzz_x[k] * ab_y + g_zz_0_zzz_xy[k];

                g_zz_0_yzzz_y[k] = -g_zz_0_zzz_y[k] * ab_y + g_zz_0_zzz_yy[k];

                g_zz_0_yzzz_z[k] = -g_zz_0_zzz_z[k] * ab_y + g_zz_0_zzz_yz[k];
            }

            /// Set up 267-270 components of targeted buffer : cbuffer.data(

            auto g_zz_0_zzzz_x = cbuffer.data(gp_geom_20_off + 267 * ccomps * dcomps);

            auto g_zz_0_zzzz_y = cbuffer.data(gp_geom_20_off + 268 * ccomps * dcomps);

            auto g_zz_0_zzzz_z = cbuffer.data(gp_geom_20_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzz_x, g_z_0_zzz_y, g_z_0_zzz_z, g_zz_0_zzz_x, g_zz_0_zzz_xz, g_zz_0_zzz_y, g_zz_0_zzz_yz, g_zz_0_zzz_z, g_zz_0_zzz_zz, g_zz_0_zzzz_x, g_zz_0_zzzz_y, g_zz_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_zzzz_x[k] = -2.0 * g_z_0_zzz_x[k] - g_zz_0_zzz_x[k] * ab_z + g_zz_0_zzz_xz[k];

                g_zz_0_zzzz_y[k] = -2.0 * g_z_0_zzz_y[k] - g_zz_0_zzz_y[k] * ab_z + g_zz_0_zzz_yz[k];

                g_zz_0_zzzz_z[k] = -2.0 * g_z_0_zzz_z[k] - g_zz_0_zzz_z[k] * ab_z + g_zz_0_zzz_zz[k];
            }
        }
    }
}

} // erirec namespace

