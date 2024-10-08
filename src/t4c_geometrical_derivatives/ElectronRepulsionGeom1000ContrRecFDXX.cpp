#include "ElectronRepulsionGeom1000ContrRecFDXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_fdxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_fdxx,
                                            const size_t idx_ddxx,
                                            const size_t idx_geom_10_ddxx,
                                            const size_t idx_geom_10_dfxx,
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
            /// Set up components of auxilary buffer : DDSS

            const auto dd_off = idx_ddxx + i * dcomps + j;

            auto g_xx_xx = cbuffer.data(dd_off + 0 * ccomps * dcomps);

            auto g_xx_xy = cbuffer.data(dd_off + 1 * ccomps * dcomps);

            auto g_xx_xz = cbuffer.data(dd_off + 2 * ccomps * dcomps);

            auto g_xx_yy = cbuffer.data(dd_off + 3 * ccomps * dcomps);

            auto g_xx_yz = cbuffer.data(dd_off + 4 * ccomps * dcomps);

            auto g_xx_zz = cbuffer.data(dd_off + 5 * ccomps * dcomps);

            auto g_xy_xx = cbuffer.data(dd_off + 6 * ccomps * dcomps);

            auto g_xy_xy = cbuffer.data(dd_off + 7 * ccomps * dcomps);

            auto g_xy_xz = cbuffer.data(dd_off + 8 * ccomps * dcomps);

            auto g_xy_yy = cbuffer.data(dd_off + 9 * ccomps * dcomps);

            auto g_xy_yz = cbuffer.data(dd_off + 10 * ccomps * dcomps);

            auto g_xy_zz = cbuffer.data(dd_off + 11 * ccomps * dcomps);

            auto g_xz_xx = cbuffer.data(dd_off + 12 * ccomps * dcomps);

            auto g_xz_xy = cbuffer.data(dd_off + 13 * ccomps * dcomps);

            auto g_xz_xz = cbuffer.data(dd_off + 14 * ccomps * dcomps);

            auto g_xz_yy = cbuffer.data(dd_off + 15 * ccomps * dcomps);

            auto g_xz_yz = cbuffer.data(dd_off + 16 * ccomps * dcomps);

            auto g_xz_zz = cbuffer.data(dd_off + 17 * ccomps * dcomps);

            auto g_yy_xx = cbuffer.data(dd_off + 18 * ccomps * dcomps);

            auto g_yy_xy = cbuffer.data(dd_off + 19 * ccomps * dcomps);

            auto g_yy_xz = cbuffer.data(dd_off + 20 * ccomps * dcomps);

            auto g_yy_yy = cbuffer.data(dd_off + 21 * ccomps * dcomps);

            auto g_yy_yz = cbuffer.data(dd_off + 22 * ccomps * dcomps);

            auto g_yy_zz = cbuffer.data(dd_off + 23 * ccomps * dcomps);

            auto g_yz_xx = cbuffer.data(dd_off + 24 * ccomps * dcomps);

            auto g_yz_xy = cbuffer.data(dd_off + 25 * ccomps * dcomps);

            auto g_yz_xz = cbuffer.data(dd_off + 26 * ccomps * dcomps);

            auto g_yz_yy = cbuffer.data(dd_off + 27 * ccomps * dcomps);

            auto g_yz_yz = cbuffer.data(dd_off + 28 * ccomps * dcomps);

            auto g_yz_zz = cbuffer.data(dd_off + 29 * ccomps * dcomps);

            auto g_zz_xx = cbuffer.data(dd_off + 30 * ccomps * dcomps);

            auto g_zz_xy = cbuffer.data(dd_off + 31 * ccomps * dcomps);

            auto g_zz_xz = cbuffer.data(dd_off + 32 * ccomps * dcomps);

            auto g_zz_yy = cbuffer.data(dd_off + 33 * ccomps * dcomps);

            auto g_zz_yz = cbuffer.data(dd_off + 34 * ccomps * dcomps);

            auto g_zz_zz = cbuffer.data(dd_off + 35 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DDSS

            const auto dd_geom_10_off = idx_geom_10_ddxx + i * dcomps + j;

            auto g_x_0_xx_xx = cbuffer.data(dd_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xx_xy = cbuffer.data(dd_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xx_xz = cbuffer.data(dd_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xx_yy = cbuffer.data(dd_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xx_yz = cbuffer.data(dd_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xx_zz = cbuffer.data(dd_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xy_xx = cbuffer.data(dd_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xy_xy = cbuffer.data(dd_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xy_xz = cbuffer.data(dd_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xy_yy = cbuffer.data(dd_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xy_yz = cbuffer.data(dd_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xy_zz = cbuffer.data(dd_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xz_xx = cbuffer.data(dd_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xz_xy = cbuffer.data(dd_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xz_xz = cbuffer.data(dd_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xz_yy = cbuffer.data(dd_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xz_yz = cbuffer.data(dd_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xz_zz = cbuffer.data(dd_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_yy_xx = cbuffer.data(dd_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_yy_xy = cbuffer.data(dd_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_yy_xz = cbuffer.data(dd_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_yy_yy = cbuffer.data(dd_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_yy_yz = cbuffer.data(dd_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_yy_zz = cbuffer.data(dd_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_yz_xx = cbuffer.data(dd_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_yz_xy = cbuffer.data(dd_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_yz_xz = cbuffer.data(dd_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_yz_yy = cbuffer.data(dd_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_yz_yz = cbuffer.data(dd_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_yz_zz = cbuffer.data(dd_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_zz_xx = cbuffer.data(dd_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_zz_xy = cbuffer.data(dd_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_zz_xz = cbuffer.data(dd_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_zz_yy = cbuffer.data(dd_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_zz_yz = cbuffer.data(dd_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_zz_zz = cbuffer.data(dd_geom_10_off + 35 * ccomps * dcomps);

            auto g_y_0_xx_xx = cbuffer.data(dd_geom_10_off + 36 * ccomps * dcomps);

            auto g_y_0_xx_xy = cbuffer.data(dd_geom_10_off + 37 * ccomps * dcomps);

            auto g_y_0_xx_xz = cbuffer.data(dd_geom_10_off + 38 * ccomps * dcomps);

            auto g_y_0_xx_yy = cbuffer.data(dd_geom_10_off + 39 * ccomps * dcomps);

            auto g_y_0_xx_yz = cbuffer.data(dd_geom_10_off + 40 * ccomps * dcomps);

            auto g_y_0_xx_zz = cbuffer.data(dd_geom_10_off + 41 * ccomps * dcomps);

            auto g_y_0_xy_xx = cbuffer.data(dd_geom_10_off + 42 * ccomps * dcomps);

            auto g_y_0_xy_xy = cbuffer.data(dd_geom_10_off + 43 * ccomps * dcomps);

            auto g_y_0_xy_xz = cbuffer.data(dd_geom_10_off + 44 * ccomps * dcomps);

            auto g_y_0_xy_yy = cbuffer.data(dd_geom_10_off + 45 * ccomps * dcomps);

            auto g_y_0_xy_yz = cbuffer.data(dd_geom_10_off + 46 * ccomps * dcomps);

            auto g_y_0_xy_zz = cbuffer.data(dd_geom_10_off + 47 * ccomps * dcomps);

            auto g_y_0_xz_xx = cbuffer.data(dd_geom_10_off + 48 * ccomps * dcomps);

            auto g_y_0_xz_xy = cbuffer.data(dd_geom_10_off + 49 * ccomps * dcomps);

            auto g_y_0_xz_xz = cbuffer.data(dd_geom_10_off + 50 * ccomps * dcomps);

            auto g_y_0_xz_yy = cbuffer.data(dd_geom_10_off + 51 * ccomps * dcomps);

            auto g_y_0_xz_yz = cbuffer.data(dd_geom_10_off + 52 * ccomps * dcomps);

            auto g_y_0_xz_zz = cbuffer.data(dd_geom_10_off + 53 * ccomps * dcomps);

            auto g_y_0_yy_xx = cbuffer.data(dd_geom_10_off + 54 * ccomps * dcomps);

            auto g_y_0_yy_xy = cbuffer.data(dd_geom_10_off + 55 * ccomps * dcomps);

            auto g_y_0_yy_xz = cbuffer.data(dd_geom_10_off + 56 * ccomps * dcomps);

            auto g_y_0_yy_yy = cbuffer.data(dd_geom_10_off + 57 * ccomps * dcomps);

            auto g_y_0_yy_yz = cbuffer.data(dd_geom_10_off + 58 * ccomps * dcomps);

            auto g_y_0_yy_zz = cbuffer.data(dd_geom_10_off + 59 * ccomps * dcomps);

            auto g_y_0_yz_xx = cbuffer.data(dd_geom_10_off + 60 * ccomps * dcomps);

            auto g_y_0_yz_xy = cbuffer.data(dd_geom_10_off + 61 * ccomps * dcomps);

            auto g_y_0_yz_xz = cbuffer.data(dd_geom_10_off + 62 * ccomps * dcomps);

            auto g_y_0_yz_yy = cbuffer.data(dd_geom_10_off + 63 * ccomps * dcomps);

            auto g_y_0_yz_yz = cbuffer.data(dd_geom_10_off + 64 * ccomps * dcomps);

            auto g_y_0_yz_zz = cbuffer.data(dd_geom_10_off + 65 * ccomps * dcomps);

            auto g_y_0_zz_xx = cbuffer.data(dd_geom_10_off + 66 * ccomps * dcomps);

            auto g_y_0_zz_xy = cbuffer.data(dd_geom_10_off + 67 * ccomps * dcomps);

            auto g_y_0_zz_xz = cbuffer.data(dd_geom_10_off + 68 * ccomps * dcomps);

            auto g_y_0_zz_yy = cbuffer.data(dd_geom_10_off + 69 * ccomps * dcomps);

            auto g_y_0_zz_yz = cbuffer.data(dd_geom_10_off + 70 * ccomps * dcomps);

            auto g_y_0_zz_zz = cbuffer.data(dd_geom_10_off + 71 * ccomps * dcomps);

            auto g_z_0_xx_xx = cbuffer.data(dd_geom_10_off + 72 * ccomps * dcomps);

            auto g_z_0_xx_xy = cbuffer.data(dd_geom_10_off + 73 * ccomps * dcomps);

            auto g_z_0_xx_xz = cbuffer.data(dd_geom_10_off + 74 * ccomps * dcomps);

            auto g_z_0_xx_yy = cbuffer.data(dd_geom_10_off + 75 * ccomps * dcomps);

            auto g_z_0_xx_yz = cbuffer.data(dd_geom_10_off + 76 * ccomps * dcomps);

            auto g_z_0_xx_zz = cbuffer.data(dd_geom_10_off + 77 * ccomps * dcomps);

            auto g_z_0_xy_xx = cbuffer.data(dd_geom_10_off + 78 * ccomps * dcomps);

            auto g_z_0_xy_xy = cbuffer.data(dd_geom_10_off + 79 * ccomps * dcomps);

            auto g_z_0_xy_xz = cbuffer.data(dd_geom_10_off + 80 * ccomps * dcomps);

            auto g_z_0_xy_yy = cbuffer.data(dd_geom_10_off + 81 * ccomps * dcomps);

            auto g_z_0_xy_yz = cbuffer.data(dd_geom_10_off + 82 * ccomps * dcomps);

            auto g_z_0_xy_zz = cbuffer.data(dd_geom_10_off + 83 * ccomps * dcomps);

            auto g_z_0_xz_xx = cbuffer.data(dd_geom_10_off + 84 * ccomps * dcomps);

            auto g_z_0_xz_xy = cbuffer.data(dd_geom_10_off + 85 * ccomps * dcomps);

            auto g_z_0_xz_xz = cbuffer.data(dd_geom_10_off + 86 * ccomps * dcomps);

            auto g_z_0_xz_yy = cbuffer.data(dd_geom_10_off + 87 * ccomps * dcomps);

            auto g_z_0_xz_yz = cbuffer.data(dd_geom_10_off + 88 * ccomps * dcomps);

            auto g_z_0_xz_zz = cbuffer.data(dd_geom_10_off + 89 * ccomps * dcomps);

            auto g_z_0_yy_xx = cbuffer.data(dd_geom_10_off + 90 * ccomps * dcomps);

            auto g_z_0_yy_xy = cbuffer.data(dd_geom_10_off + 91 * ccomps * dcomps);

            auto g_z_0_yy_xz = cbuffer.data(dd_geom_10_off + 92 * ccomps * dcomps);

            auto g_z_0_yy_yy = cbuffer.data(dd_geom_10_off + 93 * ccomps * dcomps);

            auto g_z_0_yy_yz = cbuffer.data(dd_geom_10_off + 94 * ccomps * dcomps);

            auto g_z_0_yy_zz = cbuffer.data(dd_geom_10_off + 95 * ccomps * dcomps);

            auto g_z_0_yz_xx = cbuffer.data(dd_geom_10_off + 96 * ccomps * dcomps);

            auto g_z_0_yz_xy = cbuffer.data(dd_geom_10_off + 97 * ccomps * dcomps);

            auto g_z_0_yz_xz = cbuffer.data(dd_geom_10_off + 98 * ccomps * dcomps);

            auto g_z_0_yz_yy = cbuffer.data(dd_geom_10_off + 99 * ccomps * dcomps);

            auto g_z_0_yz_yz = cbuffer.data(dd_geom_10_off + 100 * ccomps * dcomps);

            auto g_z_0_yz_zz = cbuffer.data(dd_geom_10_off + 101 * ccomps * dcomps);

            auto g_z_0_zz_xx = cbuffer.data(dd_geom_10_off + 102 * ccomps * dcomps);

            auto g_z_0_zz_xy = cbuffer.data(dd_geom_10_off + 103 * ccomps * dcomps);

            auto g_z_0_zz_xz = cbuffer.data(dd_geom_10_off + 104 * ccomps * dcomps);

            auto g_z_0_zz_yy = cbuffer.data(dd_geom_10_off + 105 * ccomps * dcomps);

            auto g_z_0_zz_yz = cbuffer.data(dd_geom_10_off + 106 * ccomps * dcomps);

            auto g_z_0_zz_zz = cbuffer.data(dd_geom_10_off + 107 * ccomps * dcomps);

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

            /// set up bra offset for contr_buffer_fdxx

            const auto fd_geom_10_off = idx_geom_10_fdxx + i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxx_xx = cbuffer.data(fd_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxx_xy = cbuffer.data(fd_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxx_xz = cbuffer.data(fd_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxx_yy = cbuffer.data(fd_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxx_yz = cbuffer.data(fd_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxx_zz = cbuffer.data(fd_geom_10_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xx_xx, g_x_0_xx_xxx, g_x_0_xx_xxy, g_x_0_xx_xxz, g_x_0_xx_xy, g_x_0_xx_xyy, g_x_0_xx_xyz, g_x_0_xx_xz, g_x_0_xx_xzz, g_x_0_xx_yy, g_x_0_xx_yz, g_x_0_xx_zz, g_x_0_xxx_xx, g_x_0_xxx_xy, g_x_0_xxx_xz, g_x_0_xxx_yy, g_x_0_xxx_yz, g_x_0_xxx_zz, g_xx_xx, g_xx_xy, g_xx_xz, g_xx_yy, g_xx_yz, g_xx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxx_xx[k] = -g_xx_xx[k] - g_x_0_xx_xx[k] * ab_x + g_x_0_xx_xxx[k];

                g_x_0_xxx_xy[k] = -g_xx_xy[k] - g_x_0_xx_xy[k] * ab_x + g_x_0_xx_xxy[k];

                g_x_0_xxx_xz[k] = -g_xx_xz[k] - g_x_0_xx_xz[k] * ab_x + g_x_0_xx_xxz[k];

                g_x_0_xxx_yy[k] = -g_xx_yy[k] - g_x_0_xx_yy[k] * ab_x + g_x_0_xx_xyy[k];

                g_x_0_xxx_yz[k] = -g_xx_yz[k] - g_x_0_xx_yz[k] * ab_x + g_x_0_xx_xyz[k];

                g_x_0_xxx_zz[k] = -g_xx_zz[k] - g_x_0_xx_zz[k] * ab_x + g_x_0_xx_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxy_xx = cbuffer.data(fd_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxy_xy = cbuffer.data(fd_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxy_xz = cbuffer.data(fd_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxy_yy = cbuffer.data(fd_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxy_yz = cbuffer.data(fd_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxy_zz = cbuffer.data(fd_geom_10_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xx_xx, g_x_0_xx_xxy, g_x_0_xx_xy, g_x_0_xx_xyy, g_x_0_xx_xyz, g_x_0_xx_xz, g_x_0_xx_yy, g_x_0_xx_yyy, g_x_0_xx_yyz, g_x_0_xx_yz, g_x_0_xx_yzz, g_x_0_xx_zz, g_x_0_xxy_xx, g_x_0_xxy_xy, g_x_0_xxy_xz, g_x_0_xxy_yy, g_x_0_xxy_yz, g_x_0_xxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxy_xx[k] = -g_x_0_xx_xx[k] * ab_y + g_x_0_xx_xxy[k];

                g_x_0_xxy_xy[k] = -g_x_0_xx_xy[k] * ab_y + g_x_0_xx_xyy[k];

                g_x_0_xxy_xz[k] = -g_x_0_xx_xz[k] * ab_y + g_x_0_xx_xyz[k];

                g_x_0_xxy_yy[k] = -g_x_0_xx_yy[k] * ab_y + g_x_0_xx_yyy[k];

                g_x_0_xxy_yz[k] = -g_x_0_xx_yz[k] * ab_y + g_x_0_xx_yyz[k];

                g_x_0_xxy_zz[k] = -g_x_0_xx_zz[k] * ab_y + g_x_0_xx_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxz_xx = cbuffer.data(fd_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxz_xy = cbuffer.data(fd_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxz_xz = cbuffer.data(fd_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxz_yy = cbuffer.data(fd_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxz_yz = cbuffer.data(fd_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxz_zz = cbuffer.data(fd_geom_10_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xx_xx, g_x_0_xx_xxz, g_x_0_xx_xy, g_x_0_xx_xyz, g_x_0_xx_xz, g_x_0_xx_xzz, g_x_0_xx_yy, g_x_0_xx_yyz, g_x_0_xx_yz, g_x_0_xx_yzz, g_x_0_xx_zz, g_x_0_xx_zzz, g_x_0_xxz_xx, g_x_0_xxz_xy, g_x_0_xxz_xz, g_x_0_xxz_yy, g_x_0_xxz_yz, g_x_0_xxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxz_xx[k] = -g_x_0_xx_xx[k] * ab_z + g_x_0_xx_xxz[k];

                g_x_0_xxz_xy[k] = -g_x_0_xx_xy[k] * ab_z + g_x_0_xx_xyz[k];

                g_x_0_xxz_xz[k] = -g_x_0_xx_xz[k] * ab_z + g_x_0_xx_xzz[k];

                g_x_0_xxz_yy[k] = -g_x_0_xx_yy[k] * ab_z + g_x_0_xx_yyz[k];

                g_x_0_xxz_yz[k] = -g_x_0_xx_yz[k] * ab_z + g_x_0_xx_yzz[k];

                g_x_0_xxz_zz[k] = -g_x_0_xx_zz[k] * ab_z + g_x_0_xx_zzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyy_xx = cbuffer.data(fd_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xyy_xy = cbuffer.data(fd_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xyy_xz = cbuffer.data(fd_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xyy_yy = cbuffer.data(fd_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xyy_yz = cbuffer.data(fd_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xyy_zz = cbuffer.data(fd_geom_10_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xy_xx, g_x_0_xy_xxy, g_x_0_xy_xy, g_x_0_xy_xyy, g_x_0_xy_xyz, g_x_0_xy_xz, g_x_0_xy_yy, g_x_0_xy_yyy, g_x_0_xy_yyz, g_x_0_xy_yz, g_x_0_xy_yzz, g_x_0_xy_zz, g_x_0_xyy_xx, g_x_0_xyy_xy, g_x_0_xyy_xz, g_x_0_xyy_yy, g_x_0_xyy_yz, g_x_0_xyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyy_xx[k] = -g_x_0_xy_xx[k] * ab_y + g_x_0_xy_xxy[k];

                g_x_0_xyy_xy[k] = -g_x_0_xy_xy[k] * ab_y + g_x_0_xy_xyy[k];

                g_x_0_xyy_xz[k] = -g_x_0_xy_xz[k] * ab_y + g_x_0_xy_xyz[k];

                g_x_0_xyy_yy[k] = -g_x_0_xy_yy[k] * ab_y + g_x_0_xy_yyy[k];

                g_x_0_xyy_yz[k] = -g_x_0_xy_yz[k] * ab_y + g_x_0_xy_yyz[k];

                g_x_0_xyy_zz[k] = -g_x_0_xy_zz[k] * ab_y + g_x_0_xy_yzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyz_xx = cbuffer.data(fd_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xyz_xy = cbuffer.data(fd_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xyz_xz = cbuffer.data(fd_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xyz_yy = cbuffer.data(fd_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xyz_yz = cbuffer.data(fd_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xyz_zz = cbuffer.data(fd_geom_10_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyz_xx, g_x_0_xyz_xy, g_x_0_xyz_xz, g_x_0_xyz_yy, g_x_0_xyz_yz, g_x_0_xyz_zz, g_x_0_xz_xx, g_x_0_xz_xxy, g_x_0_xz_xy, g_x_0_xz_xyy, g_x_0_xz_xyz, g_x_0_xz_xz, g_x_0_xz_yy, g_x_0_xz_yyy, g_x_0_xz_yyz, g_x_0_xz_yz, g_x_0_xz_yzz, g_x_0_xz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyz_xx[k] = -g_x_0_xz_xx[k] * ab_y + g_x_0_xz_xxy[k];

                g_x_0_xyz_xy[k] = -g_x_0_xz_xy[k] * ab_y + g_x_0_xz_xyy[k];

                g_x_0_xyz_xz[k] = -g_x_0_xz_xz[k] * ab_y + g_x_0_xz_xyz[k];

                g_x_0_xyz_yy[k] = -g_x_0_xz_yy[k] * ab_y + g_x_0_xz_yyy[k];

                g_x_0_xyz_yz[k] = -g_x_0_xz_yz[k] * ab_y + g_x_0_xz_yyz[k];

                g_x_0_xyz_zz[k] = -g_x_0_xz_zz[k] * ab_y + g_x_0_xz_yzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzz_xx = cbuffer.data(fd_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xzz_xy = cbuffer.data(fd_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xzz_xz = cbuffer.data(fd_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xzz_yy = cbuffer.data(fd_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xzz_yz = cbuffer.data(fd_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xzz_zz = cbuffer.data(fd_geom_10_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xz_xx, g_x_0_xz_xxz, g_x_0_xz_xy, g_x_0_xz_xyz, g_x_0_xz_xz, g_x_0_xz_xzz, g_x_0_xz_yy, g_x_0_xz_yyz, g_x_0_xz_yz, g_x_0_xz_yzz, g_x_0_xz_zz, g_x_0_xz_zzz, g_x_0_xzz_xx, g_x_0_xzz_xy, g_x_0_xzz_xz, g_x_0_xzz_yy, g_x_0_xzz_yz, g_x_0_xzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzz_xx[k] = -g_x_0_xz_xx[k] * ab_z + g_x_0_xz_xxz[k];

                g_x_0_xzz_xy[k] = -g_x_0_xz_xy[k] * ab_z + g_x_0_xz_xyz[k];

                g_x_0_xzz_xz[k] = -g_x_0_xz_xz[k] * ab_z + g_x_0_xz_xzz[k];

                g_x_0_xzz_yy[k] = -g_x_0_xz_yy[k] * ab_z + g_x_0_xz_yyz[k];

                g_x_0_xzz_yz[k] = -g_x_0_xz_yz[k] * ab_z + g_x_0_xz_yzz[k];

                g_x_0_xzz_zz[k] = -g_x_0_xz_zz[k] * ab_z + g_x_0_xz_zzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyy_xx = cbuffer.data(fd_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_yyy_xy = cbuffer.data(fd_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_yyy_xz = cbuffer.data(fd_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_yyy_yy = cbuffer.data(fd_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_yyy_yz = cbuffer.data(fd_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_yyy_zz = cbuffer.data(fd_geom_10_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yy_xx, g_x_0_yy_xxy, g_x_0_yy_xy, g_x_0_yy_xyy, g_x_0_yy_xyz, g_x_0_yy_xz, g_x_0_yy_yy, g_x_0_yy_yyy, g_x_0_yy_yyz, g_x_0_yy_yz, g_x_0_yy_yzz, g_x_0_yy_zz, g_x_0_yyy_xx, g_x_0_yyy_xy, g_x_0_yyy_xz, g_x_0_yyy_yy, g_x_0_yyy_yz, g_x_0_yyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyy_xx[k] = -g_x_0_yy_xx[k] * ab_y + g_x_0_yy_xxy[k];

                g_x_0_yyy_xy[k] = -g_x_0_yy_xy[k] * ab_y + g_x_0_yy_xyy[k];

                g_x_0_yyy_xz[k] = -g_x_0_yy_xz[k] * ab_y + g_x_0_yy_xyz[k];

                g_x_0_yyy_yy[k] = -g_x_0_yy_yy[k] * ab_y + g_x_0_yy_yyy[k];

                g_x_0_yyy_yz[k] = -g_x_0_yy_yz[k] * ab_y + g_x_0_yy_yyz[k];

                g_x_0_yyy_zz[k] = -g_x_0_yy_zz[k] * ab_y + g_x_0_yy_yzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyz_xx = cbuffer.data(fd_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_yyz_xy = cbuffer.data(fd_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_yyz_xz = cbuffer.data(fd_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_yyz_yy = cbuffer.data(fd_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_yyz_yz = cbuffer.data(fd_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_yyz_zz = cbuffer.data(fd_geom_10_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyz_xx, g_x_0_yyz_xy, g_x_0_yyz_xz, g_x_0_yyz_yy, g_x_0_yyz_yz, g_x_0_yyz_zz, g_x_0_yz_xx, g_x_0_yz_xxy, g_x_0_yz_xy, g_x_0_yz_xyy, g_x_0_yz_xyz, g_x_0_yz_xz, g_x_0_yz_yy, g_x_0_yz_yyy, g_x_0_yz_yyz, g_x_0_yz_yz, g_x_0_yz_yzz, g_x_0_yz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyz_xx[k] = -g_x_0_yz_xx[k] * ab_y + g_x_0_yz_xxy[k];

                g_x_0_yyz_xy[k] = -g_x_0_yz_xy[k] * ab_y + g_x_0_yz_xyy[k];

                g_x_0_yyz_xz[k] = -g_x_0_yz_xz[k] * ab_y + g_x_0_yz_xyz[k];

                g_x_0_yyz_yy[k] = -g_x_0_yz_yy[k] * ab_y + g_x_0_yz_yyy[k];

                g_x_0_yyz_yz[k] = -g_x_0_yz_yz[k] * ab_y + g_x_0_yz_yyz[k];

                g_x_0_yyz_zz[k] = -g_x_0_yz_zz[k] * ab_y + g_x_0_yz_yzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzz_xx = cbuffer.data(fd_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_yzz_xy = cbuffer.data(fd_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_yzz_xz = cbuffer.data(fd_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_yzz_yy = cbuffer.data(fd_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_yzz_yz = cbuffer.data(fd_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_yzz_zz = cbuffer.data(fd_geom_10_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yzz_xx, g_x_0_yzz_xy, g_x_0_yzz_xz, g_x_0_yzz_yy, g_x_0_yzz_yz, g_x_0_yzz_zz, g_x_0_zz_xx, g_x_0_zz_xxy, g_x_0_zz_xy, g_x_0_zz_xyy, g_x_0_zz_xyz, g_x_0_zz_xz, g_x_0_zz_yy, g_x_0_zz_yyy, g_x_0_zz_yyz, g_x_0_zz_yz, g_x_0_zz_yzz, g_x_0_zz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzz_xx[k] = -g_x_0_zz_xx[k] * ab_y + g_x_0_zz_xxy[k];

                g_x_0_yzz_xy[k] = -g_x_0_zz_xy[k] * ab_y + g_x_0_zz_xyy[k];

                g_x_0_yzz_xz[k] = -g_x_0_zz_xz[k] * ab_y + g_x_0_zz_xyz[k];

                g_x_0_yzz_yy[k] = -g_x_0_zz_yy[k] * ab_y + g_x_0_zz_yyy[k];

                g_x_0_yzz_yz[k] = -g_x_0_zz_yz[k] * ab_y + g_x_0_zz_yyz[k];

                g_x_0_yzz_zz[k] = -g_x_0_zz_zz[k] * ab_y + g_x_0_zz_yzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzz_xx = cbuffer.data(fd_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_zzz_xy = cbuffer.data(fd_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_zzz_xz = cbuffer.data(fd_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_zzz_yy = cbuffer.data(fd_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_zzz_yz = cbuffer.data(fd_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_zzz_zz = cbuffer.data(fd_geom_10_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zz_xx, g_x_0_zz_xxz, g_x_0_zz_xy, g_x_0_zz_xyz, g_x_0_zz_xz, g_x_0_zz_xzz, g_x_0_zz_yy, g_x_0_zz_yyz, g_x_0_zz_yz, g_x_0_zz_yzz, g_x_0_zz_zz, g_x_0_zz_zzz, g_x_0_zzz_xx, g_x_0_zzz_xy, g_x_0_zzz_xz, g_x_0_zzz_yy, g_x_0_zzz_yz, g_x_0_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzz_xx[k] = -g_x_0_zz_xx[k] * ab_z + g_x_0_zz_xxz[k];

                g_x_0_zzz_xy[k] = -g_x_0_zz_xy[k] * ab_z + g_x_0_zz_xyz[k];

                g_x_0_zzz_xz[k] = -g_x_0_zz_xz[k] * ab_z + g_x_0_zz_xzz[k];

                g_x_0_zzz_yy[k] = -g_x_0_zz_yy[k] * ab_z + g_x_0_zz_yyz[k];

                g_x_0_zzz_yz[k] = -g_x_0_zz_yz[k] * ab_z + g_x_0_zz_yzz[k];

                g_x_0_zzz_zz[k] = -g_x_0_zz_zz[k] * ab_z + g_x_0_zz_zzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxx_xx = cbuffer.data(fd_geom_10_off + 60 * ccomps * dcomps);

            auto g_y_0_xxx_xy = cbuffer.data(fd_geom_10_off + 61 * ccomps * dcomps);

            auto g_y_0_xxx_xz = cbuffer.data(fd_geom_10_off + 62 * ccomps * dcomps);

            auto g_y_0_xxx_yy = cbuffer.data(fd_geom_10_off + 63 * ccomps * dcomps);

            auto g_y_0_xxx_yz = cbuffer.data(fd_geom_10_off + 64 * ccomps * dcomps);

            auto g_y_0_xxx_zz = cbuffer.data(fd_geom_10_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xx_xx, g_y_0_xx_xxx, g_y_0_xx_xxy, g_y_0_xx_xxz, g_y_0_xx_xy, g_y_0_xx_xyy, g_y_0_xx_xyz, g_y_0_xx_xz, g_y_0_xx_xzz, g_y_0_xx_yy, g_y_0_xx_yz, g_y_0_xx_zz, g_y_0_xxx_xx, g_y_0_xxx_xy, g_y_0_xxx_xz, g_y_0_xxx_yy, g_y_0_xxx_yz, g_y_0_xxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxx_xx[k] = -g_y_0_xx_xx[k] * ab_x + g_y_0_xx_xxx[k];

                g_y_0_xxx_xy[k] = -g_y_0_xx_xy[k] * ab_x + g_y_0_xx_xxy[k];

                g_y_0_xxx_xz[k] = -g_y_0_xx_xz[k] * ab_x + g_y_0_xx_xxz[k];

                g_y_0_xxx_yy[k] = -g_y_0_xx_yy[k] * ab_x + g_y_0_xx_xyy[k];

                g_y_0_xxx_yz[k] = -g_y_0_xx_yz[k] * ab_x + g_y_0_xx_xyz[k];

                g_y_0_xxx_zz[k] = -g_y_0_xx_zz[k] * ab_x + g_y_0_xx_xzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxy_xx = cbuffer.data(fd_geom_10_off + 66 * ccomps * dcomps);

            auto g_y_0_xxy_xy = cbuffer.data(fd_geom_10_off + 67 * ccomps * dcomps);

            auto g_y_0_xxy_xz = cbuffer.data(fd_geom_10_off + 68 * ccomps * dcomps);

            auto g_y_0_xxy_yy = cbuffer.data(fd_geom_10_off + 69 * ccomps * dcomps);

            auto g_y_0_xxy_yz = cbuffer.data(fd_geom_10_off + 70 * ccomps * dcomps);

            auto g_y_0_xxy_zz = cbuffer.data(fd_geom_10_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxy_xx, g_y_0_xxy_xy, g_y_0_xxy_xz, g_y_0_xxy_yy, g_y_0_xxy_yz, g_y_0_xxy_zz, g_y_0_xy_xx, g_y_0_xy_xxx, g_y_0_xy_xxy, g_y_0_xy_xxz, g_y_0_xy_xy, g_y_0_xy_xyy, g_y_0_xy_xyz, g_y_0_xy_xz, g_y_0_xy_xzz, g_y_0_xy_yy, g_y_0_xy_yz, g_y_0_xy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxy_xx[k] = -g_y_0_xy_xx[k] * ab_x + g_y_0_xy_xxx[k];

                g_y_0_xxy_xy[k] = -g_y_0_xy_xy[k] * ab_x + g_y_0_xy_xxy[k];

                g_y_0_xxy_xz[k] = -g_y_0_xy_xz[k] * ab_x + g_y_0_xy_xxz[k];

                g_y_0_xxy_yy[k] = -g_y_0_xy_yy[k] * ab_x + g_y_0_xy_xyy[k];

                g_y_0_xxy_yz[k] = -g_y_0_xy_yz[k] * ab_x + g_y_0_xy_xyz[k];

                g_y_0_xxy_zz[k] = -g_y_0_xy_zz[k] * ab_x + g_y_0_xy_xzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxz_xx = cbuffer.data(fd_geom_10_off + 72 * ccomps * dcomps);

            auto g_y_0_xxz_xy = cbuffer.data(fd_geom_10_off + 73 * ccomps * dcomps);

            auto g_y_0_xxz_xz = cbuffer.data(fd_geom_10_off + 74 * ccomps * dcomps);

            auto g_y_0_xxz_yy = cbuffer.data(fd_geom_10_off + 75 * ccomps * dcomps);

            auto g_y_0_xxz_yz = cbuffer.data(fd_geom_10_off + 76 * ccomps * dcomps);

            auto g_y_0_xxz_zz = cbuffer.data(fd_geom_10_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxz_xx, g_y_0_xxz_xy, g_y_0_xxz_xz, g_y_0_xxz_yy, g_y_0_xxz_yz, g_y_0_xxz_zz, g_y_0_xz_xx, g_y_0_xz_xxx, g_y_0_xz_xxy, g_y_0_xz_xxz, g_y_0_xz_xy, g_y_0_xz_xyy, g_y_0_xz_xyz, g_y_0_xz_xz, g_y_0_xz_xzz, g_y_0_xz_yy, g_y_0_xz_yz, g_y_0_xz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxz_xx[k] = -g_y_0_xz_xx[k] * ab_x + g_y_0_xz_xxx[k];

                g_y_0_xxz_xy[k] = -g_y_0_xz_xy[k] * ab_x + g_y_0_xz_xxy[k];

                g_y_0_xxz_xz[k] = -g_y_0_xz_xz[k] * ab_x + g_y_0_xz_xxz[k];

                g_y_0_xxz_yy[k] = -g_y_0_xz_yy[k] * ab_x + g_y_0_xz_xyy[k];

                g_y_0_xxz_yz[k] = -g_y_0_xz_yz[k] * ab_x + g_y_0_xz_xyz[k];

                g_y_0_xxz_zz[k] = -g_y_0_xz_zz[k] * ab_x + g_y_0_xz_xzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyy_xx = cbuffer.data(fd_geom_10_off + 78 * ccomps * dcomps);

            auto g_y_0_xyy_xy = cbuffer.data(fd_geom_10_off + 79 * ccomps * dcomps);

            auto g_y_0_xyy_xz = cbuffer.data(fd_geom_10_off + 80 * ccomps * dcomps);

            auto g_y_0_xyy_yy = cbuffer.data(fd_geom_10_off + 81 * ccomps * dcomps);

            auto g_y_0_xyy_yz = cbuffer.data(fd_geom_10_off + 82 * ccomps * dcomps);

            auto g_y_0_xyy_zz = cbuffer.data(fd_geom_10_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyy_xx, g_y_0_xyy_xy, g_y_0_xyy_xz, g_y_0_xyy_yy, g_y_0_xyy_yz, g_y_0_xyy_zz, g_y_0_yy_xx, g_y_0_yy_xxx, g_y_0_yy_xxy, g_y_0_yy_xxz, g_y_0_yy_xy, g_y_0_yy_xyy, g_y_0_yy_xyz, g_y_0_yy_xz, g_y_0_yy_xzz, g_y_0_yy_yy, g_y_0_yy_yz, g_y_0_yy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyy_xx[k] = -g_y_0_yy_xx[k] * ab_x + g_y_0_yy_xxx[k];

                g_y_0_xyy_xy[k] = -g_y_0_yy_xy[k] * ab_x + g_y_0_yy_xxy[k];

                g_y_0_xyy_xz[k] = -g_y_0_yy_xz[k] * ab_x + g_y_0_yy_xxz[k];

                g_y_0_xyy_yy[k] = -g_y_0_yy_yy[k] * ab_x + g_y_0_yy_xyy[k];

                g_y_0_xyy_yz[k] = -g_y_0_yy_yz[k] * ab_x + g_y_0_yy_xyz[k];

                g_y_0_xyy_zz[k] = -g_y_0_yy_zz[k] * ab_x + g_y_0_yy_xzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyz_xx = cbuffer.data(fd_geom_10_off + 84 * ccomps * dcomps);

            auto g_y_0_xyz_xy = cbuffer.data(fd_geom_10_off + 85 * ccomps * dcomps);

            auto g_y_0_xyz_xz = cbuffer.data(fd_geom_10_off + 86 * ccomps * dcomps);

            auto g_y_0_xyz_yy = cbuffer.data(fd_geom_10_off + 87 * ccomps * dcomps);

            auto g_y_0_xyz_yz = cbuffer.data(fd_geom_10_off + 88 * ccomps * dcomps);

            auto g_y_0_xyz_zz = cbuffer.data(fd_geom_10_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyz_xx, g_y_0_xyz_xy, g_y_0_xyz_xz, g_y_0_xyz_yy, g_y_0_xyz_yz, g_y_0_xyz_zz, g_y_0_yz_xx, g_y_0_yz_xxx, g_y_0_yz_xxy, g_y_0_yz_xxz, g_y_0_yz_xy, g_y_0_yz_xyy, g_y_0_yz_xyz, g_y_0_yz_xz, g_y_0_yz_xzz, g_y_0_yz_yy, g_y_0_yz_yz, g_y_0_yz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyz_xx[k] = -g_y_0_yz_xx[k] * ab_x + g_y_0_yz_xxx[k];

                g_y_0_xyz_xy[k] = -g_y_0_yz_xy[k] * ab_x + g_y_0_yz_xxy[k];

                g_y_0_xyz_xz[k] = -g_y_0_yz_xz[k] * ab_x + g_y_0_yz_xxz[k];

                g_y_0_xyz_yy[k] = -g_y_0_yz_yy[k] * ab_x + g_y_0_yz_xyy[k];

                g_y_0_xyz_yz[k] = -g_y_0_yz_yz[k] * ab_x + g_y_0_yz_xyz[k];

                g_y_0_xyz_zz[k] = -g_y_0_yz_zz[k] * ab_x + g_y_0_yz_xzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzz_xx = cbuffer.data(fd_geom_10_off + 90 * ccomps * dcomps);

            auto g_y_0_xzz_xy = cbuffer.data(fd_geom_10_off + 91 * ccomps * dcomps);

            auto g_y_0_xzz_xz = cbuffer.data(fd_geom_10_off + 92 * ccomps * dcomps);

            auto g_y_0_xzz_yy = cbuffer.data(fd_geom_10_off + 93 * ccomps * dcomps);

            auto g_y_0_xzz_yz = cbuffer.data(fd_geom_10_off + 94 * ccomps * dcomps);

            auto g_y_0_xzz_zz = cbuffer.data(fd_geom_10_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xzz_xx, g_y_0_xzz_xy, g_y_0_xzz_xz, g_y_0_xzz_yy, g_y_0_xzz_yz, g_y_0_xzz_zz, g_y_0_zz_xx, g_y_0_zz_xxx, g_y_0_zz_xxy, g_y_0_zz_xxz, g_y_0_zz_xy, g_y_0_zz_xyy, g_y_0_zz_xyz, g_y_0_zz_xz, g_y_0_zz_xzz, g_y_0_zz_yy, g_y_0_zz_yz, g_y_0_zz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzz_xx[k] = -g_y_0_zz_xx[k] * ab_x + g_y_0_zz_xxx[k];

                g_y_0_xzz_xy[k] = -g_y_0_zz_xy[k] * ab_x + g_y_0_zz_xxy[k];

                g_y_0_xzz_xz[k] = -g_y_0_zz_xz[k] * ab_x + g_y_0_zz_xxz[k];

                g_y_0_xzz_yy[k] = -g_y_0_zz_yy[k] * ab_x + g_y_0_zz_xyy[k];

                g_y_0_xzz_yz[k] = -g_y_0_zz_yz[k] * ab_x + g_y_0_zz_xyz[k];

                g_y_0_xzz_zz[k] = -g_y_0_zz_zz[k] * ab_x + g_y_0_zz_xzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyy_xx = cbuffer.data(fd_geom_10_off + 96 * ccomps * dcomps);

            auto g_y_0_yyy_xy = cbuffer.data(fd_geom_10_off + 97 * ccomps * dcomps);

            auto g_y_0_yyy_xz = cbuffer.data(fd_geom_10_off + 98 * ccomps * dcomps);

            auto g_y_0_yyy_yy = cbuffer.data(fd_geom_10_off + 99 * ccomps * dcomps);

            auto g_y_0_yyy_yz = cbuffer.data(fd_geom_10_off + 100 * ccomps * dcomps);

            auto g_y_0_yyy_zz = cbuffer.data(fd_geom_10_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yy_xx, g_y_0_yy_xxy, g_y_0_yy_xy, g_y_0_yy_xyy, g_y_0_yy_xyz, g_y_0_yy_xz, g_y_0_yy_yy, g_y_0_yy_yyy, g_y_0_yy_yyz, g_y_0_yy_yz, g_y_0_yy_yzz, g_y_0_yy_zz, g_y_0_yyy_xx, g_y_0_yyy_xy, g_y_0_yyy_xz, g_y_0_yyy_yy, g_y_0_yyy_yz, g_y_0_yyy_zz, g_yy_xx, g_yy_xy, g_yy_xz, g_yy_yy, g_yy_yz, g_yy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyy_xx[k] = -g_yy_xx[k] - g_y_0_yy_xx[k] * ab_y + g_y_0_yy_xxy[k];

                g_y_0_yyy_xy[k] = -g_yy_xy[k] - g_y_0_yy_xy[k] * ab_y + g_y_0_yy_xyy[k];

                g_y_0_yyy_xz[k] = -g_yy_xz[k] - g_y_0_yy_xz[k] * ab_y + g_y_0_yy_xyz[k];

                g_y_0_yyy_yy[k] = -g_yy_yy[k] - g_y_0_yy_yy[k] * ab_y + g_y_0_yy_yyy[k];

                g_y_0_yyy_yz[k] = -g_yy_yz[k] - g_y_0_yy_yz[k] * ab_y + g_y_0_yy_yyz[k];

                g_y_0_yyy_zz[k] = -g_yy_zz[k] - g_y_0_yy_zz[k] * ab_y + g_y_0_yy_yzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyz_xx = cbuffer.data(fd_geom_10_off + 102 * ccomps * dcomps);

            auto g_y_0_yyz_xy = cbuffer.data(fd_geom_10_off + 103 * ccomps * dcomps);

            auto g_y_0_yyz_xz = cbuffer.data(fd_geom_10_off + 104 * ccomps * dcomps);

            auto g_y_0_yyz_yy = cbuffer.data(fd_geom_10_off + 105 * ccomps * dcomps);

            auto g_y_0_yyz_yz = cbuffer.data(fd_geom_10_off + 106 * ccomps * dcomps);

            auto g_y_0_yyz_zz = cbuffer.data(fd_geom_10_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yy_xx, g_y_0_yy_xxz, g_y_0_yy_xy, g_y_0_yy_xyz, g_y_0_yy_xz, g_y_0_yy_xzz, g_y_0_yy_yy, g_y_0_yy_yyz, g_y_0_yy_yz, g_y_0_yy_yzz, g_y_0_yy_zz, g_y_0_yy_zzz, g_y_0_yyz_xx, g_y_0_yyz_xy, g_y_0_yyz_xz, g_y_0_yyz_yy, g_y_0_yyz_yz, g_y_0_yyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyz_xx[k] = -g_y_0_yy_xx[k] * ab_z + g_y_0_yy_xxz[k];

                g_y_0_yyz_xy[k] = -g_y_0_yy_xy[k] * ab_z + g_y_0_yy_xyz[k];

                g_y_0_yyz_xz[k] = -g_y_0_yy_xz[k] * ab_z + g_y_0_yy_xzz[k];

                g_y_0_yyz_yy[k] = -g_y_0_yy_yy[k] * ab_z + g_y_0_yy_yyz[k];

                g_y_0_yyz_yz[k] = -g_y_0_yy_yz[k] * ab_z + g_y_0_yy_yzz[k];

                g_y_0_yyz_zz[k] = -g_y_0_yy_zz[k] * ab_z + g_y_0_yy_zzz[k];
            }

            /// Set up 108-114 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzz_xx = cbuffer.data(fd_geom_10_off + 108 * ccomps * dcomps);

            auto g_y_0_yzz_xy = cbuffer.data(fd_geom_10_off + 109 * ccomps * dcomps);

            auto g_y_0_yzz_xz = cbuffer.data(fd_geom_10_off + 110 * ccomps * dcomps);

            auto g_y_0_yzz_yy = cbuffer.data(fd_geom_10_off + 111 * ccomps * dcomps);

            auto g_y_0_yzz_yz = cbuffer.data(fd_geom_10_off + 112 * ccomps * dcomps);

            auto g_y_0_yzz_zz = cbuffer.data(fd_geom_10_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yz_xx, g_y_0_yz_xxz, g_y_0_yz_xy, g_y_0_yz_xyz, g_y_0_yz_xz, g_y_0_yz_xzz, g_y_0_yz_yy, g_y_0_yz_yyz, g_y_0_yz_yz, g_y_0_yz_yzz, g_y_0_yz_zz, g_y_0_yz_zzz, g_y_0_yzz_xx, g_y_0_yzz_xy, g_y_0_yzz_xz, g_y_0_yzz_yy, g_y_0_yzz_yz, g_y_0_yzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzz_xx[k] = -g_y_0_yz_xx[k] * ab_z + g_y_0_yz_xxz[k];

                g_y_0_yzz_xy[k] = -g_y_0_yz_xy[k] * ab_z + g_y_0_yz_xyz[k];

                g_y_0_yzz_xz[k] = -g_y_0_yz_xz[k] * ab_z + g_y_0_yz_xzz[k];

                g_y_0_yzz_yy[k] = -g_y_0_yz_yy[k] * ab_z + g_y_0_yz_yyz[k];

                g_y_0_yzz_yz[k] = -g_y_0_yz_yz[k] * ab_z + g_y_0_yz_yzz[k];

                g_y_0_yzz_zz[k] = -g_y_0_yz_zz[k] * ab_z + g_y_0_yz_zzz[k];
            }

            /// Set up 114-120 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzz_xx = cbuffer.data(fd_geom_10_off + 114 * ccomps * dcomps);

            auto g_y_0_zzz_xy = cbuffer.data(fd_geom_10_off + 115 * ccomps * dcomps);

            auto g_y_0_zzz_xz = cbuffer.data(fd_geom_10_off + 116 * ccomps * dcomps);

            auto g_y_0_zzz_yy = cbuffer.data(fd_geom_10_off + 117 * ccomps * dcomps);

            auto g_y_0_zzz_yz = cbuffer.data(fd_geom_10_off + 118 * ccomps * dcomps);

            auto g_y_0_zzz_zz = cbuffer.data(fd_geom_10_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zz_xx, g_y_0_zz_xxz, g_y_0_zz_xy, g_y_0_zz_xyz, g_y_0_zz_xz, g_y_0_zz_xzz, g_y_0_zz_yy, g_y_0_zz_yyz, g_y_0_zz_yz, g_y_0_zz_yzz, g_y_0_zz_zz, g_y_0_zz_zzz, g_y_0_zzz_xx, g_y_0_zzz_xy, g_y_0_zzz_xz, g_y_0_zzz_yy, g_y_0_zzz_yz, g_y_0_zzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzz_xx[k] = -g_y_0_zz_xx[k] * ab_z + g_y_0_zz_xxz[k];

                g_y_0_zzz_xy[k] = -g_y_0_zz_xy[k] * ab_z + g_y_0_zz_xyz[k];

                g_y_0_zzz_xz[k] = -g_y_0_zz_xz[k] * ab_z + g_y_0_zz_xzz[k];

                g_y_0_zzz_yy[k] = -g_y_0_zz_yy[k] * ab_z + g_y_0_zz_yyz[k];

                g_y_0_zzz_yz[k] = -g_y_0_zz_yz[k] * ab_z + g_y_0_zz_yzz[k];

                g_y_0_zzz_zz[k] = -g_y_0_zz_zz[k] * ab_z + g_y_0_zz_zzz[k];
            }

            /// Set up 120-126 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxx_xx = cbuffer.data(fd_geom_10_off + 120 * ccomps * dcomps);

            auto g_z_0_xxx_xy = cbuffer.data(fd_geom_10_off + 121 * ccomps * dcomps);

            auto g_z_0_xxx_xz = cbuffer.data(fd_geom_10_off + 122 * ccomps * dcomps);

            auto g_z_0_xxx_yy = cbuffer.data(fd_geom_10_off + 123 * ccomps * dcomps);

            auto g_z_0_xxx_yz = cbuffer.data(fd_geom_10_off + 124 * ccomps * dcomps);

            auto g_z_0_xxx_zz = cbuffer.data(fd_geom_10_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xx_xx, g_z_0_xx_xxx, g_z_0_xx_xxy, g_z_0_xx_xxz, g_z_0_xx_xy, g_z_0_xx_xyy, g_z_0_xx_xyz, g_z_0_xx_xz, g_z_0_xx_xzz, g_z_0_xx_yy, g_z_0_xx_yz, g_z_0_xx_zz, g_z_0_xxx_xx, g_z_0_xxx_xy, g_z_0_xxx_xz, g_z_0_xxx_yy, g_z_0_xxx_yz, g_z_0_xxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxx_xx[k] = -g_z_0_xx_xx[k] * ab_x + g_z_0_xx_xxx[k];

                g_z_0_xxx_xy[k] = -g_z_0_xx_xy[k] * ab_x + g_z_0_xx_xxy[k];

                g_z_0_xxx_xz[k] = -g_z_0_xx_xz[k] * ab_x + g_z_0_xx_xxz[k];

                g_z_0_xxx_yy[k] = -g_z_0_xx_yy[k] * ab_x + g_z_0_xx_xyy[k];

                g_z_0_xxx_yz[k] = -g_z_0_xx_yz[k] * ab_x + g_z_0_xx_xyz[k];

                g_z_0_xxx_zz[k] = -g_z_0_xx_zz[k] * ab_x + g_z_0_xx_xzz[k];
            }

            /// Set up 126-132 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxy_xx = cbuffer.data(fd_geom_10_off + 126 * ccomps * dcomps);

            auto g_z_0_xxy_xy = cbuffer.data(fd_geom_10_off + 127 * ccomps * dcomps);

            auto g_z_0_xxy_xz = cbuffer.data(fd_geom_10_off + 128 * ccomps * dcomps);

            auto g_z_0_xxy_yy = cbuffer.data(fd_geom_10_off + 129 * ccomps * dcomps);

            auto g_z_0_xxy_yz = cbuffer.data(fd_geom_10_off + 130 * ccomps * dcomps);

            auto g_z_0_xxy_zz = cbuffer.data(fd_geom_10_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxy_xx, g_z_0_xxy_xy, g_z_0_xxy_xz, g_z_0_xxy_yy, g_z_0_xxy_yz, g_z_0_xxy_zz, g_z_0_xy_xx, g_z_0_xy_xxx, g_z_0_xy_xxy, g_z_0_xy_xxz, g_z_0_xy_xy, g_z_0_xy_xyy, g_z_0_xy_xyz, g_z_0_xy_xz, g_z_0_xy_xzz, g_z_0_xy_yy, g_z_0_xy_yz, g_z_0_xy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxy_xx[k] = -g_z_0_xy_xx[k] * ab_x + g_z_0_xy_xxx[k];

                g_z_0_xxy_xy[k] = -g_z_0_xy_xy[k] * ab_x + g_z_0_xy_xxy[k];

                g_z_0_xxy_xz[k] = -g_z_0_xy_xz[k] * ab_x + g_z_0_xy_xxz[k];

                g_z_0_xxy_yy[k] = -g_z_0_xy_yy[k] * ab_x + g_z_0_xy_xyy[k];

                g_z_0_xxy_yz[k] = -g_z_0_xy_yz[k] * ab_x + g_z_0_xy_xyz[k];

                g_z_0_xxy_zz[k] = -g_z_0_xy_zz[k] * ab_x + g_z_0_xy_xzz[k];
            }

            /// Set up 132-138 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxz_xx = cbuffer.data(fd_geom_10_off + 132 * ccomps * dcomps);

            auto g_z_0_xxz_xy = cbuffer.data(fd_geom_10_off + 133 * ccomps * dcomps);

            auto g_z_0_xxz_xz = cbuffer.data(fd_geom_10_off + 134 * ccomps * dcomps);

            auto g_z_0_xxz_yy = cbuffer.data(fd_geom_10_off + 135 * ccomps * dcomps);

            auto g_z_0_xxz_yz = cbuffer.data(fd_geom_10_off + 136 * ccomps * dcomps);

            auto g_z_0_xxz_zz = cbuffer.data(fd_geom_10_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxz_xx, g_z_0_xxz_xy, g_z_0_xxz_xz, g_z_0_xxz_yy, g_z_0_xxz_yz, g_z_0_xxz_zz, g_z_0_xz_xx, g_z_0_xz_xxx, g_z_0_xz_xxy, g_z_0_xz_xxz, g_z_0_xz_xy, g_z_0_xz_xyy, g_z_0_xz_xyz, g_z_0_xz_xz, g_z_0_xz_xzz, g_z_0_xz_yy, g_z_0_xz_yz, g_z_0_xz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxz_xx[k] = -g_z_0_xz_xx[k] * ab_x + g_z_0_xz_xxx[k];

                g_z_0_xxz_xy[k] = -g_z_0_xz_xy[k] * ab_x + g_z_0_xz_xxy[k];

                g_z_0_xxz_xz[k] = -g_z_0_xz_xz[k] * ab_x + g_z_0_xz_xxz[k];

                g_z_0_xxz_yy[k] = -g_z_0_xz_yy[k] * ab_x + g_z_0_xz_xyy[k];

                g_z_0_xxz_yz[k] = -g_z_0_xz_yz[k] * ab_x + g_z_0_xz_xyz[k];

                g_z_0_xxz_zz[k] = -g_z_0_xz_zz[k] * ab_x + g_z_0_xz_xzz[k];
            }

            /// Set up 138-144 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyy_xx = cbuffer.data(fd_geom_10_off + 138 * ccomps * dcomps);

            auto g_z_0_xyy_xy = cbuffer.data(fd_geom_10_off + 139 * ccomps * dcomps);

            auto g_z_0_xyy_xz = cbuffer.data(fd_geom_10_off + 140 * ccomps * dcomps);

            auto g_z_0_xyy_yy = cbuffer.data(fd_geom_10_off + 141 * ccomps * dcomps);

            auto g_z_0_xyy_yz = cbuffer.data(fd_geom_10_off + 142 * ccomps * dcomps);

            auto g_z_0_xyy_zz = cbuffer.data(fd_geom_10_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyy_xx, g_z_0_xyy_xy, g_z_0_xyy_xz, g_z_0_xyy_yy, g_z_0_xyy_yz, g_z_0_xyy_zz, g_z_0_yy_xx, g_z_0_yy_xxx, g_z_0_yy_xxy, g_z_0_yy_xxz, g_z_0_yy_xy, g_z_0_yy_xyy, g_z_0_yy_xyz, g_z_0_yy_xz, g_z_0_yy_xzz, g_z_0_yy_yy, g_z_0_yy_yz, g_z_0_yy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyy_xx[k] = -g_z_0_yy_xx[k] * ab_x + g_z_0_yy_xxx[k];

                g_z_0_xyy_xy[k] = -g_z_0_yy_xy[k] * ab_x + g_z_0_yy_xxy[k];

                g_z_0_xyy_xz[k] = -g_z_0_yy_xz[k] * ab_x + g_z_0_yy_xxz[k];

                g_z_0_xyy_yy[k] = -g_z_0_yy_yy[k] * ab_x + g_z_0_yy_xyy[k];

                g_z_0_xyy_yz[k] = -g_z_0_yy_yz[k] * ab_x + g_z_0_yy_xyz[k];

                g_z_0_xyy_zz[k] = -g_z_0_yy_zz[k] * ab_x + g_z_0_yy_xzz[k];
            }

            /// Set up 144-150 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyz_xx = cbuffer.data(fd_geom_10_off + 144 * ccomps * dcomps);

            auto g_z_0_xyz_xy = cbuffer.data(fd_geom_10_off + 145 * ccomps * dcomps);

            auto g_z_0_xyz_xz = cbuffer.data(fd_geom_10_off + 146 * ccomps * dcomps);

            auto g_z_0_xyz_yy = cbuffer.data(fd_geom_10_off + 147 * ccomps * dcomps);

            auto g_z_0_xyz_yz = cbuffer.data(fd_geom_10_off + 148 * ccomps * dcomps);

            auto g_z_0_xyz_zz = cbuffer.data(fd_geom_10_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyz_xx, g_z_0_xyz_xy, g_z_0_xyz_xz, g_z_0_xyz_yy, g_z_0_xyz_yz, g_z_0_xyz_zz, g_z_0_yz_xx, g_z_0_yz_xxx, g_z_0_yz_xxy, g_z_0_yz_xxz, g_z_0_yz_xy, g_z_0_yz_xyy, g_z_0_yz_xyz, g_z_0_yz_xz, g_z_0_yz_xzz, g_z_0_yz_yy, g_z_0_yz_yz, g_z_0_yz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyz_xx[k] = -g_z_0_yz_xx[k] * ab_x + g_z_0_yz_xxx[k];

                g_z_0_xyz_xy[k] = -g_z_0_yz_xy[k] * ab_x + g_z_0_yz_xxy[k];

                g_z_0_xyz_xz[k] = -g_z_0_yz_xz[k] * ab_x + g_z_0_yz_xxz[k];

                g_z_0_xyz_yy[k] = -g_z_0_yz_yy[k] * ab_x + g_z_0_yz_xyy[k];

                g_z_0_xyz_yz[k] = -g_z_0_yz_yz[k] * ab_x + g_z_0_yz_xyz[k];

                g_z_0_xyz_zz[k] = -g_z_0_yz_zz[k] * ab_x + g_z_0_yz_xzz[k];
            }

            /// Set up 150-156 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzz_xx = cbuffer.data(fd_geom_10_off + 150 * ccomps * dcomps);

            auto g_z_0_xzz_xy = cbuffer.data(fd_geom_10_off + 151 * ccomps * dcomps);

            auto g_z_0_xzz_xz = cbuffer.data(fd_geom_10_off + 152 * ccomps * dcomps);

            auto g_z_0_xzz_yy = cbuffer.data(fd_geom_10_off + 153 * ccomps * dcomps);

            auto g_z_0_xzz_yz = cbuffer.data(fd_geom_10_off + 154 * ccomps * dcomps);

            auto g_z_0_xzz_zz = cbuffer.data(fd_geom_10_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xzz_xx, g_z_0_xzz_xy, g_z_0_xzz_xz, g_z_0_xzz_yy, g_z_0_xzz_yz, g_z_0_xzz_zz, g_z_0_zz_xx, g_z_0_zz_xxx, g_z_0_zz_xxy, g_z_0_zz_xxz, g_z_0_zz_xy, g_z_0_zz_xyy, g_z_0_zz_xyz, g_z_0_zz_xz, g_z_0_zz_xzz, g_z_0_zz_yy, g_z_0_zz_yz, g_z_0_zz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzz_xx[k] = -g_z_0_zz_xx[k] * ab_x + g_z_0_zz_xxx[k];

                g_z_0_xzz_xy[k] = -g_z_0_zz_xy[k] * ab_x + g_z_0_zz_xxy[k];

                g_z_0_xzz_xz[k] = -g_z_0_zz_xz[k] * ab_x + g_z_0_zz_xxz[k];

                g_z_0_xzz_yy[k] = -g_z_0_zz_yy[k] * ab_x + g_z_0_zz_xyy[k];

                g_z_0_xzz_yz[k] = -g_z_0_zz_yz[k] * ab_x + g_z_0_zz_xyz[k];

                g_z_0_xzz_zz[k] = -g_z_0_zz_zz[k] * ab_x + g_z_0_zz_xzz[k];
            }

            /// Set up 156-162 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyy_xx = cbuffer.data(fd_geom_10_off + 156 * ccomps * dcomps);

            auto g_z_0_yyy_xy = cbuffer.data(fd_geom_10_off + 157 * ccomps * dcomps);

            auto g_z_0_yyy_xz = cbuffer.data(fd_geom_10_off + 158 * ccomps * dcomps);

            auto g_z_0_yyy_yy = cbuffer.data(fd_geom_10_off + 159 * ccomps * dcomps);

            auto g_z_0_yyy_yz = cbuffer.data(fd_geom_10_off + 160 * ccomps * dcomps);

            auto g_z_0_yyy_zz = cbuffer.data(fd_geom_10_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yy_xx, g_z_0_yy_xxy, g_z_0_yy_xy, g_z_0_yy_xyy, g_z_0_yy_xyz, g_z_0_yy_xz, g_z_0_yy_yy, g_z_0_yy_yyy, g_z_0_yy_yyz, g_z_0_yy_yz, g_z_0_yy_yzz, g_z_0_yy_zz, g_z_0_yyy_xx, g_z_0_yyy_xy, g_z_0_yyy_xz, g_z_0_yyy_yy, g_z_0_yyy_yz, g_z_0_yyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyy_xx[k] = -g_z_0_yy_xx[k] * ab_y + g_z_0_yy_xxy[k];

                g_z_0_yyy_xy[k] = -g_z_0_yy_xy[k] * ab_y + g_z_0_yy_xyy[k];

                g_z_0_yyy_xz[k] = -g_z_0_yy_xz[k] * ab_y + g_z_0_yy_xyz[k];

                g_z_0_yyy_yy[k] = -g_z_0_yy_yy[k] * ab_y + g_z_0_yy_yyy[k];

                g_z_0_yyy_yz[k] = -g_z_0_yy_yz[k] * ab_y + g_z_0_yy_yyz[k];

                g_z_0_yyy_zz[k] = -g_z_0_yy_zz[k] * ab_y + g_z_0_yy_yzz[k];
            }

            /// Set up 162-168 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyz_xx = cbuffer.data(fd_geom_10_off + 162 * ccomps * dcomps);

            auto g_z_0_yyz_xy = cbuffer.data(fd_geom_10_off + 163 * ccomps * dcomps);

            auto g_z_0_yyz_xz = cbuffer.data(fd_geom_10_off + 164 * ccomps * dcomps);

            auto g_z_0_yyz_yy = cbuffer.data(fd_geom_10_off + 165 * ccomps * dcomps);

            auto g_z_0_yyz_yz = cbuffer.data(fd_geom_10_off + 166 * ccomps * dcomps);

            auto g_z_0_yyz_zz = cbuffer.data(fd_geom_10_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyz_xx, g_z_0_yyz_xy, g_z_0_yyz_xz, g_z_0_yyz_yy, g_z_0_yyz_yz, g_z_0_yyz_zz, g_z_0_yz_xx, g_z_0_yz_xxy, g_z_0_yz_xy, g_z_0_yz_xyy, g_z_0_yz_xyz, g_z_0_yz_xz, g_z_0_yz_yy, g_z_0_yz_yyy, g_z_0_yz_yyz, g_z_0_yz_yz, g_z_0_yz_yzz, g_z_0_yz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyz_xx[k] = -g_z_0_yz_xx[k] * ab_y + g_z_0_yz_xxy[k];

                g_z_0_yyz_xy[k] = -g_z_0_yz_xy[k] * ab_y + g_z_0_yz_xyy[k];

                g_z_0_yyz_xz[k] = -g_z_0_yz_xz[k] * ab_y + g_z_0_yz_xyz[k];

                g_z_0_yyz_yy[k] = -g_z_0_yz_yy[k] * ab_y + g_z_0_yz_yyy[k];

                g_z_0_yyz_yz[k] = -g_z_0_yz_yz[k] * ab_y + g_z_0_yz_yyz[k];

                g_z_0_yyz_zz[k] = -g_z_0_yz_zz[k] * ab_y + g_z_0_yz_yzz[k];
            }

            /// Set up 168-174 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzz_xx = cbuffer.data(fd_geom_10_off + 168 * ccomps * dcomps);

            auto g_z_0_yzz_xy = cbuffer.data(fd_geom_10_off + 169 * ccomps * dcomps);

            auto g_z_0_yzz_xz = cbuffer.data(fd_geom_10_off + 170 * ccomps * dcomps);

            auto g_z_0_yzz_yy = cbuffer.data(fd_geom_10_off + 171 * ccomps * dcomps);

            auto g_z_0_yzz_yz = cbuffer.data(fd_geom_10_off + 172 * ccomps * dcomps);

            auto g_z_0_yzz_zz = cbuffer.data(fd_geom_10_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yzz_xx, g_z_0_yzz_xy, g_z_0_yzz_xz, g_z_0_yzz_yy, g_z_0_yzz_yz, g_z_0_yzz_zz, g_z_0_zz_xx, g_z_0_zz_xxy, g_z_0_zz_xy, g_z_0_zz_xyy, g_z_0_zz_xyz, g_z_0_zz_xz, g_z_0_zz_yy, g_z_0_zz_yyy, g_z_0_zz_yyz, g_z_0_zz_yz, g_z_0_zz_yzz, g_z_0_zz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzz_xx[k] = -g_z_0_zz_xx[k] * ab_y + g_z_0_zz_xxy[k];

                g_z_0_yzz_xy[k] = -g_z_0_zz_xy[k] * ab_y + g_z_0_zz_xyy[k];

                g_z_0_yzz_xz[k] = -g_z_0_zz_xz[k] * ab_y + g_z_0_zz_xyz[k];

                g_z_0_yzz_yy[k] = -g_z_0_zz_yy[k] * ab_y + g_z_0_zz_yyy[k];

                g_z_0_yzz_yz[k] = -g_z_0_zz_yz[k] * ab_y + g_z_0_zz_yyz[k];

                g_z_0_yzz_zz[k] = -g_z_0_zz_zz[k] * ab_y + g_z_0_zz_yzz[k];
            }

            /// Set up 174-180 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzz_xx = cbuffer.data(fd_geom_10_off + 174 * ccomps * dcomps);

            auto g_z_0_zzz_xy = cbuffer.data(fd_geom_10_off + 175 * ccomps * dcomps);

            auto g_z_0_zzz_xz = cbuffer.data(fd_geom_10_off + 176 * ccomps * dcomps);

            auto g_z_0_zzz_yy = cbuffer.data(fd_geom_10_off + 177 * ccomps * dcomps);

            auto g_z_0_zzz_yz = cbuffer.data(fd_geom_10_off + 178 * ccomps * dcomps);

            auto g_z_0_zzz_zz = cbuffer.data(fd_geom_10_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zz_xx, g_z_0_zz_xxz, g_z_0_zz_xy, g_z_0_zz_xyz, g_z_0_zz_xz, g_z_0_zz_xzz, g_z_0_zz_yy, g_z_0_zz_yyz, g_z_0_zz_yz, g_z_0_zz_yzz, g_z_0_zz_zz, g_z_0_zz_zzz, g_z_0_zzz_xx, g_z_0_zzz_xy, g_z_0_zzz_xz, g_z_0_zzz_yy, g_z_0_zzz_yz, g_z_0_zzz_zz, g_zz_xx, g_zz_xy, g_zz_xz, g_zz_yy, g_zz_yz, g_zz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzz_xx[k] = -g_zz_xx[k] - g_z_0_zz_xx[k] * ab_z + g_z_0_zz_xxz[k];

                g_z_0_zzz_xy[k] = -g_zz_xy[k] - g_z_0_zz_xy[k] * ab_z + g_z_0_zz_xyz[k];

                g_z_0_zzz_xz[k] = -g_zz_xz[k] - g_z_0_zz_xz[k] * ab_z + g_z_0_zz_xzz[k];

                g_z_0_zzz_yy[k] = -g_zz_yy[k] - g_z_0_zz_yy[k] * ab_z + g_z_0_zz_yyz[k];

                g_z_0_zzz_yz[k] = -g_zz_yz[k] - g_z_0_zz_yz[k] * ab_z + g_z_0_zz_yzz[k];

                g_z_0_zzz_zz[k] = -g_zz_zz[k] - g_z_0_zz_zz[k] * ab_z + g_z_0_zz_zzz[k];
            }
        }
    }
}

} // erirec namespace

