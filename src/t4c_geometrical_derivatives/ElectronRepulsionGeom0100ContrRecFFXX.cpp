#include "ElectronRepulsionGeom0100ContrRecFFXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_ffxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_ffxx,
                                            const size_t idx_dfxx,
                                            const size_t idx_geom_01_dfxx,
                                            const size_t idx_geom_01_dgxx,
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

            const auto df_off = idx_dfxx + i * dcomps + j;

            auto g_xx_xxx = cbuffer.data(df_off + 0 * ccomps * dcomps);

            auto g_xx_xxy = cbuffer.data(df_off + 1 * ccomps * dcomps);

            auto g_xx_xxz = cbuffer.data(df_off + 2 * ccomps * dcomps);

            auto g_xx_xyy = cbuffer.data(df_off + 3 * ccomps * dcomps);

            auto g_xx_xyz = cbuffer.data(df_off + 4 * ccomps * dcomps);

            auto g_xx_xzz = cbuffer.data(df_off + 5 * ccomps * dcomps);

            auto g_xx_yyy = cbuffer.data(df_off + 6 * ccomps * dcomps);

            auto g_xx_yyz = cbuffer.data(df_off + 7 * ccomps * dcomps);

            auto g_xx_yzz = cbuffer.data(df_off + 8 * ccomps * dcomps);

            auto g_xx_zzz = cbuffer.data(df_off + 9 * ccomps * dcomps);

            auto g_xy_xxx = cbuffer.data(df_off + 10 * ccomps * dcomps);

            auto g_xy_xxy = cbuffer.data(df_off + 11 * ccomps * dcomps);

            auto g_xy_xxz = cbuffer.data(df_off + 12 * ccomps * dcomps);

            auto g_xy_xyy = cbuffer.data(df_off + 13 * ccomps * dcomps);

            auto g_xy_xyz = cbuffer.data(df_off + 14 * ccomps * dcomps);

            auto g_xy_xzz = cbuffer.data(df_off + 15 * ccomps * dcomps);

            auto g_xy_yyy = cbuffer.data(df_off + 16 * ccomps * dcomps);

            auto g_xy_yyz = cbuffer.data(df_off + 17 * ccomps * dcomps);

            auto g_xy_yzz = cbuffer.data(df_off + 18 * ccomps * dcomps);

            auto g_xy_zzz = cbuffer.data(df_off + 19 * ccomps * dcomps);

            auto g_xz_xxx = cbuffer.data(df_off + 20 * ccomps * dcomps);

            auto g_xz_xxy = cbuffer.data(df_off + 21 * ccomps * dcomps);

            auto g_xz_xxz = cbuffer.data(df_off + 22 * ccomps * dcomps);

            auto g_xz_xyy = cbuffer.data(df_off + 23 * ccomps * dcomps);

            auto g_xz_xyz = cbuffer.data(df_off + 24 * ccomps * dcomps);

            auto g_xz_xzz = cbuffer.data(df_off + 25 * ccomps * dcomps);

            auto g_xz_yyy = cbuffer.data(df_off + 26 * ccomps * dcomps);

            auto g_xz_yyz = cbuffer.data(df_off + 27 * ccomps * dcomps);

            auto g_xz_yzz = cbuffer.data(df_off + 28 * ccomps * dcomps);

            auto g_xz_zzz = cbuffer.data(df_off + 29 * ccomps * dcomps);

            auto g_yy_xxx = cbuffer.data(df_off + 30 * ccomps * dcomps);

            auto g_yy_xxy = cbuffer.data(df_off + 31 * ccomps * dcomps);

            auto g_yy_xxz = cbuffer.data(df_off + 32 * ccomps * dcomps);

            auto g_yy_xyy = cbuffer.data(df_off + 33 * ccomps * dcomps);

            auto g_yy_xyz = cbuffer.data(df_off + 34 * ccomps * dcomps);

            auto g_yy_xzz = cbuffer.data(df_off + 35 * ccomps * dcomps);

            auto g_yy_yyy = cbuffer.data(df_off + 36 * ccomps * dcomps);

            auto g_yy_yyz = cbuffer.data(df_off + 37 * ccomps * dcomps);

            auto g_yy_yzz = cbuffer.data(df_off + 38 * ccomps * dcomps);

            auto g_yy_zzz = cbuffer.data(df_off + 39 * ccomps * dcomps);

            auto g_yz_xxx = cbuffer.data(df_off + 40 * ccomps * dcomps);

            auto g_yz_xxy = cbuffer.data(df_off + 41 * ccomps * dcomps);

            auto g_yz_xxz = cbuffer.data(df_off + 42 * ccomps * dcomps);

            auto g_yz_xyy = cbuffer.data(df_off + 43 * ccomps * dcomps);

            auto g_yz_xyz = cbuffer.data(df_off + 44 * ccomps * dcomps);

            auto g_yz_xzz = cbuffer.data(df_off + 45 * ccomps * dcomps);

            auto g_yz_yyy = cbuffer.data(df_off + 46 * ccomps * dcomps);

            auto g_yz_yyz = cbuffer.data(df_off + 47 * ccomps * dcomps);

            auto g_yz_yzz = cbuffer.data(df_off + 48 * ccomps * dcomps);

            auto g_yz_zzz = cbuffer.data(df_off + 49 * ccomps * dcomps);

            auto g_zz_xxx = cbuffer.data(df_off + 50 * ccomps * dcomps);

            auto g_zz_xxy = cbuffer.data(df_off + 51 * ccomps * dcomps);

            auto g_zz_xxz = cbuffer.data(df_off + 52 * ccomps * dcomps);

            auto g_zz_xyy = cbuffer.data(df_off + 53 * ccomps * dcomps);

            auto g_zz_xyz = cbuffer.data(df_off + 54 * ccomps * dcomps);

            auto g_zz_xzz = cbuffer.data(df_off + 55 * ccomps * dcomps);

            auto g_zz_yyy = cbuffer.data(df_off + 56 * ccomps * dcomps);

            auto g_zz_yyz = cbuffer.data(df_off + 57 * ccomps * dcomps);

            auto g_zz_yzz = cbuffer.data(df_off + 58 * ccomps * dcomps);

            auto g_zz_zzz = cbuffer.data(df_off + 59 * ccomps * dcomps);

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

            /// Set up components of auxilary buffer : DGSS

            const auto dg_geom_01_off = idx_geom_01_dgxx + i * dcomps + j;

            auto g_0_x_xx_xxxx = cbuffer.data(dg_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xx_xxxy = cbuffer.data(dg_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xx_xxxz = cbuffer.data(dg_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xx_xxyy = cbuffer.data(dg_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xx_xxyz = cbuffer.data(dg_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xx_xxzz = cbuffer.data(dg_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xx_xyyy = cbuffer.data(dg_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xx_xyyz = cbuffer.data(dg_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xx_xyzz = cbuffer.data(dg_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xx_xzzz = cbuffer.data(dg_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xx_yyyy = cbuffer.data(dg_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xx_yyyz = cbuffer.data(dg_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xx_yyzz = cbuffer.data(dg_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xx_yzzz = cbuffer.data(dg_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xx_zzzz = cbuffer.data(dg_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xy_xxxx = cbuffer.data(dg_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xy_xxxy = cbuffer.data(dg_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xy_xxxz = cbuffer.data(dg_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xy_xxyy = cbuffer.data(dg_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xy_xxyz = cbuffer.data(dg_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xy_xxzz = cbuffer.data(dg_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xy_xyyy = cbuffer.data(dg_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xy_xyyz = cbuffer.data(dg_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xy_xyzz = cbuffer.data(dg_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xy_xzzz = cbuffer.data(dg_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xy_yyyy = cbuffer.data(dg_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xy_yyyz = cbuffer.data(dg_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xy_yyzz = cbuffer.data(dg_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xy_yzzz = cbuffer.data(dg_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xy_zzzz = cbuffer.data(dg_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xz_xxxx = cbuffer.data(dg_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xz_xxxy = cbuffer.data(dg_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xz_xxxz = cbuffer.data(dg_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xz_xxyy = cbuffer.data(dg_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xz_xxyz = cbuffer.data(dg_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xz_xxzz = cbuffer.data(dg_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xz_xyyy = cbuffer.data(dg_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xz_xyyz = cbuffer.data(dg_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xz_xyzz = cbuffer.data(dg_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xz_xzzz = cbuffer.data(dg_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xz_yyyy = cbuffer.data(dg_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xz_yyyz = cbuffer.data(dg_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xz_yyzz = cbuffer.data(dg_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xz_yzzz = cbuffer.data(dg_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xz_zzzz = cbuffer.data(dg_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_yy_xxxx = cbuffer.data(dg_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_yy_xxxy = cbuffer.data(dg_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_yy_xxxz = cbuffer.data(dg_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_yy_xxyy = cbuffer.data(dg_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_yy_xxyz = cbuffer.data(dg_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_yy_xxzz = cbuffer.data(dg_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_yy_xyyy = cbuffer.data(dg_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_yy_xyyz = cbuffer.data(dg_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_yy_xyzz = cbuffer.data(dg_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_yy_xzzz = cbuffer.data(dg_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_yy_yyyy = cbuffer.data(dg_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_yy_yyyz = cbuffer.data(dg_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_yy_yyzz = cbuffer.data(dg_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_yy_yzzz = cbuffer.data(dg_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_yy_zzzz = cbuffer.data(dg_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_yz_xxxx = cbuffer.data(dg_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_yz_xxxy = cbuffer.data(dg_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_yz_xxxz = cbuffer.data(dg_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_yz_xxyy = cbuffer.data(dg_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_yz_xxyz = cbuffer.data(dg_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_yz_xxzz = cbuffer.data(dg_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_yz_xyyy = cbuffer.data(dg_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_yz_xyyz = cbuffer.data(dg_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_yz_xyzz = cbuffer.data(dg_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_yz_xzzz = cbuffer.data(dg_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_yz_yyyy = cbuffer.data(dg_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_yz_yyyz = cbuffer.data(dg_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_yz_yyzz = cbuffer.data(dg_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_yz_yzzz = cbuffer.data(dg_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_yz_zzzz = cbuffer.data(dg_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_zz_xxxx = cbuffer.data(dg_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_zz_xxxy = cbuffer.data(dg_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_zz_xxxz = cbuffer.data(dg_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_zz_xxyy = cbuffer.data(dg_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_zz_xxyz = cbuffer.data(dg_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_zz_xxzz = cbuffer.data(dg_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_zz_xyyy = cbuffer.data(dg_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_zz_xyyz = cbuffer.data(dg_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_zz_xyzz = cbuffer.data(dg_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_zz_xzzz = cbuffer.data(dg_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_zz_yyyy = cbuffer.data(dg_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_zz_yyyz = cbuffer.data(dg_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_zz_yyzz = cbuffer.data(dg_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_zz_yzzz = cbuffer.data(dg_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_zz_zzzz = cbuffer.data(dg_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_y_xx_xxxx = cbuffer.data(dg_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_y_xx_xxxy = cbuffer.data(dg_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_y_xx_xxxz = cbuffer.data(dg_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_y_xx_xxyy = cbuffer.data(dg_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_y_xx_xxyz = cbuffer.data(dg_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_y_xx_xxzz = cbuffer.data(dg_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_y_xx_xyyy = cbuffer.data(dg_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_y_xx_xyyz = cbuffer.data(dg_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_y_xx_xyzz = cbuffer.data(dg_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_y_xx_xzzz = cbuffer.data(dg_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_y_xx_yyyy = cbuffer.data(dg_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_y_xx_yyyz = cbuffer.data(dg_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_y_xx_yyzz = cbuffer.data(dg_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_y_xx_yzzz = cbuffer.data(dg_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_y_xx_zzzz = cbuffer.data(dg_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_y_xy_xxxx = cbuffer.data(dg_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_y_xy_xxxy = cbuffer.data(dg_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_y_xy_xxxz = cbuffer.data(dg_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_y_xy_xxyy = cbuffer.data(dg_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_xy_xxyz = cbuffer.data(dg_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_xy_xxzz = cbuffer.data(dg_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_y_xy_xyyy = cbuffer.data(dg_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_y_xy_xyyz = cbuffer.data(dg_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_xy_xyzz = cbuffer.data(dg_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_y_xy_xzzz = cbuffer.data(dg_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_xy_yyyy = cbuffer.data(dg_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_xy_yyyz = cbuffer.data(dg_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_y_xy_yyzz = cbuffer.data(dg_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_xy_yzzz = cbuffer.data(dg_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_xy_zzzz = cbuffer.data(dg_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_y_xz_xxxx = cbuffer.data(dg_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_y_xz_xxxy = cbuffer.data(dg_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_y_xz_xxxz = cbuffer.data(dg_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_y_xz_xxyy = cbuffer.data(dg_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_y_xz_xxyz = cbuffer.data(dg_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_y_xz_xxzz = cbuffer.data(dg_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_y_xz_xyyy = cbuffer.data(dg_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_y_xz_xyyz = cbuffer.data(dg_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_y_xz_xyzz = cbuffer.data(dg_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_y_xz_xzzz = cbuffer.data(dg_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_y_xz_yyyy = cbuffer.data(dg_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_y_xz_yyyz = cbuffer.data(dg_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_y_xz_yyzz = cbuffer.data(dg_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_y_xz_yzzz = cbuffer.data(dg_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_y_xz_zzzz = cbuffer.data(dg_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_y_yy_xxxx = cbuffer.data(dg_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_yy_xxxy = cbuffer.data(dg_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_yy_xxxz = cbuffer.data(dg_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_y_yy_xxyy = cbuffer.data(dg_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_yy_xxyz = cbuffer.data(dg_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_y_yy_xxzz = cbuffer.data(dg_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_y_yy_xyyy = cbuffer.data(dg_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_yy_xyyz = cbuffer.data(dg_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_yy_xyzz = cbuffer.data(dg_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_y_yy_xzzz = cbuffer.data(dg_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_yy_yyyy = cbuffer.data(dg_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_yy_yyyz = cbuffer.data(dg_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_y_yy_yyzz = cbuffer.data(dg_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_yy_yzzz = cbuffer.data(dg_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_yy_zzzz = cbuffer.data(dg_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_y_yz_xxxx = cbuffer.data(dg_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_yz_xxxy = cbuffer.data(dg_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_yz_xxxz = cbuffer.data(dg_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_yz_xxyy = cbuffer.data(dg_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_yz_xxyz = cbuffer.data(dg_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_yz_xxzz = cbuffer.data(dg_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_yz_xyyy = cbuffer.data(dg_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_yz_xyyz = cbuffer.data(dg_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_yz_xyzz = cbuffer.data(dg_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_yz_xzzz = cbuffer.data(dg_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_yz_yyyy = cbuffer.data(dg_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_yz_yyyz = cbuffer.data(dg_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_yz_yyzz = cbuffer.data(dg_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_yz_yzzz = cbuffer.data(dg_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_yz_zzzz = cbuffer.data(dg_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_zz_xxxx = cbuffer.data(dg_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_zz_xxxy = cbuffer.data(dg_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_zz_xxxz = cbuffer.data(dg_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_y_zz_xxyy = cbuffer.data(dg_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_zz_xxyz = cbuffer.data(dg_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_zz_xxzz = cbuffer.data(dg_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_zz_xyyy = cbuffer.data(dg_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_zz_xyyz = cbuffer.data(dg_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_zz_xyzz = cbuffer.data(dg_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_zz_xzzz = cbuffer.data(dg_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_zz_yyyy = cbuffer.data(dg_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_zz_yyyz = cbuffer.data(dg_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_zz_yyzz = cbuffer.data(dg_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_zz_yzzz = cbuffer.data(dg_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_zz_zzzz = cbuffer.data(dg_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_z_xx_xxxx = cbuffer.data(dg_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_z_xx_xxxy = cbuffer.data(dg_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_z_xx_xxxz = cbuffer.data(dg_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_z_xx_xxyy = cbuffer.data(dg_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_z_xx_xxyz = cbuffer.data(dg_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_z_xx_xxzz = cbuffer.data(dg_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_z_xx_xyyy = cbuffer.data(dg_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_z_xx_xyyz = cbuffer.data(dg_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_z_xx_xyzz = cbuffer.data(dg_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_z_xx_xzzz = cbuffer.data(dg_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_z_xx_yyyy = cbuffer.data(dg_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_z_xx_yyyz = cbuffer.data(dg_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_z_xx_yyzz = cbuffer.data(dg_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_z_xx_yzzz = cbuffer.data(dg_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_z_xx_zzzz = cbuffer.data(dg_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_z_xy_xxxx = cbuffer.data(dg_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_z_xy_xxxy = cbuffer.data(dg_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_z_xy_xxxz = cbuffer.data(dg_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_z_xy_xxyy = cbuffer.data(dg_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_z_xy_xxyz = cbuffer.data(dg_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_z_xy_xxzz = cbuffer.data(dg_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_z_xy_xyyy = cbuffer.data(dg_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_z_xy_xyyz = cbuffer.data(dg_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_z_xy_xyzz = cbuffer.data(dg_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_z_xy_xzzz = cbuffer.data(dg_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_z_xy_yyyy = cbuffer.data(dg_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_z_xy_yyyz = cbuffer.data(dg_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_z_xy_yyzz = cbuffer.data(dg_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_z_xy_yzzz = cbuffer.data(dg_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_z_xy_zzzz = cbuffer.data(dg_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_z_xz_xxxx = cbuffer.data(dg_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_z_xz_xxxy = cbuffer.data(dg_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_z_xz_xxxz = cbuffer.data(dg_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_z_xz_xxyy = cbuffer.data(dg_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_z_xz_xxyz = cbuffer.data(dg_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_z_xz_xxzz = cbuffer.data(dg_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_z_xz_xyyy = cbuffer.data(dg_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_z_xz_xyyz = cbuffer.data(dg_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_z_xz_xyzz = cbuffer.data(dg_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_z_xz_xzzz = cbuffer.data(dg_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_z_xz_yyyy = cbuffer.data(dg_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_z_xz_yyyz = cbuffer.data(dg_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_z_xz_yyzz = cbuffer.data(dg_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_z_xz_yzzz = cbuffer.data(dg_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_z_xz_zzzz = cbuffer.data(dg_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_z_yy_xxxx = cbuffer.data(dg_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_z_yy_xxxy = cbuffer.data(dg_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_z_yy_xxxz = cbuffer.data(dg_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_z_yy_xxyy = cbuffer.data(dg_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_z_yy_xxyz = cbuffer.data(dg_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_z_yy_xxzz = cbuffer.data(dg_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_z_yy_xyyy = cbuffer.data(dg_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_z_yy_xyyz = cbuffer.data(dg_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_z_yy_xyzz = cbuffer.data(dg_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_z_yy_xzzz = cbuffer.data(dg_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_z_yy_yyyy = cbuffer.data(dg_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_z_yy_yyyz = cbuffer.data(dg_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_z_yy_yyzz = cbuffer.data(dg_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_z_yy_yzzz = cbuffer.data(dg_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_z_yy_zzzz = cbuffer.data(dg_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_z_yz_xxxx = cbuffer.data(dg_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_z_yz_xxxy = cbuffer.data(dg_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_z_yz_xxxz = cbuffer.data(dg_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_z_yz_xxyy = cbuffer.data(dg_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_z_yz_xxyz = cbuffer.data(dg_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_z_yz_xxzz = cbuffer.data(dg_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_z_yz_xyyy = cbuffer.data(dg_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_z_yz_xyyz = cbuffer.data(dg_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_z_yz_xyzz = cbuffer.data(dg_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_z_yz_xzzz = cbuffer.data(dg_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_z_yz_yyyy = cbuffer.data(dg_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_z_yz_yyyz = cbuffer.data(dg_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_z_yz_yyzz = cbuffer.data(dg_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_z_yz_yzzz = cbuffer.data(dg_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_z_yz_zzzz = cbuffer.data(dg_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_z_zz_xxxx = cbuffer.data(dg_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_z_zz_xxxy = cbuffer.data(dg_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_z_zz_xxxz = cbuffer.data(dg_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_z_zz_xxyy = cbuffer.data(dg_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_z_zz_xxyz = cbuffer.data(dg_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_z_zz_xxzz = cbuffer.data(dg_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_z_zz_xyyy = cbuffer.data(dg_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_z_zz_xyyz = cbuffer.data(dg_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_z_zz_xyzz = cbuffer.data(dg_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_z_zz_xzzz = cbuffer.data(dg_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_z_zz_yyyy = cbuffer.data(dg_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_z_zz_yyyz = cbuffer.data(dg_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_z_zz_yyzz = cbuffer.data(dg_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_z_zz_yzzz = cbuffer.data(dg_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_z_zz_zzzz = cbuffer.data(dg_geom_01_off + 269 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_ffxx

            const auto ff_geom_01_off = idx_geom_01_ffxx + i * dcomps + j;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xx_xxx, g_0_x_xx_xxxx, g_0_x_xx_xxxy, g_0_x_xx_xxxz, g_0_x_xx_xxy, g_0_x_xx_xxyy, g_0_x_xx_xxyz, g_0_x_xx_xxz, g_0_x_xx_xxzz, g_0_x_xx_xyy, g_0_x_xx_xyyy, g_0_x_xx_xyyz, g_0_x_xx_xyz, g_0_x_xx_xyzz, g_0_x_xx_xzz, g_0_x_xx_xzzz, g_0_x_xx_yyy, g_0_x_xx_yyz, g_0_x_xx_yzz, g_0_x_xx_zzz, g_0_x_xxx_xxx, g_0_x_xxx_xxy, g_0_x_xxx_xxz, g_0_x_xxx_xyy, g_0_x_xxx_xyz, g_0_x_xxx_xzz, g_0_x_xxx_yyy, g_0_x_xxx_yyz, g_0_x_xxx_yzz, g_0_x_xxx_zzz, g_xx_xxx, g_xx_xxy, g_xx_xxz, g_xx_xyy, g_xx_xyz, g_xx_xzz, g_xx_yyy, g_xx_yyz, g_xx_yzz, g_xx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxx_xxx[k] = g_xx_xxx[k] - g_0_x_xx_xxx[k] * ab_x + g_0_x_xx_xxxx[k];

                g_0_x_xxx_xxy[k] = g_xx_xxy[k] - g_0_x_xx_xxy[k] * ab_x + g_0_x_xx_xxxy[k];

                g_0_x_xxx_xxz[k] = g_xx_xxz[k] - g_0_x_xx_xxz[k] * ab_x + g_0_x_xx_xxxz[k];

                g_0_x_xxx_xyy[k] = g_xx_xyy[k] - g_0_x_xx_xyy[k] * ab_x + g_0_x_xx_xxyy[k];

                g_0_x_xxx_xyz[k] = g_xx_xyz[k] - g_0_x_xx_xyz[k] * ab_x + g_0_x_xx_xxyz[k];

                g_0_x_xxx_xzz[k] = g_xx_xzz[k] - g_0_x_xx_xzz[k] * ab_x + g_0_x_xx_xxzz[k];

                g_0_x_xxx_yyy[k] = g_xx_yyy[k] - g_0_x_xx_yyy[k] * ab_x + g_0_x_xx_xyyy[k];

                g_0_x_xxx_yyz[k] = g_xx_yyz[k] - g_0_x_xx_yyz[k] * ab_x + g_0_x_xx_xyyz[k];

                g_0_x_xxx_yzz[k] = g_xx_yzz[k] - g_0_x_xx_yzz[k] * ab_x + g_0_x_xx_xyzz[k];

                g_0_x_xxx_zzz[k] = g_xx_zzz[k] - g_0_x_xx_zzz[k] * ab_x + g_0_x_xx_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xx_xxx, g_0_x_xx_xxxy, g_0_x_xx_xxy, g_0_x_xx_xxyy, g_0_x_xx_xxyz, g_0_x_xx_xxz, g_0_x_xx_xyy, g_0_x_xx_xyyy, g_0_x_xx_xyyz, g_0_x_xx_xyz, g_0_x_xx_xyzz, g_0_x_xx_xzz, g_0_x_xx_yyy, g_0_x_xx_yyyy, g_0_x_xx_yyyz, g_0_x_xx_yyz, g_0_x_xx_yyzz, g_0_x_xx_yzz, g_0_x_xx_yzzz, g_0_x_xx_zzz, g_0_x_xxy_xxx, g_0_x_xxy_xxy, g_0_x_xxy_xxz, g_0_x_xxy_xyy, g_0_x_xxy_xyz, g_0_x_xxy_xzz, g_0_x_xxy_yyy, g_0_x_xxy_yyz, g_0_x_xxy_yzz, g_0_x_xxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxy_xxx[k] = -g_0_x_xx_xxx[k] * ab_y + g_0_x_xx_xxxy[k];

                g_0_x_xxy_xxy[k] = -g_0_x_xx_xxy[k] * ab_y + g_0_x_xx_xxyy[k];

                g_0_x_xxy_xxz[k] = -g_0_x_xx_xxz[k] * ab_y + g_0_x_xx_xxyz[k];

                g_0_x_xxy_xyy[k] = -g_0_x_xx_xyy[k] * ab_y + g_0_x_xx_xyyy[k];

                g_0_x_xxy_xyz[k] = -g_0_x_xx_xyz[k] * ab_y + g_0_x_xx_xyyz[k];

                g_0_x_xxy_xzz[k] = -g_0_x_xx_xzz[k] * ab_y + g_0_x_xx_xyzz[k];

                g_0_x_xxy_yyy[k] = -g_0_x_xx_yyy[k] * ab_y + g_0_x_xx_yyyy[k];

                g_0_x_xxy_yyz[k] = -g_0_x_xx_yyz[k] * ab_y + g_0_x_xx_yyyz[k];

                g_0_x_xxy_yzz[k] = -g_0_x_xx_yzz[k] * ab_y + g_0_x_xx_yyzz[k];

                g_0_x_xxy_zzz[k] = -g_0_x_xx_zzz[k] * ab_y + g_0_x_xx_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xx_xxx, g_0_x_xx_xxxz, g_0_x_xx_xxy, g_0_x_xx_xxyz, g_0_x_xx_xxz, g_0_x_xx_xxzz, g_0_x_xx_xyy, g_0_x_xx_xyyz, g_0_x_xx_xyz, g_0_x_xx_xyzz, g_0_x_xx_xzz, g_0_x_xx_xzzz, g_0_x_xx_yyy, g_0_x_xx_yyyz, g_0_x_xx_yyz, g_0_x_xx_yyzz, g_0_x_xx_yzz, g_0_x_xx_yzzz, g_0_x_xx_zzz, g_0_x_xx_zzzz, g_0_x_xxz_xxx, g_0_x_xxz_xxy, g_0_x_xxz_xxz, g_0_x_xxz_xyy, g_0_x_xxz_xyz, g_0_x_xxz_xzz, g_0_x_xxz_yyy, g_0_x_xxz_yyz, g_0_x_xxz_yzz, g_0_x_xxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxz_xxx[k] = -g_0_x_xx_xxx[k] * ab_z + g_0_x_xx_xxxz[k];

                g_0_x_xxz_xxy[k] = -g_0_x_xx_xxy[k] * ab_z + g_0_x_xx_xxyz[k];

                g_0_x_xxz_xxz[k] = -g_0_x_xx_xxz[k] * ab_z + g_0_x_xx_xxzz[k];

                g_0_x_xxz_xyy[k] = -g_0_x_xx_xyy[k] * ab_z + g_0_x_xx_xyyz[k];

                g_0_x_xxz_xyz[k] = -g_0_x_xx_xyz[k] * ab_z + g_0_x_xx_xyzz[k];

                g_0_x_xxz_xzz[k] = -g_0_x_xx_xzz[k] * ab_z + g_0_x_xx_xzzz[k];

                g_0_x_xxz_yyy[k] = -g_0_x_xx_yyy[k] * ab_z + g_0_x_xx_yyyz[k];

                g_0_x_xxz_yyz[k] = -g_0_x_xx_yyz[k] * ab_z + g_0_x_xx_yyzz[k];

                g_0_x_xxz_yzz[k] = -g_0_x_xx_yzz[k] * ab_z + g_0_x_xx_yzzz[k];

                g_0_x_xxz_zzz[k] = -g_0_x_xx_zzz[k] * ab_z + g_0_x_xx_zzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xy_xxx, g_0_x_xy_xxxy, g_0_x_xy_xxy, g_0_x_xy_xxyy, g_0_x_xy_xxyz, g_0_x_xy_xxz, g_0_x_xy_xyy, g_0_x_xy_xyyy, g_0_x_xy_xyyz, g_0_x_xy_xyz, g_0_x_xy_xyzz, g_0_x_xy_xzz, g_0_x_xy_yyy, g_0_x_xy_yyyy, g_0_x_xy_yyyz, g_0_x_xy_yyz, g_0_x_xy_yyzz, g_0_x_xy_yzz, g_0_x_xy_yzzz, g_0_x_xy_zzz, g_0_x_xyy_xxx, g_0_x_xyy_xxy, g_0_x_xyy_xxz, g_0_x_xyy_xyy, g_0_x_xyy_xyz, g_0_x_xyy_xzz, g_0_x_xyy_yyy, g_0_x_xyy_yyz, g_0_x_xyy_yzz, g_0_x_xyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyy_xxx[k] = -g_0_x_xy_xxx[k] * ab_y + g_0_x_xy_xxxy[k];

                g_0_x_xyy_xxy[k] = -g_0_x_xy_xxy[k] * ab_y + g_0_x_xy_xxyy[k];

                g_0_x_xyy_xxz[k] = -g_0_x_xy_xxz[k] * ab_y + g_0_x_xy_xxyz[k];

                g_0_x_xyy_xyy[k] = -g_0_x_xy_xyy[k] * ab_y + g_0_x_xy_xyyy[k];

                g_0_x_xyy_xyz[k] = -g_0_x_xy_xyz[k] * ab_y + g_0_x_xy_xyyz[k];

                g_0_x_xyy_xzz[k] = -g_0_x_xy_xzz[k] * ab_y + g_0_x_xy_xyzz[k];

                g_0_x_xyy_yyy[k] = -g_0_x_xy_yyy[k] * ab_y + g_0_x_xy_yyyy[k];

                g_0_x_xyy_yyz[k] = -g_0_x_xy_yyz[k] * ab_y + g_0_x_xy_yyyz[k];

                g_0_x_xyy_yzz[k] = -g_0_x_xy_yzz[k] * ab_y + g_0_x_xy_yyzz[k];

                g_0_x_xyy_zzz[k] = -g_0_x_xy_zzz[k] * ab_y + g_0_x_xy_yzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xyz_xxx, g_0_x_xyz_xxy, g_0_x_xyz_xxz, g_0_x_xyz_xyy, g_0_x_xyz_xyz, g_0_x_xyz_xzz, g_0_x_xyz_yyy, g_0_x_xyz_yyz, g_0_x_xyz_yzz, g_0_x_xyz_zzz, g_0_x_xz_xxx, g_0_x_xz_xxxy, g_0_x_xz_xxy, g_0_x_xz_xxyy, g_0_x_xz_xxyz, g_0_x_xz_xxz, g_0_x_xz_xyy, g_0_x_xz_xyyy, g_0_x_xz_xyyz, g_0_x_xz_xyz, g_0_x_xz_xyzz, g_0_x_xz_xzz, g_0_x_xz_yyy, g_0_x_xz_yyyy, g_0_x_xz_yyyz, g_0_x_xz_yyz, g_0_x_xz_yyzz, g_0_x_xz_yzz, g_0_x_xz_yzzz, g_0_x_xz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyz_xxx[k] = -g_0_x_xz_xxx[k] * ab_y + g_0_x_xz_xxxy[k];

                g_0_x_xyz_xxy[k] = -g_0_x_xz_xxy[k] * ab_y + g_0_x_xz_xxyy[k];

                g_0_x_xyz_xxz[k] = -g_0_x_xz_xxz[k] * ab_y + g_0_x_xz_xxyz[k];

                g_0_x_xyz_xyy[k] = -g_0_x_xz_xyy[k] * ab_y + g_0_x_xz_xyyy[k];

                g_0_x_xyz_xyz[k] = -g_0_x_xz_xyz[k] * ab_y + g_0_x_xz_xyyz[k];

                g_0_x_xyz_xzz[k] = -g_0_x_xz_xzz[k] * ab_y + g_0_x_xz_xyzz[k];

                g_0_x_xyz_yyy[k] = -g_0_x_xz_yyy[k] * ab_y + g_0_x_xz_yyyy[k];

                g_0_x_xyz_yyz[k] = -g_0_x_xz_yyz[k] * ab_y + g_0_x_xz_yyyz[k];

                g_0_x_xyz_yzz[k] = -g_0_x_xz_yzz[k] * ab_y + g_0_x_xz_yyzz[k];

                g_0_x_xyz_zzz[k] = -g_0_x_xz_zzz[k] * ab_y + g_0_x_xz_yzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_xz_xxx, g_0_x_xz_xxxz, g_0_x_xz_xxy, g_0_x_xz_xxyz, g_0_x_xz_xxz, g_0_x_xz_xxzz, g_0_x_xz_xyy, g_0_x_xz_xyyz, g_0_x_xz_xyz, g_0_x_xz_xyzz, g_0_x_xz_xzz, g_0_x_xz_xzzz, g_0_x_xz_yyy, g_0_x_xz_yyyz, g_0_x_xz_yyz, g_0_x_xz_yyzz, g_0_x_xz_yzz, g_0_x_xz_yzzz, g_0_x_xz_zzz, g_0_x_xz_zzzz, g_0_x_xzz_xxx, g_0_x_xzz_xxy, g_0_x_xzz_xxz, g_0_x_xzz_xyy, g_0_x_xzz_xyz, g_0_x_xzz_xzz, g_0_x_xzz_yyy, g_0_x_xzz_yyz, g_0_x_xzz_yzz, g_0_x_xzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzz_xxx[k] = -g_0_x_xz_xxx[k] * ab_z + g_0_x_xz_xxxz[k];

                g_0_x_xzz_xxy[k] = -g_0_x_xz_xxy[k] * ab_z + g_0_x_xz_xxyz[k];

                g_0_x_xzz_xxz[k] = -g_0_x_xz_xxz[k] * ab_z + g_0_x_xz_xxzz[k];

                g_0_x_xzz_xyy[k] = -g_0_x_xz_xyy[k] * ab_z + g_0_x_xz_xyyz[k];

                g_0_x_xzz_xyz[k] = -g_0_x_xz_xyz[k] * ab_z + g_0_x_xz_xyzz[k];

                g_0_x_xzz_xzz[k] = -g_0_x_xz_xzz[k] * ab_z + g_0_x_xz_xzzz[k];

                g_0_x_xzz_yyy[k] = -g_0_x_xz_yyy[k] * ab_z + g_0_x_xz_yyyz[k];

                g_0_x_xzz_yyz[k] = -g_0_x_xz_yyz[k] * ab_z + g_0_x_xz_yyzz[k];

                g_0_x_xzz_yzz[k] = -g_0_x_xz_yzz[k] * ab_z + g_0_x_xz_yzzz[k];

                g_0_x_xzz_zzz[k] = -g_0_x_xz_zzz[k] * ab_z + g_0_x_xz_zzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_yy_xxx, g_0_x_yy_xxxy, g_0_x_yy_xxy, g_0_x_yy_xxyy, g_0_x_yy_xxyz, g_0_x_yy_xxz, g_0_x_yy_xyy, g_0_x_yy_xyyy, g_0_x_yy_xyyz, g_0_x_yy_xyz, g_0_x_yy_xyzz, g_0_x_yy_xzz, g_0_x_yy_yyy, g_0_x_yy_yyyy, g_0_x_yy_yyyz, g_0_x_yy_yyz, g_0_x_yy_yyzz, g_0_x_yy_yzz, g_0_x_yy_yzzz, g_0_x_yy_zzz, g_0_x_yyy_xxx, g_0_x_yyy_xxy, g_0_x_yyy_xxz, g_0_x_yyy_xyy, g_0_x_yyy_xyz, g_0_x_yyy_xzz, g_0_x_yyy_yyy, g_0_x_yyy_yyz, g_0_x_yyy_yzz, g_0_x_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyy_xxx[k] = -g_0_x_yy_xxx[k] * ab_y + g_0_x_yy_xxxy[k];

                g_0_x_yyy_xxy[k] = -g_0_x_yy_xxy[k] * ab_y + g_0_x_yy_xxyy[k];

                g_0_x_yyy_xxz[k] = -g_0_x_yy_xxz[k] * ab_y + g_0_x_yy_xxyz[k];

                g_0_x_yyy_xyy[k] = -g_0_x_yy_xyy[k] * ab_y + g_0_x_yy_xyyy[k];

                g_0_x_yyy_xyz[k] = -g_0_x_yy_xyz[k] * ab_y + g_0_x_yy_xyyz[k];

                g_0_x_yyy_xzz[k] = -g_0_x_yy_xzz[k] * ab_y + g_0_x_yy_xyzz[k];

                g_0_x_yyy_yyy[k] = -g_0_x_yy_yyy[k] * ab_y + g_0_x_yy_yyyy[k];

                g_0_x_yyy_yyz[k] = -g_0_x_yy_yyz[k] * ab_y + g_0_x_yy_yyyz[k];

                g_0_x_yyy_yzz[k] = -g_0_x_yy_yzz[k] * ab_y + g_0_x_yy_yyzz[k];

                g_0_x_yyy_zzz[k] = -g_0_x_yy_zzz[k] * ab_y + g_0_x_yy_yzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_yyz_xxx, g_0_x_yyz_xxy, g_0_x_yyz_xxz, g_0_x_yyz_xyy, g_0_x_yyz_xyz, g_0_x_yyz_xzz, g_0_x_yyz_yyy, g_0_x_yyz_yyz, g_0_x_yyz_yzz, g_0_x_yyz_zzz, g_0_x_yz_xxx, g_0_x_yz_xxxy, g_0_x_yz_xxy, g_0_x_yz_xxyy, g_0_x_yz_xxyz, g_0_x_yz_xxz, g_0_x_yz_xyy, g_0_x_yz_xyyy, g_0_x_yz_xyyz, g_0_x_yz_xyz, g_0_x_yz_xyzz, g_0_x_yz_xzz, g_0_x_yz_yyy, g_0_x_yz_yyyy, g_0_x_yz_yyyz, g_0_x_yz_yyz, g_0_x_yz_yyzz, g_0_x_yz_yzz, g_0_x_yz_yzzz, g_0_x_yz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyz_xxx[k] = -g_0_x_yz_xxx[k] * ab_y + g_0_x_yz_xxxy[k];

                g_0_x_yyz_xxy[k] = -g_0_x_yz_xxy[k] * ab_y + g_0_x_yz_xxyy[k];

                g_0_x_yyz_xxz[k] = -g_0_x_yz_xxz[k] * ab_y + g_0_x_yz_xxyz[k];

                g_0_x_yyz_xyy[k] = -g_0_x_yz_xyy[k] * ab_y + g_0_x_yz_xyyy[k];

                g_0_x_yyz_xyz[k] = -g_0_x_yz_xyz[k] * ab_y + g_0_x_yz_xyyz[k];

                g_0_x_yyz_xzz[k] = -g_0_x_yz_xzz[k] * ab_y + g_0_x_yz_xyzz[k];

                g_0_x_yyz_yyy[k] = -g_0_x_yz_yyy[k] * ab_y + g_0_x_yz_yyyy[k];

                g_0_x_yyz_yyz[k] = -g_0_x_yz_yyz[k] * ab_y + g_0_x_yz_yyyz[k];

                g_0_x_yyz_yzz[k] = -g_0_x_yz_yzz[k] * ab_y + g_0_x_yz_yyzz[k];

                g_0_x_yyz_zzz[k] = -g_0_x_yz_zzz[k] * ab_y + g_0_x_yz_yzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_yzz_xxx, g_0_x_yzz_xxy, g_0_x_yzz_xxz, g_0_x_yzz_xyy, g_0_x_yzz_xyz, g_0_x_yzz_xzz, g_0_x_yzz_yyy, g_0_x_yzz_yyz, g_0_x_yzz_yzz, g_0_x_yzz_zzz, g_0_x_zz_xxx, g_0_x_zz_xxxy, g_0_x_zz_xxy, g_0_x_zz_xxyy, g_0_x_zz_xxyz, g_0_x_zz_xxz, g_0_x_zz_xyy, g_0_x_zz_xyyy, g_0_x_zz_xyyz, g_0_x_zz_xyz, g_0_x_zz_xyzz, g_0_x_zz_xzz, g_0_x_zz_yyy, g_0_x_zz_yyyy, g_0_x_zz_yyyz, g_0_x_zz_yyz, g_0_x_zz_yyzz, g_0_x_zz_yzz, g_0_x_zz_yzzz, g_0_x_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzz_xxx[k] = -g_0_x_zz_xxx[k] * ab_y + g_0_x_zz_xxxy[k];

                g_0_x_yzz_xxy[k] = -g_0_x_zz_xxy[k] * ab_y + g_0_x_zz_xxyy[k];

                g_0_x_yzz_xxz[k] = -g_0_x_zz_xxz[k] * ab_y + g_0_x_zz_xxyz[k];

                g_0_x_yzz_xyy[k] = -g_0_x_zz_xyy[k] * ab_y + g_0_x_zz_xyyy[k];

                g_0_x_yzz_xyz[k] = -g_0_x_zz_xyz[k] * ab_y + g_0_x_zz_xyyz[k];

                g_0_x_yzz_xzz[k] = -g_0_x_zz_xzz[k] * ab_y + g_0_x_zz_xyzz[k];

                g_0_x_yzz_yyy[k] = -g_0_x_zz_yyy[k] * ab_y + g_0_x_zz_yyyy[k];

                g_0_x_yzz_yyz[k] = -g_0_x_zz_yyz[k] * ab_y + g_0_x_zz_yyyz[k];

                g_0_x_yzz_yzz[k] = -g_0_x_zz_yzz[k] * ab_y + g_0_x_zz_yyzz[k];

                g_0_x_yzz_zzz[k] = -g_0_x_zz_zzz[k] * ab_y + g_0_x_zz_yzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_zz_xxx, g_0_x_zz_xxxz, g_0_x_zz_xxy, g_0_x_zz_xxyz, g_0_x_zz_xxz, g_0_x_zz_xxzz, g_0_x_zz_xyy, g_0_x_zz_xyyz, g_0_x_zz_xyz, g_0_x_zz_xyzz, g_0_x_zz_xzz, g_0_x_zz_xzzz, g_0_x_zz_yyy, g_0_x_zz_yyyz, g_0_x_zz_yyz, g_0_x_zz_yyzz, g_0_x_zz_yzz, g_0_x_zz_yzzz, g_0_x_zz_zzz, g_0_x_zz_zzzz, g_0_x_zzz_xxx, g_0_x_zzz_xxy, g_0_x_zzz_xxz, g_0_x_zzz_xyy, g_0_x_zzz_xyz, g_0_x_zzz_xzz, g_0_x_zzz_yyy, g_0_x_zzz_yyz, g_0_x_zzz_yzz, g_0_x_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzz_xxx[k] = -g_0_x_zz_xxx[k] * ab_z + g_0_x_zz_xxxz[k];

                g_0_x_zzz_xxy[k] = -g_0_x_zz_xxy[k] * ab_z + g_0_x_zz_xxyz[k];

                g_0_x_zzz_xxz[k] = -g_0_x_zz_xxz[k] * ab_z + g_0_x_zz_xxzz[k];

                g_0_x_zzz_xyy[k] = -g_0_x_zz_xyy[k] * ab_z + g_0_x_zz_xyyz[k];

                g_0_x_zzz_xyz[k] = -g_0_x_zz_xyz[k] * ab_z + g_0_x_zz_xyzz[k];

                g_0_x_zzz_xzz[k] = -g_0_x_zz_xzz[k] * ab_z + g_0_x_zz_xzzz[k];

                g_0_x_zzz_yyy[k] = -g_0_x_zz_yyy[k] * ab_z + g_0_x_zz_yyyz[k];

                g_0_x_zzz_yyz[k] = -g_0_x_zz_yyz[k] * ab_z + g_0_x_zz_yyzz[k];

                g_0_x_zzz_yzz[k] = -g_0_x_zz_yzz[k] * ab_z + g_0_x_zz_yzzz[k];

                g_0_x_zzz_zzz[k] = -g_0_x_zz_zzz[k] * ab_z + g_0_x_zz_zzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xx_xxx, g_0_y_xx_xxxx, g_0_y_xx_xxxy, g_0_y_xx_xxxz, g_0_y_xx_xxy, g_0_y_xx_xxyy, g_0_y_xx_xxyz, g_0_y_xx_xxz, g_0_y_xx_xxzz, g_0_y_xx_xyy, g_0_y_xx_xyyy, g_0_y_xx_xyyz, g_0_y_xx_xyz, g_0_y_xx_xyzz, g_0_y_xx_xzz, g_0_y_xx_xzzz, g_0_y_xx_yyy, g_0_y_xx_yyz, g_0_y_xx_yzz, g_0_y_xx_zzz, g_0_y_xxx_xxx, g_0_y_xxx_xxy, g_0_y_xxx_xxz, g_0_y_xxx_xyy, g_0_y_xxx_xyz, g_0_y_xxx_xzz, g_0_y_xxx_yyy, g_0_y_xxx_yyz, g_0_y_xxx_yzz, g_0_y_xxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxx_xxx[k] = -g_0_y_xx_xxx[k] * ab_x + g_0_y_xx_xxxx[k];

                g_0_y_xxx_xxy[k] = -g_0_y_xx_xxy[k] * ab_x + g_0_y_xx_xxxy[k];

                g_0_y_xxx_xxz[k] = -g_0_y_xx_xxz[k] * ab_x + g_0_y_xx_xxxz[k];

                g_0_y_xxx_xyy[k] = -g_0_y_xx_xyy[k] * ab_x + g_0_y_xx_xxyy[k];

                g_0_y_xxx_xyz[k] = -g_0_y_xx_xyz[k] * ab_x + g_0_y_xx_xxyz[k];

                g_0_y_xxx_xzz[k] = -g_0_y_xx_xzz[k] * ab_x + g_0_y_xx_xxzz[k];

                g_0_y_xxx_yyy[k] = -g_0_y_xx_yyy[k] * ab_x + g_0_y_xx_xyyy[k];

                g_0_y_xxx_yyz[k] = -g_0_y_xx_yyz[k] * ab_x + g_0_y_xx_xyyz[k];

                g_0_y_xxx_yzz[k] = -g_0_y_xx_yzz[k] * ab_x + g_0_y_xx_xyzz[k];

                g_0_y_xxx_zzz[k] = -g_0_y_xx_zzz[k] * ab_x + g_0_y_xx_xzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xxy_xxx, g_0_y_xxy_xxy, g_0_y_xxy_xxz, g_0_y_xxy_xyy, g_0_y_xxy_xyz, g_0_y_xxy_xzz, g_0_y_xxy_yyy, g_0_y_xxy_yyz, g_0_y_xxy_yzz, g_0_y_xxy_zzz, g_0_y_xy_xxx, g_0_y_xy_xxxx, g_0_y_xy_xxxy, g_0_y_xy_xxxz, g_0_y_xy_xxy, g_0_y_xy_xxyy, g_0_y_xy_xxyz, g_0_y_xy_xxz, g_0_y_xy_xxzz, g_0_y_xy_xyy, g_0_y_xy_xyyy, g_0_y_xy_xyyz, g_0_y_xy_xyz, g_0_y_xy_xyzz, g_0_y_xy_xzz, g_0_y_xy_xzzz, g_0_y_xy_yyy, g_0_y_xy_yyz, g_0_y_xy_yzz, g_0_y_xy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxy_xxx[k] = -g_0_y_xy_xxx[k] * ab_x + g_0_y_xy_xxxx[k];

                g_0_y_xxy_xxy[k] = -g_0_y_xy_xxy[k] * ab_x + g_0_y_xy_xxxy[k];

                g_0_y_xxy_xxz[k] = -g_0_y_xy_xxz[k] * ab_x + g_0_y_xy_xxxz[k];

                g_0_y_xxy_xyy[k] = -g_0_y_xy_xyy[k] * ab_x + g_0_y_xy_xxyy[k];

                g_0_y_xxy_xyz[k] = -g_0_y_xy_xyz[k] * ab_x + g_0_y_xy_xxyz[k];

                g_0_y_xxy_xzz[k] = -g_0_y_xy_xzz[k] * ab_x + g_0_y_xy_xxzz[k];

                g_0_y_xxy_yyy[k] = -g_0_y_xy_yyy[k] * ab_x + g_0_y_xy_xyyy[k];

                g_0_y_xxy_yyz[k] = -g_0_y_xy_yyz[k] * ab_x + g_0_y_xy_xyyz[k];

                g_0_y_xxy_yzz[k] = -g_0_y_xy_yzz[k] * ab_x + g_0_y_xy_xyzz[k];

                g_0_y_xxy_zzz[k] = -g_0_y_xy_zzz[k] * ab_x + g_0_y_xy_xzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xxz_xxx, g_0_y_xxz_xxy, g_0_y_xxz_xxz, g_0_y_xxz_xyy, g_0_y_xxz_xyz, g_0_y_xxz_xzz, g_0_y_xxz_yyy, g_0_y_xxz_yyz, g_0_y_xxz_yzz, g_0_y_xxz_zzz, g_0_y_xz_xxx, g_0_y_xz_xxxx, g_0_y_xz_xxxy, g_0_y_xz_xxxz, g_0_y_xz_xxy, g_0_y_xz_xxyy, g_0_y_xz_xxyz, g_0_y_xz_xxz, g_0_y_xz_xxzz, g_0_y_xz_xyy, g_0_y_xz_xyyy, g_0_y_xz_xyyz, g_0_y_xz_xyz, g_0_y_xz_xyzz, g_0_y_xz_xzz, g_0_y_xz_xzzz, g_0_y_xz_yyy, g_0_y_xz_yyz, g_0_y_xz_yzz, g_0_y_xz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxz_xxx[k] = -g_0_y_xz_xxx[k] * ab_x + g_0_y_xz_xxxx[k];

                g_0_y_xxz_xxy[k] = -g_0_y_xz_xxy[k] * ab_x + g_0_y_xz_xxxy[k];

                g_0_y_xxz_xxz[k] = -g_0_y_xz_xxz[k] * ab_x + g_0_y_xz_xxxz[k];

                g_0_y_xxz_xyy[k] = -g_0_y_xz_xyy[k] * ab_x + g_0_y_xz_xxyy[k];

                g_0_y_xxz_xyz[k] = -g_0_y_xz_xyz[k] * ab_x + g_0_y_xz_xxyz[k];

                g_0_y_xxz_xzz[k] = -g_0_y_xz_xzz[k] * ab_x + g_0_y_xz_xxzz[k];

                g_0_y_xxz_yyy[k] = -g_0_y_xz_yyy[k] * ab_x + g_0_y_xz_xyyy[k];

                g_0_y_xxz_yyz[k] = -g_0_y_xz_yyz[k] * ab_x + g_0_y_xz_xyyz[k];

                g_0_y_xxz_yzz[k] = -g_0_y_xz_yzz[k] * ab_x + g_0_y_xz_xyzz[k];

                g_0_y_xxz_zzz[k] = -g_0_y_xz_zzz[k] * ab_x + g_0_y_xz_xzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xyy_xxx, g_0_y_xyy_xxy, g_0_y_xyy_xxz, g_0_y_xyy_xyy, g_0_y_xyy_xyz, g_0_y_xyy_xzz, g_0_y_xyy_yyy, g_0_y_xyy_yyz, g_0_y_xyy_yzz, g_0_y_xyy_zzz, g_0_y_yy_xxx, g_0_y_yy_xxxx, g_0_y_yy_xxxy, g_0_y_yy_xxxz, g_0_y_yy_xxy, g_0_y_yy_xxyy, g_0_y_yy_xxyz, g_0_y_yy_xxz, g_0_y_yy_xxzz, g_0_y_yy_xyy, g_0_y_yy_xyyy, g_0_y_yy_xyyz, g_0_y_yy_xyz, g_0_y_yy_xyzz, g_0_y_yy_xzz, g_0_y_yy_xzzz, g_0_y_yy_yyy, g_0_y_yy_yyz, g_0_y_yy_yzz, g_0_y_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyy_xxx[k] = -g_0_y_yy_xxx[k] * ab_x + g_0_y_yy_xxxx[k];

                g_0_y_xyy_xxy[k] = -g_0_y_yy_xxy[k] * ab_x + g_0_y_yy_xxxy[k];

                g_0_y_xyy_xxz[k] = -g_0_y_yy_xxz[k] * ab_x + g_0_y_yy_xxxz[k];

                g_0_y_xyy_xyy[k] = -g_0_y_yy_xyy[k] * ab_x + g_0_y_yy_xxyy[k];

                g_0_y_xyy_xyz[k] = -g_0_y_yy_xyz[k] * ab_x + g_0_y_yy_xxyz[k];

                g_0_y_xyy_xzz[k] = -g_0_y_yy_xzz[k] * ab_x + g_0_y_yy_xxzz[k];

                g_0_y_xyy_yyy[k] = -g_0_y_yy_yyy[k] * ab_x + g_0_y_yy_xyyy[k];

                g_0_y_xyy_yyz[k] = -g_0_y_yy_yyz[k] * ab_x + g_0_y_yy_xyyz[k];

                g_0_y_xyy_yzz[k] = -g_0_y_yy_yzz[k] * ab_x + g_0_y_yy_xyzz[k];

                g_0_y_xyy_zzz[k] = -g_0_y_yy_zzz[k] * ab_x + g_0_y_yy_xzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xyz_xxx, g_0_y_xyz_xxy, g_0_y_xyz_xxz, g_0_y_xyz_xyy, g_0_y_xyz_xyz, g_0_y_xyz_xzz, g_0_y_xyz_yyy, g_0_y_xyz_yyz, g_0_y_xyz_yzz, g_0_y_xyz_zzz, g_0_y_yz_xxx, g_0_y_yz_xxxx, g_0_y_yz_xxxy, g_0_y_yz_xxxz, g_0_y_yz_xxy, g_0_y_yz_xxyy, g_0_y_yz_xxyz, g_0_y_yz_xxz, g_0_y_yz_xxzz, g_0_y_yz_xyy, g_0_y_yz_xyyy, g_0_y_yz_xyyz, g_0_y_yz_xyz, g_0_y_yz_xyzz, g_0_y_yz_xzz, g_0_y_yz_xzzz, g_0_y_yz_yyy, g_0_y_yz_yyz, g_0_y_yz_yzz, g_0_y_yz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyz_xxx[k] = -g_0_y_yz_xxx[k] * ab_x + g_0_y_yz_xxxx[k];

                g_0_y_xyz_xxy[k] = -g_0_y_yz_xxy[k] * ab_x + g_0_y_yz_xxxy[k];

                g_0_y_xyz_xxz[k] = -g_0_y_yz_xxz[k] * ab_x + g_0_y_yz_xxxz[k];

                g_0_y_xyz_xyy[k] = -g_0_y_yz_xyy[k] * ab_x + g_0_y_yz_xxyy[k];

                g_0_y_xyz_xyz[k] = -g_0_y_yz_xyz[k] * ab_x + g_0_y_yz_xxyz[k];

                g_0_y_xyz_xzz[k] = -g_0_y_yz_xzz[k] * ab_x + g_0_y_yz_xxzz[k];

                g_0_y_xyz_yyy[k] = -g_0_y_yz_yyy[k] * ab_x + g_0_y_yz_xyyy[k];

                g_0_y_xyz_yyz[k] = -g_0_y_yz_yyz[k] * ab_x + g_0_y_yz_xyyz[k];

                g_0_y_xyz_yzz[k] = -g_0_y_yz_yzz[k] * ab_x + g_0_y_yz_xyzz[k];

                g_0_y_xyz_zzz[k] = -g_0_y_yz_zzz[k] * ab_x + g_0_y_yz_xzzz[k];
            }

            /// Set up 150-160 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_xzz_xxx, g_0_y_xzz_xxy, g_0_y_xzz_xxz, g_0_y_xzz_xyy, g_0_y_xzz_xyz, g_0_y_xzz_xzz, g_0_y_xzz_yyy, g_0_y_xzz_yyz, g_0_y_xzz_yzz, g_0_y_xzz_zzz, g_0_y_zz_xxx, g_0_y_zz_xxxx, g_0_y_zz_xxxy, g_0_y_zz_xxxz, g_0_y_zz_xxy, g_0_y_zz_xxyy, g_0_y_zz_xxyz, g_0_y_zz_xxz, g_0_y_zz_xxzz, g_0_y_zz_xyy, g_0_y_zz_xyyy, g_0_y_zz_xyyz, g_0_y_zz_xyz, g_0_y_zz_xyzz, g_0_y_zz_xzz, g_0_y_zz_xzzz, g_0_y_zz_yyy, g_0_y_zz_yyz, g_0_y_zz_yzz, g_0_y_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzz_xxx[k] = -g_0_y_zz_xxx[k] * ab_x + g_0_y_zz_xxxx[k];

                g_0_y_xzz_xxy[k] = -g_0_y_zz_xxy[k] * ab_x + g_0_y_zz_xxxy[k];

                g_0_y_xzz_xxz[k] = -g_0_y_zz_xxz[k] * ab_x + g_0_y_zz_xxxz[k];

                g_0_y_xzz_xyy[k] = -g_0_y_zz_xyy[k] * ab_x + g_0_y_zz_xxyy[k];

                g_0_y_xzz_xyz[k] = -g_0_y_zz_xyz[k] * ab_x + g_0_y_zz_xxyz[k];

                g_0_y_xzz_xzz[k] = -g_0_y_zz_xzz[k] * ab_x + g_0_y_zz_xxzz[k];

                g_0_y_xzz_yyy[k] = -g_0_y_zz_yyy[k] * ab_x + g_0_y_zz_xyyy[k];

                g_0_y_xzz_yyz[k] = -g_0_y_zz_yyz[k] * ab_x + g_0_y_zz_xyyz[k];

                g_0_y_xzz_yzz[k] = -g_0_y_zz_yzz[k] * ab_x + g_0_y_zz_xyzz[k];

                g_0_y_xzz_zzz[k] = -g_0_y_zz_zzz[k] * ab_x + g_0_y_zz_xzzz[k];
            }

            /// Set up 160-170 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_yy_xxx, g_0_y_yy_xxxy, g_0_y_yy_xxy, g_0_y_yy_xxyy, g_0_y_yy_xxyz, g_0_y_yy_xxz, g_0_y_yy_xyy, g_0_y_yy_xyyy, g_0_y_yy_xyyz, g_0_y_yy_xyz, g_0_y_yy_xyzz, g_0_y_yy_xzz, g_0_y_yy_yyy, g_0_y_yy_yyyy, g_0_y_yy_yyyz, g_0_y_yy_yyz, g_0_y_yy_yyzz, g_0_y_yy_yzz, g_0_y_yy_yzzz, g_0_y_yy_zzz, g_0_y_yyy_xxx, g_0_y_yyy_xxy, g_0_y_yyy_xxz, g_0_y_yyy_xyy, g_0_y_yyy_xyz, g_0_y_yyy_xzz, g_0_y_yyy_yyy, g_0_y_yyy_yyz, g_0_y_yyy_yzz, g_0_y_yyy_zzz, g_yy_xxx, g_yy_xxy, g_yy_xxz, g_yy_xyy, g_yy_xyz, g_yy_xzz, g_yy_yyy, g_yy_yyz, g_yy_yzz, g_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyy_xxx[k] = g_yy_xxx[k] - g_0_y_yy_xxx[k] * ab_y + g_0_y_yy_xxxy[k];

                g_0_y_yyy_xxy[k] = g_yy_xxy[k] - g_0_y_yy_xxy[k] * ab_y + g_0_y_yy_xxyy[k];

                g_0_y_yyy_xxz[k] = g_yy_xxz[k] - g_0_y_yy_xxz[k] * ab_y + g_0_y_yy_xxyz[k];

                g_0_y_yyy_xyy[k] = g_yy_xyy[k] - g_0_y_yy_xyy[k] * ab_y + g_0_y_yy_xyyy[k];

                g_0_y_yyy_xyz[k] = g_yy_xyz[k] - g_0_y_yy_xyz[k] * ab_y + g_0_y_yy_xyyz[k];

                g_0_y_yyy_xzz[k] = g_yy_xzz[k] - g_0_y_yy_xzz[k] * ab_y + g_0_y_yy_xyzz[k];

                g_0_y_yyy_yyy[k] = g_yy_yyy[k] - g_0_y_yy_yyy[k] * ab_y + g_0_y_yy_yyyy[k];

                g_0_y_yyy_yyz[k] = g_yy_yyz[k] - g_0_y_yy_yyz[k] * ab_y + g_0_y_yy_yyyz[k];

                g_0_y_yyy_yzz[k] = g_yy_yzz[k] - g_0_y_yy_yzz[k] * ab_y + g_0_y_yy_yyzz[k];

                g_0_y_yyy_zzz[k] = g_yy_zzz[k] - g_0_y_yy_zzz[k] * ab_y + g_0_y_yy_yzzz[k];
            }

            /// Set up 170-180 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_yy_xxx, g_0_y_yy_xxxz, g_0_y_yy_xxy, g_0_y_yy_xxyz, g_0_y_yy_xxz, g_0_y_yy_xxzz, g_0_y_yy_xyy, g_0_y_yy_xyyz, g_0_y_yy_xyz, g_0_y_yy_xyzz, g_0_y_yy_xzz, g_0_y_yy_xzzz, g_0_y_yy_yyy, g_0_y_yy_yyyz, g_0_y_yy_yyz, g_0_y_yy_yyzz, g_0_y_yy_yzz, g_0_y_yy_yzzz, g_0_y_yy_zzz, g_0_y_yy_zzzz, g_0_y_yyz_xxx, g_0_y_yyz_xxy, g_0_y_yyz_xxz, g_0_y_yyz_xyy, g_0_y_yyz_xyz, g_0_y_yyz_xzz, g_0_y_yyz_yyy, g_0_y_yyz_yyz, g_0_y_yyz_yzz, g_0_y_yyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyz_xxx[k] = -g_0_y_yy_xxx[k] * ab_z + g_0_y_yy_xxxz[k];

                g_0_y_yyz_xxy[k] = -g_0_y_yy_xxy[k] * ab_z + g_0_y_yy_xxyz[k];

                g_0_y_yyz_xxz[k] = -g_0_y_yy_xxz[k] * ab_z + g_0_y_yy_xxzz[k];

                g_0_y_yyz_xyy[k] = -g_0_y_yy_xyy[k] * ab_z + g_0_y_yy_xyyz[k];

                g_0_y_yyz_xyz[k] = -g_0_y_yy_xyz[k] * ab_z + g_0_y_yy_xyzz[k];

                g_0_y_yyz_xzz[k] = -g_0_y_yy_xzz[k] * ab_z + g_0_y_yy_xzzz[k];

                g_0_y_yyz_yyy[k] = -g_0_y_yy_yyy[k] * ab_z + g_0_y_yy_yyyz[k];

                g_0_y_yyz_yyz[k] = -g_0_y_yy_yyz[k] * ab_z + g_0_y_yy_yyzz[k];

                g_0_y_yyz_yzz[k] = -g_0_y_yy_yzz[k] * ab_z + g_0_y_yy_yzzz[k];

                g_0_y_yyz_zzz[k] = -g_0_y_yy_zzz[k] * ab_z + g_0_y_yy_zzzz[k];
            }

            /// Set up 180-190 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_yz_xxx, g_0_y_yz_xxxz, g_0_y_yz_xxy, g_0_y_yz_xxyz, g_0_y_yz_xxz, g_0_y_yz_xxzz, g_0_y_yz_xyy, g_0_y_yz_xyyz, g_0_y_yz_xyz, g_0_y_yz_xyzz, g_0_y_yz_xzz, g_0_y_yz_xzzz, g_0_y_yz_yyy, g_0_y_yz_yyyz, g_0_y_yz_yyz, g_0_y_yz_yyzz, g_0_y_yz_yzz, g_0_y_yz_yzzz, g_0_y_yz_zzz, g_0_y_yz_zzzz, g_0_y_yzz_xxx, g_0_y_yzz_xxy, g_0_y_yzz_xxz, g_0_y_yzz_xyy, g_0_y_yzz_xyz, g_0_y_yzz_xzz, g_0_y_yzz_yyy, g_0_y_yzz_yyz, g_0_y_yzz_yzz, g_0_y_yzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzz_xxx[k] = -g_0_y_yz_xxx[k] * ab_z + g_0_y_yz_xxxz[k];

                g_0_y_yzz_xxy[k] = -g_0_y_yz_xxy[k] * ab_z + g_0_y_yz_xxyz[k];

                g_0_y_yzz_xxz[k] = -g_0_y_yz_xxz[k] * ab_z + g_0_y_yz_xxzz[k];

                g_0_y_yzz_xyy[k] = -g_0_y_yz_xyy[k] * ab_z + g_0_y_yz_xyyz[k];

                g_0_y_yzz_xyz[k] = -g_0_y_yz_xyz[k] * ab_z + g_0_y_yz_xyzz[k];

                g_0_y_yzz_xzz[k] = -g_0_y_yz_xzz[k] * ab_z + g_0_y_yz_xzzz[k];

                g_0_y_yzz_yyy[k] = -g_0_y_yz_yyy[k] * ab_z + g_0_y_yz_yyyz[k];

                g_0_y_yzz_yyz[k] = -g_0_y_yz_yyz[k] * ab_z + g_0_y_yz_yyzz[k];

                g_0_y_yzz_yzz[k] = -g_0_y_yz_yzz[k] * ab_z + g_0_y_yz_yzzz[k];

                g_0_y_yzz_zzz[k] = -g_0_y_yz_zzz[k] * ab_z + g_0_y_yz_zzzz[k];
            }

            /// Set up 190-200 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_y_zz_xxx, g_0_y_zz_xxxz, g_0_y_zz_xxy, g_0_y_zz_xxyz, g_0_y_zz_xxz, g_0_y_zz_xxzz, g_0_y_zz_xyy, g_0_y_zz_xyyz, g_0_y_zz_xyz, g_0_y_zz_xyzz, g_0_y_zz_xzz, g_0_y_zz_xzzz, g_0_y_zz_yyy, g_0_y_zz_yyyz, g_0_y_zz_yyz, g_0_y_zz_yyzz, g_0_y_zz_yzz, g_0_y_zz_yzzz, g_0_y_zz_zzz, g_0_y_zz_zzzz, g_0_y_zzz_xxx, g_0_y_zzz_xxy, g_0_y_zzz_xxz, g_0_y_zzz_xyy, g_0_y_zzz_xyz, g_0_y_zzz_xzz, g_0_y_zzz_yyy, g_0_y_zzz_yyz, g_0_y_zzz_yzz, g_0_y_zzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzz_xxx[k] = -g_0_y_zz_xxx[k] * ab_z + g_0_y_zz_xxxz[k];

                g_0_y_zzz_xxy[k] = -g_0_y_zz_xxy[k] * ab_z + g_0_y_zz_xxyz[k];

                g_0_y_zzz_xxz[k] = -g_0_y_zz_xxz[k] * ab_z + g_0_y_zz_xxzz[k];

                g_0_y_zzz_xyy[k] = -g_0_y_zz_xyy[k] * ab_z + g_0_y_zz_xyyz[k];

                g_0_y_zzz_xyz[k] = -g_0_y_zz_xyz[k] * ab_z + g_0_y_zz_xyzz[k];

                g_0_y_zzz_xzz[k] = -g_0_y_zz_xzz[k] * ab_z + g_0_y_zz_xzzz[k];

                g_0_y_zzz_yyy[k] = -g_0_y_zz_yyy[k] * ab_z + g_0_y_zz_yyyz[k];

                g_0_y_zzz_yyz[k] = -g_0_y_zz_yyz[k] * ab_z + g_0_y_zz_yyzz[k];

                g_0_y_zzz_yzz[k] = -g_0_y_zz_yzz[k] * ab_z + g_0_y_zz_yzzz[k];

                g_0_y_zzz_zzz[k] = -g_0_y_zz_zzz[k] * ab_z + g_0_y_zz_zzzz[k];
            }

            /// Set up 200-210 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xx_xxx, g_0_z_xx_xxxx, g_0_z_xx_xxxy, g_0_z_xx_xxxz, g_0_z_xx_xxy, g_0_z_xx_xxyy, g_0_z_xx_xxyz, g_0_z_xx_xxz, g_0_z_xx_xxzz, g_0_z_xx_xyy, g_0_z_xx_xyyy, g_0_z_xx_xyyz, g_0_z_xx_xyz, g_0_z_xx_xyzz, g_0_z_xx_xzz, g_0_z_xx_xzzz, g_0_z_xx_yyy, g_0_z_xx_yyz, g_0_z_xx_yzz, g_0_z_xx_zzz, g_0_z_xxx_xxx, g_0_z_xxx_xxy, g_0_z_xxx_xxz, g_0_z_xxx_xyy, g_0_z_xxx_xyz, g_0_z_xxx_xzz, g_0_z_xxx_yyy, g_0_z_xxx_yyz, g_0_z_xxx_yzz, g_0_z_xxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxx_xxx[k] = -g_0_z_xx_xxx[k] * ab_x + g_0_z_xx_xxxx[k];

                g_0_z_xxx_xxy[k] = -g_0_z_xx_xxy[k] * ab_x + g_0_z_xx_xxxy[k];

                g_0_z_xxx_xxz[k] = -g_0_z_xx_xxz[k] * ab_x + g_0_z_xx_xxxz[k];

                g_0_z_xxx_xyy[k] = -g_0_z_xx_xyy[k] * ab_x + g_0_z_xx_xxyy[k];

                g_0_z_xxx_xyz[k] = -g_0_z_xx_xyz[k] * ab_x + g_0_z_xx_xxyz[k];

                g_0_z_xxx_xzz[k] = -g_0_z_xx_xzz[k] * ab_x + g_0_z_xx_xxzz[k];

                g_0_z_xxx_yyy[k] = -g_0_z_xx_yyy[k] * ab_x + g_0_z_xx_xyyy[k];

                g_0_z_xxx_yyz[k] = -g_0_z_xx_yyz[k] * ab_x + g_0_z_xx_xyyz[k];

                g_0_z_xxx_yzz[k] = -g_0_z_xx_yzz[k] * ab_x + g_0_z_xx_xyzz[k];

                g_0_z_xxx_zzz[k] = -g_0_z_xx_zzz[k] * ab_x + g_0_z_xx_xzzz[k];
            }

            /// Set up 210-220 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xxy_xxx, g_0_z_xxy_xxy, g_0_z_xxy_xxz, g_0_z_xxy_xyy, g_0_z_xxy_xyz, g_0_z_xxy_xzz, g_0_z_xxy_yyy, g_0_z_xxy_yyz, g_0_z_xxy_yzz, g_0_z_xxy_zzz, g_0_z_xy_xxx, g_0_z_xy_xxxx, g_0_z_xy_xxxy, g_0_z_xy_xxxz, g_0_z_xy_xxy, g_0_z_xy_xxyy, g_0_z_xy_xxyz, g_0_z_xy_xxz, g_0_z_xy_xxzz, g_0_z_xy_xyy, g_0_z_xy_xyyy, g_0_z_xy_xyyz, g_0_z_xy_xyz, g_0_z_xy_xyzz, g_0_z_xy_xzz, g_0_z_xy_xzzz, g_0_z_xy_yyy, g_0_z_xy_yyz, g_0_z_xy_yzz, g_0_z_xy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxy_xxx[k] = -g_0_z_xy_xxx[k] * ab_x + g_0_z_xy_xxxx[k];

                g_0_z_xxy_xxy[k] = -g_0_z_xy_xxy[k] * ab_x + g_0_z_xy_xxxy[k];

                g_0_z_xxy_xxz[k] = -g_0_z_xy_xxz[k] * ab_x + g_0_z_xy_xxxz[k];

                g_0_z_xxy_xyy[k] = -g_0_z_xy_xyy[k] * ab_x + g_0_z_xy_xxyy[k];

                g_0_z_xxy_xyz[k] = -g_0_z_xy_xyz[k] * ab_x + g_0_z_xy_xxyz[k];

                g_0_z_xxy_xzz[k] = -g_0_z_xy_xzz[k] * ab_x + g_0_z_xy_xxzz[k];

                g_0_z_xxy_yyy[k] = -g_0_z_xy_yyy[k] * ab_x + g_0_z_xy_xyyy[k];

                g_0_z_xxy_yyz[k] = -g_0_z_xy_yyz[k] * ab_x + g_0_z_xy_xyyz[k];

                g_0_z_xxy_yzz[k] = -g_0_z_xy_yzz[k] * ab_x + g_0_z_xy_xyzz[k];

                g_0_z_xxy_zzz[k] = -g_0_z_xy_zzz[k] * ab_x + g_0_z_xy_xzzz[k];
            }

            /// Set up 220-230 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xxz_xxx, g_0_z_xxz_xxy, g_0_z_xxz_xxz, g_0_z_xxz_xyy, g_0_z_xxz_xyz, g_0_z_xxz_xzz, g_0_z_xxz_yyy, g_0_z_xxz_yyz, g_0_z_xxz_yzz, g_0_z_xxz_zzz, g_0_z_xz_xxx, g_0_z_xz_xxxx, g_0_z_xz_xxxy, g_0_z_xz_xxxz, g_0_z_xz_xxy, g_0_z_xz_xxyy, g_0_z_xz_xxyz, g_0_z_xz_xxz, g_0_z_xz_xxzz, g_0_z_xz_xyy, g_0_z_xz_xyyy, g_0_z_xz_xyyz, g_0_z_xz_xyz, g_0_z_xz_xyzz, g_0_z_xz_xzz, g_0_z_xz_xzzz, g_0_z_xz_yyy, g_0_z_xz_yyz, g_0_z_xz_yzz, g_0_z_xz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxz_xxx[k] = -g_0_z_xz_xxx[k] * ab_x + g_0_z_xz_xxxx[k];

                g_0_z_xxz_xxy[k] = -g_0_z_xz_xxy[k] * ab_x + g_0_z_xz_xxxy[k];

                g_0_z_xxz_xxz[k] = -g_0_z_xz_xxz[k] * ab_x + g_0_z_xz_xxxz[k];

                g_0_z_xxz_xyy[k] = -g_0_z_xz_xyy[k] * ab_x + g_0_z_xz_xxyy[k];

                g_0_z_xxz_xyz[k] = -g_0_z_xz_xyz[k] * ab_x + g_0_z_xz_xxyz[k];

                g_0_z_xxz_xzz[k] = -g_0_z_xz_xzz[k] * ab_x + g_0_z_xz_xxzz[k];

                g_0_z_xxz_yyy[k] = -g_0_z_xz_yyy[k] * ab_x + g_0_z_xz_xyyy[k];

                g_0_z_xxz_yyz[k] = -g_0_z_xz_yyz[k] * ab_x + g_0_z_xz_xyyz[k];

                g_0_z_xxz_yzz[k] = -g_0_z_xz_yzz[k] * ab_x + g_0_z_xz_xyzz[k];

                g_0_z_xxz_zzz[k] = -g_0_z_xz_zzz[k] * ab_x + g_0_z_xz_xzzz[k];
            }

            /// Set up 230-240 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xyy_xxx, g_0_z_xyy_xxy, g_0_z_xyy_xxz, g_0_z_xyy_xyy, g_0_z_xyy_xyz, g_0_z_xyy_xzz, g_0_z_xyy_yyy, g_0_z_xyy_yyz, g_0_z_xyy_yzz, g_0_z_xyy_zzz, g_0_z_yy_xxx, g_0_z_yy_xxxx, g_0_z_yy_xxxy, g_0_z_yy_xxxz, g_0_z_yy_xxy, g_0_z_yy_xxyy, g_0_z_yy_xxyz, g_0_z_yy_xxz, g_0_z_yy_xxzz, g_0_z_yy_xyy, g_0_z_yy_xyyy, g_0_z_yy_xyyz, g_0_z_yy_xyz, g_0_z_yy_xyzz, g_0_z_yy_xzz, g_0_z_yy_xzzz, g_0_z_yy_yyy, g_0_z_yy_yyz, g_0_z_yy_yzz, g_0_z_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyy_xxx[k] = -g_0_z_yy_xxx[k] * ab_x + g_0_z_yy_xxxx[k];

                g_0_z_xyy_xxy[k] = -g_0_z_yy_xxy[k] * ab_x + g_0_z_yy_xxxy[k];

                g_0_z_xyy_xxz[k] = -g_0_z_yy_xxz[k] * ab_x + g_0_z_yy_xxxz[k];

                g_0_z_xyy_xyy[k] = -g_0_z_yy_xyy[k] * ab_x + g_0_z_yy_xxyy[k];

                g_0_z_xyy_xyz[k] = -g_0_z_yy_xyz[k] * ab_x + g_0_z_yy_xxyz[k];

                g_0_z_xyy_xzz[k] = -g_0_z_yy_xzz[k] * ab_x + g_0_z_yy_xxzz[k];

                g_0_z_xyy_yyy[k] = -g_0_z_yy_yyy[k] * ab_x + g_0_z_yy_xyyy[k];

                g_0_z_xyy_yyz[k] = -g_0_z_yy_yyz[k] * ab_x + g_0_z_yy_xyyz[k];

                g_0_z_xyy_yzz[k] = -g_0_z_yy_yzz[k] * ab_x + g_0_z_yy_xyzz[k];

                g_0_z_xyy_zzz[k] = -g_0_z_yy_zzz[k] * ab_x + g_0_z_yy_xzzz[k];
            }

            /// Set up 240-250 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xyz_xxx, g_0_z_xyz_xxy, g_0_z_xyz_xxz, g_0_z_xyz_xyy, g_0_z_xyz_xyz, g_0_z_xyz_xzz, g_0_z_xyz_yyy, g_0_z_xyz_yyz, g_0_z_xyz_yzz, g_0_z_xyz_zzz, g_0_z_yz_xxx, g_0_z_yz_xxxx, g_0_z_yz_xxxy, g_0_z_yz_xxxz, g_0_z_yz_xxy, g_0_z_yz_xxyy, g_0_z_yz_xxyz, g_0_z_yz_xxz, g_0_z_yz_xxzz, g_0_z_yz_xyy, g_0_z_yz_xyyy, g_0_z_yz_xyyz, g_0_z_yz_xyz, g_0_z_yz_xyzz, g_0_z_yz_xzz, g_0_z_yz_xzzz, g_0_z_yz_yyy, g_0_z_yz_yyz, g_0_z_yz_yzz, g_0_z_yz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyz_xxx[k] = -g_0_z_yz_xxx[k] * ab_x + g_0_z_yz_xxxx[k];

                g_0_z_xyz_xxy[k] = -g_0_z_yz_xxy[k] * ab_x + g_0_z_yz_xxxy[k];

                g_0_z_xyz_xxz[k] = -g_0_z_yz_xxz[k] * ab_x + g_0_z_yz_xxxz[k];

                g_0_z_xyz_xyy[k] = -g_0_z_yz_xyy[k] * ab_x + g_0_z_yz_xxyy[k];

                g_0_z_xyz_xyz[k] = -g_0_z_yz_xyz[k] * ab_x + g_0_z_yz_xxyz[k];

                g_0_z_xyz_xzz[k] = -g_0_z_yz_xzz[k] * ab_x + g_0_z_yz_xxzz[k];

                g_0_z_xyz_yyy[k] = -g_0_z_yz_yyy[k] * ab_x + g_0_z_yz_xyyy[k];

                g_0_z_xyz_yyz[k] = -g_0_z_yz_yyz[k] * ab_x + g_0_z_yz_xyyz[k];

                g_0_z_xyz_yzz[k] = -g_0_z_yz_yzz[k] * ab_x + g_0_z_yz_xyzz[k];

                g_0_z_xyz_zzz[k] = -g_0_z_yz_zzz[k] * ab_x + g_0_z_yz_xzzz[k];
            }

            /// Set up 250-260 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_xzz_xxx, g_0_z_xzz_xxy, g_0_z_xzz_xxz, g_0_z_xzz_xyy, g_0_z_xzz_xyz, g_0_z_xzz_xzz, g_0_z_xzz_yyy, g_0_z_xzz_yyz, g_0_z_xzz_yzz, g_0_z_xzz_zzz, g_0_z_zz_xxx, g_0_z_zz_xxxx, g_0_z_zz_xxxy, g_0_z_zz_xxxz, g_0_z_zz_xxy, g_0_z_zz_xxyy, g_0_z_zz_xxyz, g_0_z_zz_xxz, g_0_z_zz_xxzz, g_0_z_zz_xyy, g_0_z_zz_xyyy, g_0_z_zz_xyyz, g_0_z_zz_xyz, g_0_z_zz_xyzz, g_0_z_zz_xzz, g_0_z_zz_xzzz, g_0_z_zz_yyy, g_0_z_zz_yyz, g_0_z_zz_yzz, g_0_z_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzz_xxx[k] = -g_0_z_zz_xxx[k] * ab_x + g_0_z_zz_xxxx[k];

                g_0_z_xzz_xxy[k] = -g_0_z_zz_xxy[k] * ab_x + g_0_z_zz_xxxy[k];

                g_0_z_xzz_xxz[k] = -g_0_z_zz_xxz[k] * ab_x + g_0_z_zz_xxxz[k];

                g_0_z_xzz_xyy[k] = -g_0_z_zz_xyy[k] * ab_x + g_0_z_zz_xxyy[k];

                g_0_z_xzz_xyz[k] = -g_0_z_zz_xyz[k] * ab_x + g_0_z_zz_xxyz[k];

                g_0_z_xzz_xzz[k] = -g_0_z_zz_xzz[k] * ab_x + g_0_z_zz_xxzz[k];

                g_0_z_xzz_yyy[k] = -g_0_z_zz_yyy[k] * ab_x + g_0_z_zz_xyyy[k];

                g_0_z_xzz_yyz[k] = -g_0_z_zz_yyz[k] * ab_x + g_0_z_zz_xyyz[k];

                g_0_z_xzz_yzz[k] = -g_0_z_zz_yzz[k] * ab_x + g_0_z_zz_xyzz[k];

                g_0_z_xzz_zzz[k] = -g_0_z_zz_zzz[k] * ab_x + g_0_z_zz_xzzz[k];
            }

            /// Set up 260-270 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_yy_xxx, g_0_z_yy_xxxy, g_0_z_yy_xxy, g_0_z_yy_xxyy, g_0_z_yy_xxyz, g_0_z_yy_xxz, g_0_z_yy_xyy, g_0_z_yy_xyyy, g_0_z_yy_xyyz, g_0_z_yy_xyz, g_0_z_yy_xyzz, g_0_z_yy_xzz, g_0_z_yy_yyy, g_0_z_yy_yyyy, g_0_z_yy_yyyz, g_0_z_yy_yyz, g_0_z_yy_yyzz, g_0_z_yy_yzz, g_0_z_yy_yzzz, g_0_z_yy_zzz, g_0_z_yyy_xxx, g_0_z_yyy_xxy, g_0_z_yyy_xxz, g_0_z_yyy_xyy, g_0_z_yyy_xyz, g_0_z_yyy_xzz, g_0_z_yyy_yyy, g_0_z_yyy_yyz, g_0_z_yyy_yzz, g_0_z_yyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyy_xxx[k] = -g_0_z_yy_xxx[k] * ab_y + g_0_z_yy_xxxy[k];

                g_0_z_yyy_xxy[k] = -g_0_z_yy_xxy[k] * ab_y + g_0_z_yy_xxyy[k];

                g_0_z_yyy_xxz[k] = -g_0_z_yy_xxz[k] * ab_y + g_0_z_yy_xxyz[k];

                g_0_z_yyy_xyy[k] = -g_0_z_yy_xyy[k] * ab_y + g_0_z_yy_xyyy[k];

                g_0_z_yyy_xyz[k] = -g_0_z_yy_xyz[k] * ab_y + g_0_z_yy_xyyz[k];

                g_0_z_yyy_xzz[k] = -g_0_z_yy_xzz[k] * ab_y + g_0_z_yy_xyzz[k];

                g_0_z_yyy_yyy[k] = -g_0_z_yy_yyy[k] * ab_y + g_0_z_yy_yyyy[k];

                g_0_z_yyy_yyz[k] = -g_0_z_yy_yyz[k] * ab_y + g_0_z_yy_yyyz[k];

                g_0_z_yyy_yzz[k] = -g_0_z_yy_yzz[k] * ab_y + g_0_z_yy_yyzz[k];

                g_0_z_yyy_zzz[k] = -g_0_z_yy_zzz[k] * ab_y + g_0_z_yy_yzzz[k];
            }

            /// Set up 270-280 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_yyz_xxx, g_0_z_yyz_xxy, g_0_z_yyz_xxz, g_0_z_yyz_xyy, g_0_z_yyz_xyz, g_0_z_yyz_xzz, g_0_z_yyz_yyy, g_0_z_yyz_yyz, g_0_z_yyz_yzz, g_0_z_yyz_zzz, g_0_z_yz_xxx, g_0_z_yz_xxxy, g_0_z_yz_xxy, g_0_z_yz_xxyy, g_0_z_yz_xxyz, g_0_z_yz_xxz, g_0_z_yz_xyy, g_0_z_yz_xyyy, g_0_z_yz_xyyz, g_0_z_yz_xyz, g_0_z_yz_xyzz, g_0_z_yz_xzz, g_0_z_yz_yyy, g_0_z_yz_yyyy, g_0_z_yz_yyyz, g_0_z_yz_yyz, g_0_z_yz_yyzz, g_0_z_yz_yzz, g_0_z_yz_yzzz, g_0_z_yz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyz_xxx[k] = -g_0_z_yz_xxx[k] * ab_y + g_0_z_yz_xxxy[k];

                g_0_z_yyz_xxy[k] = -g_0_z_yz_xxy[k] * ab_y + g_0_z_yz_xxyy[k];

                g_0_z_yyz_xxz[k] = -g_0_z_yz_xxz[k] * ab_y + g_0_z_yz_xxyz[k];

                g_0_z_yyz_xyy[k] = -g_0_z_yz_xyy[k] * ab_y + g_0_z_yz_xyyy[k];

                g_0_z_yyz_xyz[k] = -g_0_z_yz_xyz[k] * ab_y + g_0_z_yz_xyyz[k];

                g_0_z_yyz_xzz[k] = -g_0_z_yz_xzz[k] * ab_y + g_0_z_yz_xyzz[k];

                g_0_z_yyz_yyy[k] = -g_0_z_yz_yyy[k] * ab_y + g_0_z_yz_yyyy[k];

                g_0_z_yyz_yyz[k] = -g_0_z_yz_yyz[k] * ab_y + g_0_z_yz_yyyz[k];

                g_0_z_yyz_yzz[k] = -g_0_z_yz_yzz[k] * ab_y + g_0_z_yz_yyzz[k];

                g_0_z_yyz_zzz[k] = -g_0_z_yz_zzz[k] * ab_y + g_0_z_yz_yzzz[k];
            }

            /// Set up 280-290 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_yzz_xxx, g_0_z_yzz_xxy, g_0_z_yzz_xxz, g_0_z_yzz_xyy, g_0_z_yzz_xyz, g_0_z_yzz_xzz, g_0_z_yzz_yyy, g_0_z_yzz_yyz, g_0_z_yzz_yzz, g_0_z_yzz_zzz, g_0_z_zz_xxx, g_0_z_zz_xxxy, g_0_z_zz_xxy, g_0_z_zz_xxyy, g_0_z_zz_xxyz, g_0_z_zz_xxz, g_0_z_zz_xyy, g_0_z_zz_xyyy, g_0_z_zz_xyyz, g_0_z_zz_xyz, g_0_z_zz_xyzz, g_0_z_zz_xzz, g_0_z_zz_yyy, g_0_z_zz_yyyy, g_0_z_zz_yyyz, g_0_z_zz_yyz, g_0_z_zz_yyzz, g_0_z_zz_yzz, g_0_z_zz_yzzz, g_0_z_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzz_xxx[k] = -g_0_z_zz_xxx[k] * ab_y + g_0_z_zz_xxxy[k];

                g_0_z_yzz_xxy[k] = -g_0_z_zz_xxy[k] * ab_y + g_0_z_zz_xxyy[k];

                g_0_z_yzz_xxz[k] = -g_0_z_zz_xxz[k] * ab_y + g_0_z_zz_xxyz[k];

                g_0_z_yzz_xyy[k] = -g_0_z_zz_xyy[k] * ab_y + g_0_z_zz_xyyy[k];

                g_0_z_yzz_xyz[k] = -g_0_z_zz_xyz[k] * ab_y + g_0_z_zz_xyyz[k];

                g_0_z_yzz_xzz[k] = -g_0_z_zz_xzz[k] * ab_y + g_0_z_zz_xyzz[k];

                g_0_z_yzz_yyy[k] = -g_0_z_zz_yyy[k] * ab_y + g_0_z_zz_yyyy[k];

                g_0_z_yzz_yyz[k] = -g_0_z_zz_yyz[k] * ab_y + g_0_z_zz_yyyz[k];

                g_0_z_yzz_yzz[k] = -g_0_z_zz_yzz[k] * ab_y + g_0_z_zz_yyzz[k];

                g_0_z_yzz_zzz[k] = -g_0_z_zz_zzz[k] * ab_y + g_0_z_zz_yzzz[k];
            }

            /// Set up 290-300 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_z_zz_xxx, g_0_z_zz_xxxz, g_0_z_zz_xxy, g_0_z_zz_xxyz, g_0_z_zz_xxz, g_0_z_zz_xxzz, g_0_z_zz_xyy, g_0_z_zz_xyyz, g_0_z_zz_xyz, g_0_z_zz_xyzz, g_0_z_zz_xzz, g_0_z_zz_xzzz, g_0_z_zz_yyy, g_0_z_zz_yyyz, g_0_z_zz_yyz, g_0_z_zz_yyzz, g_0_z_zz_yzz, g_0_z_zz_yzzz, g_0_z_zz_zzz, g_0_z_zz_zzzz, g_0_z_zzz_xxx, g_0_z_zzz_xxy, g_0_z_zzz_xxz, g_0_z_zzz_xyy, g_0_z_zzz_xyz, g_0_z_zzz_xzz, g_0_z_zzz_yyy, g_0_z_zzz_yyz, g_0_z_zzz_yzz, g_0_z_zzz_zzz, g_zz_xxx, g_zz_xxy, g_zz_xxz, g_zz_xyy, g_zz_xyz, g_zz_xzz, g_zz_yyy, g_zz_yyz, g_zz_yzz, g_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzz_xxx[k] = g_zz_xxx[k] - g_0_z_zz_xxx[k] * ab_z + g_0_z_zz_xxxz[k];

                g_0_z_zzz_xxy[k] = g_zz_xxy[k] - g_0_z_zz_xxy[k] * ab_z + g_0_z_zz_xxyz[k];

                g_0_z_zzz_xxz[k] = g_zz_xxz[k] - g_0_z_zz_xxz[k] * ab_z + g_0_z_zz_xxzz[k];

                g_0_z_zzz_xyy[k] = g_zz_xyy[k] - g_0_z_zz_xyy[k] * ab_z + g_0_z_zz_xyyz[k];

                g_0_z_zzz_xyz[k] = g_zz_xyz[k] - g_0_z_zz_xyz[k] * ab_z + g_0_z_zz_xyzz[k];

                g_0_z_zzz_xzz[k] = g_zz_xzz[k] - g_0_z_zz_xzz[k] * ab_z + g_0_z_zz_xzzz[k];

                g_0_z_zzz_yyy[k] = g_zz_yyy[k] - g_0_z_zz_yyy[k] * ab_z + g_0_z_zz_yyyz[k];

                g_0_z_zzz_yyz[k] = g_zz_yyz[k] - g_0_z_zz_yyz[k] * ab_z + g_0_z_zz_yyzz[k];

                g_0_z_zzz_yzz[k] = g_zz_yzz[k] - g_0_z_zz_yzz[k] * ab_z + g_0_z_zz_yzzz[k];

                g_0_z_zzz_zzz[k] = g_zz_zzz[k] - g_0_z_zz_zzz[k] * ab_z + g_0_z_zz_zzzz[k];
            }
        }
    }
}

} // erirec namespace

