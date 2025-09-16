#include "ElectronRepulsionGeom1100ContrRecGPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom11_hrr_electron_repulsion_gpxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_gpxx,
                                            const size_t idx_geom_01_fpxx,
                                            const size_t idx_geom_10_fpxx,
                                            const size_t idx_geom_11_fpxx,
                                            const size_t idx_geom_11_fdxx,
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

            const auto fp_geom_11_off = idx_geom_11_fpxx + i * dcomps + j;

            auto g_x_x_xxx_x = cbuffer.data(fp_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xxx_y = cbuffer.data(fp_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xxx_z = cbuffer.data(fp_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_xxy_x = cbuffer.data(fp_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xxy_y = cbuffer.data(fp_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xxy_z = cbuffer.data(fp_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_xxz_x = cbuffer.data(fp_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xxz_y = cbuffer.data(fp_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xxz_z = cbuffer.data(fp_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_xyy_x = cbuffer.data(fp_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_xyy_y = cbuffer.data(fp_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_xyy_z = cbuffer.data(fp_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_xyz_x = cbuffer.data(fp_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_xyz_y = cbuffer.data(fp_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_xyz_z = cbuffer.data(fp_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_xzz_x = cbuffer.data(fp_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_xzz_y = cbuffer.data(fp_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_xzz_z = cbuffer.data(fp_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_yyy_x = cbuffer.data(fp_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_yyy_y = cbuffer.data(fp_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_yyy_z = cbuffer.data(fp_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_yyz_x = cbuffer.data(fp_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_yyz_y = cbuffer.data(fp_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_yyz_z = cbuffer.data(fp_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_yzz_x = cbuffer.data(fp_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_yzz_y = cbuffer.data(fp_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_yzz_z = cbuffer.data(fp_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_zzz_x = cbuffer.data(fp_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_zzz_y = cbuffer.data(fp_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_zzz_z = cbuffer.data(fp_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_y_xxx_x = cbuffer.data(fp_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_y_xxx_y = cbuffer.data(fp_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_y_xxx_z = cbuffer.data(fp_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_y_xxy_x = cbuffer.data(fp_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_y_xxy_y = cbuffer.data(fp_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_y_xxy_z = cbuffer.data(fp_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_y_xxz_x = cbuffer.data(fp_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_y_xxz_y = cbuffer.data(fp_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_y_xxz_z = cbuffer.data(fp_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_y_xyy_x = cbuffer.data(fp_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_y_xyy_y = cbuffer.data(fp_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_y_xyy_z = cbuffer.data(fp_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_y_xyz_x = cbuffer.data(fp_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_y_xyz_y = cbuffer.data(fp_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_y_xyz_z = cbuffer.data(fp_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_y_xzz_x = cbuffer.data(fp_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_y_xzz_y = cbuffer.data(fp_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_y_xzz_z = cbuffer.data(fp_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_y_yyy_x = cbuffer.data(fp_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_y_yyy_y = cbuffer.data(fp_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_y_yyy_z = cbuffer.data(fp_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_y_yyz_x = cbuffer.data(fp_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_y_yyz_y = cbuffer.data(fp_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_y_yyz_z = cbuffer.data(fp_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_y_yzz_x = cbuffer.data(fp_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_y_yzz_y = cbuffer.data(fp_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_y_yzz_z = cbuffer.data(fp_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_y_zzz_x = cbuffer.data(fp_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_y_zzz_y = cbuffer.data(fp_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_y_zzz_z = cbuffer.data(fp_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_z_xxx_x = cbuffer.data(fp_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_z_xxx_y = cbuffer.data(fp_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_z_xxx_z = cbuffer.data(fp_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_z_xxy_x = cbuffer.data(fp_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_z_xxy_y = cbuffer.data(fp_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_z_xxy_z = cbuffer.data(fp_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_z_xxz_x = cbuffer.data(fp_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_z_xxz_y = cbuffer.data(fp_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_z_xxz_z = cbuffer.data(fp_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_z_xyy_x = cbuffer.data(fp_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_z_xyy_y = cbuffer.data(fp_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_z_xyy_z = cbuffer.data(fp_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_z_xyz_x = cbuffer.data(fp_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_z_xyz_y = cbuffer.data(fp_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_z_xyz_z = cbuffer.data(fp_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_z_xzz_x = cbuffer.data(fp_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_z_xzz_y = cbuffer.data(fp_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_z_xzz_z = cbuffer.data(fp_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_z_yyy_x = cbuffer.data(fp_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_z_yyy_y = cbuffer.data(fp_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_z_yyy_z = cbuffer.data(fp_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_z_yyz_x = cbuffer.data(fp_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_z_yyz_y = cbuffer.data(fp_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_z_yyz_z = cbuffer.data(fp_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_z_yzz_x = cbuffer.data(fp_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_z_yzz_y = cbuffer.data(fp_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_z_yzz_z = cbuffer.data(fp_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_z_zzz_x = cbuffer.data(fp_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_z_zzz_y = cbuffer.data(fp_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_z_zzz_z = cbuffer.data(fp_geom_11_off + 89 * ccomps * dcomps);

            auto g_y_x_xxx_x = cbuffer.data(fp_geom_11_off + 90 * ccomps * dcomps);

            auto g_y_x_xxx_y = cbuffer.data(fp_geom_11_off + 91 * ccomps * dcomps);

            auto g_y_x_xxx_z = cbuffer.data(fp_geom_11_off + 92 * ccomps * dcomps);

            auto g_y_x_xxy_x = cbuffer.data(fp_geom_11_off + 93 * ccomps * dcomps);

            auto g_y_x_xxy_y = cbuffer.data(fp_geom_11_off + 94 * ccomps * dcomps);

            auto g_y_x_xxy_z = cbuffer.data(fp_geom_11_off + 95 * ccomps * dcomps);

            auto g_y_x_xxz_x = cbuffer.data(fp_geom_11_off + 96 * ccomps * dcomps);

            auto g_y_x_xxz_y = cbuffer.data(fp_geom_11_off + 97 * ccomps * dcomps);

            auto g_y_x_xxz_z = cbuffer.data(fp_geom_11_off + 98 * ccomps * dcomps);

            auto g_y_x_xyy_x = cbuffer.data(fp_geom_11_off + 99 * ccomps * dcomps);

            auto g_y_x_xyy_y = cbuffer.data(fp_geom_11_off + 100 * ccomps * dcomps);

            auto g_y_x_xyy_z = cbuffer.data(fp_geom_11_off + 101 * ccomps * dcomps);

            auto g_y_x_xyz_x = cbuffer.data(fp_geom_11_off + 102 * ccomps * dcomps);

            auto g_y_x_xyz_y = cbuffer.data(fp_geom_11_off + 103 * ccomps * dcomps);

            auto g_y_x_xyz_z = cbuffer.data(fp_geom_11_off + 104 * ccomps * dcomps);

            auto g_y_x_xzz_x = cbuffer.data(fp_geom_11_off + 105 * ccomps * dcomps);

            auto g_y_x_xzz_y = cbuffer.data(fp_geom_11_off + 106 * ccomps * dcomps);

            auto g_y_x_xzz_z = cbuffer.data(fp_geom_11_off + 107 * ccomps * dcomps);

            auto g_y_x_yyy_x = cbuffer.data(fp_geom_11_off + 108 * ccomps * dcomps);

            auto g_y_x_yyy_y = cbuffer.data(fp_geom_11_off + 109 * ccomps * dcomps);

            auto g_y_x_yyy_z = cbuffer.data(fp_geom_11_off + 110 * ccomps * dcomps);

            auto g_y_x_yyz_x = cbuffer.data(fp_geom_11_off + 111 * ccomps * dcomps);

            auto g_y_x_yyz_y = cbuffer.data(fp_geom_11_off + 112 * ccomps * dcomps);

            auto g_y_x_yyz_z = cbuffer.data(fp_geom_11_off + 113 * ccomps * dcomps);

            auto g_y_x_yzz_x = cbuffer.data(fp_geom_11_off + 114 * ccomps * dcomps);

            auto g_y_x_yzz_y = cbuffer.data(fp_geom_11_off + 115 * ccomps * dcomps);

            auto g_y_x_yzz_z = cbuffer.data(fp_geom_11_off + 116 * ccomps * dcomps);

            auto g_y_x_zzz_x = cbuffer.data(fp_geom_11_off + 117 * ccomps * dcomps);

            auto g_y_x_zzz_y = cbuffer.data(fp_geom_11_off + 118 * ccomps * dcomps);

            auto g_y_x_zzz_z = cbuffer.data(fp_geom_11_off + 119 * ccomps * dcomps);

            auto g_y_y_xxx_x = cbuffer.data(fp_geom_11_off + 120 * ccomps * dcomps);

            auto g_y_y_xxx_y = cbuffer.data(fp_geom_11_off + 121 * ccomps * dcomps);

            auto g_y_y_xxx_z = cbuffer.data(fp_geom_11_off + 122 * ccomps * dcomps);

            auto g_y_y_xxy_x = cbuffer.data(fp_geom_11_off + 123 * ccomps * dcomps);

            auto g_y_y_xxy_y = cbuffer.data(fp_geom_11_off + 124 * ccomps * dcomps);

            auto g_y_y_xxy_z = cbuffer.data(fp_geom_11_off + 125 * ccomps * dcomps);

            auto g_y_y_xxz_x = cbuffer.data(fp_geom_11_off + 126 * ccomps * dcomps);

            auto g_y_y_xxz_y = cbuffer.data(fp_geom_11_off + 127 * ccomps * dcomps);

            auto g_y_y_xxz_z = cbuffer.data(fp_geom_11_off + 128 * ccomps * dcomps);

            auto g_y_y_xyy_x = cbuffer.data(fp_geom_11_off + 129 * ccomps * dcomps);

            auto g_y_y_xyy_y = cbuffer.data(fp_geom_11_off + 130 * ccomps * dcomps);

            auto g_y_y_xyy_z = cbuffer.data(fp_geom_11_off + 131 * ccomps * dcomps);

            auto g_y_y_xyz_x = cbuffer.data(fp_geom_11_off + 132 * ccomps * dcomps);

            auto g_y_y_xyz_y = cbuffer.data(fp_geom_11_off + 133 * ccomps * dcomps);

            auto g_y_y_xyz_z = cbuffer.data(fp_geom_11_off + 134 * ccomps * dcomps);

            auto g_y_y_xzz_x = cbuffer.data(fp_geom_11_off + 135 * ccomps * dcomps);

            auto g_y_y_xzz_y = cbuffer.data(fp_geom_11_off + 136 * ccomps * dcomps);

            auto g_y_y_xzz_z = cbuffer.data(fp_geom_11_off + 137 * ccomps * dcomps);

            auto g_y_y_yyy_x = cbuffer.data(fp_geom_11_off + 138 * ccomps * dcomps);

            auto g_y_y_yyy_y = cbuffer.data(fp_geom_11_off + 139 * ccomps * dcomps);

            auto g_y_y_yyy_z = cbuffer.data(fp_geom_11_off + 140 * ccomps * dcomps);

            auto g_y_y_yyz_x = cbuffer.data(fp_geom_11_off + 141 * ccomps * dcomps);

            auto g_y_y_yyz_y = cbuffer.data(fp_geom_11_off + 142 * ccomps * dcomps);

            auto g_y_y_yyz_z = cbuffer.data(fp_geom_11_off + 143 * ccomps * dcomps);

            auto g_y_y_yzz_x = cbuffer.data(fp_geom_11_off + 144 * ccomps * dcomps);

            auto g_y_y_yzz_y = cbuffer.data(fp_geom_11_off + 145 * ccomps * dcomps);

            auto g_y_y_yzz_z = cbuffer.data(fp_geom_11_off + 146 * ccomps * dcomps);

            auto g_y_y_zzz_x = cbuffer.data(fp_geom_11_off + 147 * ccomps * dcomps);

            auto g_y_y_zzz_y = cbuffer.data(fp_geom_11_off + 148 * ccomps * dcomps);

            auto g_y_y_zzz_z = cbuffer.data(fp_geom_11_off + 149 * ccomps * dcomps);

            auto g_y_z_xxx_x = cbuffer.data(fp_geom_11_off + 150 * ccomps * dcomps);

            auto g_y_z_xxx_y = cbuffer.data(fp_geom_11_off + 151 * ccomps * dcomps);

            auto g_y_z_xxx_z = cbuffer.data(fp_geom_11_off + 152 * ccomps * dcomps);

            auto g_y_z_xxy_x = cbuffer.data(fp_geom_11_off + 153 * ccomps * dcomps);

            auto g_y_z_xxy_y = cbuffer.data(fp_geom_11_off + 154 * ccomps * dcomps);

            auto g_y_z_xxy_z = cbuffer.data(fp_geom_11_off + 155 * ccomps * dcomps);

            auto g_y_z_xxz_x = cbuffer.data(fp_geom_11_off + 156 * ccomps * dcomps);

            auto g_y_z_xxz_y = cbuffer.data(fp_geom_11_off + 157 * ccomps * dcomps);

            auto g_y_z_xxz_z = cbuffer.data(fp_geom_11_off + 158 * ccomps * dcomps);

            auto g_y_z_xyy_x = cbuffer.data(fp_geom_11_off + 159 * ccomps * dcomps);

            auto g_y_z_xyy_y = cbuffer.data(fp_geom_11_off + 160 * ccomps * dcomps);

            auto g_y_z_xyy_z = cbuffer.data(fp_geom_11_off + 161 * ccomps * dcomps);

            auto g_y_z_xyz_x = cbuffer.data(fp_geom_11_off + 162 * ccomps * dcomps);

            auto g_y_z_xyz_y = cbuffer.data(fp_geom_11_off + 163 * ccomps * dcomps);

            auto g_y_z_xyz_z = cbuffer.data(fp_geom_11_off + 164 * ccomps * dcomps);

            auto g_y_z_xzz_x = cbuffer.data(fp_geom_11_off + 165 * ccomps * dcomps);

            auto g_y_z_xzz_y = cbuffer.data(fp_geom_11_off + 166 * ccomps * dcomps);

            auto g_y_z_xzz_z = cbuffer.data(fp_geom_11_off + 167 * ccomps * dcomps);

            auto g_y_z_yyy_x = cbuffer.data(fp_geom_11_off + 168 * ccomps * dcomps);

            auto g_y_z_yyy_y = cbuffer.data(fp_geom_11_off + 169 * ccomps * dcomps);

            auto g_y_z_yyy_z = cbuffer.data(fp_geom_11_off + 170 * ccomps * dcomps);

            auto g_y_z_yyz_x = cbuffer.data(fp_geom_11_off + 171 * ccomps * dcomps);

            auto g_y_z_yyz_y = cbuffer.data(fp_geom_11_off + 172 * ccomps * dcomps);

            auto g_y_z_yyz_z = cbuffer.data(fp_geom_11_off + 173 * ccomps * dcomps);

            auto g_y_z_yzz_x = cbuffer.data(fp_geom_11_off + 174 * ccomps * dcomps);

            auto g_y_z_yzz_y = cbuffer.data(fp_geom_11_off + 175 * ccomps * dcomps);

            auto g_y_z_yzz_z = cbuffer.data(fp_geom_11_off + 176 * ccomps * dcomps);

            auto g_y_z_zzz_x = cbuffer.data(fp_geom_11_off + 177 * ccomps * dcomps);

            auto g_y_z_zzz_y = cbuffer.data(fp_geom_11_off + 178 * ccomps * dcomps);

            auto g_y_z_zzz_z = cbuffer.data(fp_geom_11_off + 179 * ccomps * dcomps);

            auto g_z_x_xxx_x = cbuffer.data(fp_geom_11_off + 180 * ccomps * dcomps);

            auto g_z_x_xxx_y = cbuffer.data(fp_geom_11_off + 181 * ccomps * dcomps);

            auto g_z_x_xxx_z = cbuffer.data(fp_geom_11_off + 182 * ccomps * dcomps);

            auto g_z_x_xxy_x = cbuffer.data(fp_geom_11_off + 183 * ccomps * dcomps);

            auto g_z_x_xxy_y = cbuffer.data(fp_geom_11_off + 184 * ccomps * dcomps);

            auto g_z_x_xxy_z = cbuffer.data(fp_geom_11_off + 185 * ccomps * dcomps);

            auto g_z_x_xxz_x = cbuffer.data(fp_geom_11_off + 186 * ccomps * dcomps);

            auto g_z_x_xxz_y = cbuffer.data(fp_geom_11_off + 187 * ccomps * dcomps);

            auto g_z_x_xxz_z = cbuffer.data(fp_geom_11_off + 188 * ccomps * dcomps);

            auto g_z_x_xyy_x = cbuffer.data(fp_geom_11_off + 189 * ccomps * dcomps);

            auto g_z_x_xyy_y = cbuffer.data(fp_geom_11_off + 190 * ccomps * dcomps);

            auto g_z_x_xyy_z = cbuffer.data(fp_geom_11_off + 191 * ccomps * dcomps);

            auto g_z_x_xyz_x = cbuffer.data(fp_geom_11_off + 192 * ccomps * dcomps);

            auto g_z_x_xyz_y = cbuffer.data(fp_geom_11_off + 193 * ccomps * dcomps);

            auto g_z_x_xyz_z = cbuffer.data(fp_geom_11_off + 194 * ccomps * dcomps);

            auto g_z_x_xzz_x = cbuffer.data(fp_geom_11_off + 195 * ccomps * dcomps);

            auto g_z_x_xzz_y = cbuffer.data(fp_geom_11_off + 196 * ccomps * dcomps);

            auto g_z_x_xzz_z = cbuffer.data(fp_geom_11_off + 197 * ccomps * dcomps);

            auto g_z_x_yyy_x = cbuffer.data(fp_geom_11_off + 198 * ccomps * dcomps);

            auto g_z_x_yyy_y = cbuffer.data(fp_geom_11_off + 199 * ccomps * dcomps);

            auto g_z_x_yyy_z = cbuffer.data(fp_geom_11_off + 200 * ccomps * dcomps);

            auto g_z_x_yyz_x = cbuffer.data(fp_geom_11_off + 201 * ccomps * dcomps);

            auto g_z_x_yyz_y = cbuffer.data(fp_geom_11_off + 202 * ccomps * dcomps);

            auto g_z_x_yyz_z = cbuffer.data(fp_geom_11_off + 203 * ccomps * dcomps);

            auto g_z_x_yzz_x = cbuffer.data(fp_geom_11_off + 204 * ccomps * dcomps);

            auto g_z_x_yzz_y = cbuffer.data(fp_geom_11_off + 205 * ccomps * dcomps);

            auto g_z_x_yzz_z = cbuffer.data(fp_geom_11_off + 206 * ccomps * dcomps);

            auto g_z_x_zzz_x = cbuffer.data(fp_geom_11_off + 207 * ccomps * dcomps);

            auto g_z_x_zzz_y = cbuffer.data(fp_geom_11_off + 208 * ccomps * dcomps);

            auto g_z_x_zzz_z = cbuffer.data(fp_geom_11_off + 209 * ccomps * dcomps);

            auto g_z_y_xxx_x = cbuffer.data(fp_geom_11_off + 210 * ccomps * dcomps);

            auto g_z_y_xxx_y = cbuffer.data(fp_geom_11_off + 211 * ccomps * dcomps);

            auto g_z_y_xxx_z = cbuffer.data(fp_geom_11_off + 212 * ccomps * dcomps);

            auto g_z_y_xxy_x = cbuffer.data(fp_geom_11_off + 213 * ccomps * dcomps);

            auto g_z_y_xxy_y = cbuffer.data(fp_geom_11_off + 214 * ccomps * dcomps);

            auto g_z_y_xxy_z = cbuffer.data(fp_geom_11_off + 215 * ccomps * dcomps);

            auto g_z_y_xxz_x = cbuffer.data(fp_geom_11_off + 216 * ccomps * dcomps);

            auto g_z_y_xxz_y = cbuffer.data(fp_geom_11_off + 217 * ccomps * dcomps);

            auto g_z_y_xxz_z = cbuffer.data(fp_geom_11_off + 218 * ccomps * dcomps);

            auto g_z_y_xyy_x = cbuffer.data(fp_geom_11_off + 219 * ccomps * dcomps);

            auto g_z_y_xyy_y = cbuffer.data(fp_geom_11_off + 220 * ccomps * dcomps);

            auto g_z_y_xyy_z = cbuffer.data(fp_geom_11_off + 221 * ccomps * dcomps);

            auto g_z_y_xyz_x = cbuffer.data(fp_geom_11_off + 222 * ccomps * dcomps);

            auto g_z_y_xyz_y = cbuffer.data(fp_geom_11_off + 223 * ccomps * dcomps);

            auto g_z_y_xyz_z = cbuffer.data(fp_geom_11_off + 224 * ccomps * dcomps);

            auto g_z_y_xzz_x = cbuffer.data(fp_geom_11_off + 225 * ccomps * dcomps);

            auto g_z_y_xzz_y = cbuffer.data(fp_geom_11_off + 226 * ccomps * dcomps);

            auto g_z_y_xzz_z = cbuffer.data(fp_geom_11_off + 227 * ccomps * dcomps);

            auto g_z_y_yyy_x = cbuffer.data(fp_geom_11_off + 228 * ccomps * dcomps);

            auto g_z_y_yyy_y = cbuffer.data(fp_geom_11_off + 229 * ccomps * dcomps);

            auto g_z_y_yyy_z = cbuffer.data(fp_geom_11_off + 230 * ccomps * dcomps);

            auto g_z_y_yyz_x = cbuffer.data(fp_geom_11_off + 231 * ccomps * dcomps);

            auto g_z_y_yyz_y = cbuffer.data(fp_geom_11_off + 232 * ccomps * dcomps);

            auto g_z_y_yyz_z = cbuffer.data(fp_geom_11_off + 233 * ccomps * dcomps);

            auto g_z_y_yzz_x = cbuffer.data(fp_geom_11_off + 234 * ccomps * dcomps);

            auto g_z_y_yzz_y = cbuffer.data(fp_geom_11_off + 235 * ccomps * dcomps);

            auto g_z_y_yzz_z = cbuffer.data(fp_geom_11_off + 236 * ccomps * dcomps);

            auto g_z_y_zzz_x = cbuffer.data(fp_geom_11_off + 237 * ccomps * dcomps);

            auto g_z_y_zzz_y = cbuffer.data(fp_geom_11_off + 238 * ccomps * dcomps);

            auto g_z_y_zzz_z = cbuffer.data(fp_geom_11_off + 239 * ccomps * dcomps);

            auto g_z_z_xxx_x = cbuffer.data(fp_geom_11_off + 240 * ccomps * dcomps);

            auto g_z_z_xxx_y = cbuffer.data(fp_geom_11_off + 241 * ccomps * dcomps);

            auto g_z_z_xxx_z = cbuffer.data(fp_geom_11_off + 242 * ccomps * dcomps);

            auto g_z_z_xxy_x = cbuffer.data(fp_geom_11_off + 243 * ccomps * dcomps);

            auto g_z_z_xxy_y = cbuffer.data(fp_geom_11_off + 244 * ccomps * dcomps);

            auto g_z_z_xxy_z = cbuffer.data(fp_geom_11_off + 245 * ccomps * dcomps);

            auto g_z_z_xxz_x = cbuffer.data(fp_geom_11_off + 246 * ccomps * dcomps);

            auto g_z_z_xxz_y = cbuffer.data(fp_geom_11_off + 247 * ccomps * dcomps);

            auto g_z_z_xxz_z = cbuffer.data(fp_geom_11_off + 248 * ccomps * dcomps);

            auto g_z_z_xyy_x = cbuffer.data(fp_geom_11_off + 249 * ccomps * dcomps);

            auto g_z_z_xyy_y = cbuffer.data(fp_geom_11_off + 250 * ccomps * dcomps);

            auto g_z_z_xyy_z = cbuffer.data(fp_geom_11_off + 251 * ccomps * dcomps);

            auto g_z_z_xyz_x = cbuffer.data(fp_geom_11_off + 252 * ccomps * dcomps);

            auto g_z_z_xyz_y = cbuffer.data(fp_geom_11_off + 253 * ccomps * dcomps);

            auto g_z_z_xyz_z = cbuffer.data(fp_geom_11_off + 254 * ccomps * dcomps);

            auto g_z_z_xzz_x = cbuffer.data(fp_geom_11_off + 255 * ccomps * dcomps);

            auto g_z_z_xzz_y = cbuffer.data(fp_geom_11_off + 256 * ccomps * dcomps);

            auto g_z_z_xzz_z = cbuffer.data(fp_geom_11_off + 257 * ccomps * dcomps);

            auto g_z_z_yyy_x = cbuffer.data(fp_geom_11_off + 258 * ccomps * dcomps);

            auto g_z_z_yyy_y = cbuffer.data(fp_geom_11_off + 259 * ccomps * dcomps);

            auto g_z_z_yyy_z = cbuffer.data(fp_geom_11_off + 260 * ccomps * dcomps);

            auto g_z_z_yyz_x = cbuffer.data(fp_geom_11_off + 261 * ccomps * dcomps);

            auto g_z_z_yyz_y = cbuffer.data(fp_geom_11_off + 262 * ccomps * dcomps);

            auto g_z_z_yyz_z = cbuffer.data(fp_geom_11_off + 263 * ccomps * dcomps);

            auto g_z_z_yzz_x = cbuffer.data(fp_geom_11_off + 264 * ccomps * dcomps);

            auto g_z_z_yzz_y = cbuffer.data(fp_geom_11_off + 265 * ccomps * dcomps);

            auto g_z_z_yzz_z = cbuffer.data(fp_geom_11_off + 266 * ccomps * dcomps);

            auto g_z_z_zzz_x = cbuffer.data(fp_geom_11_off + 267 * ccomps * dcomps);

            auto g_z_z_zzz_y = cbuffer.data(fp_geom_11_off + 268 * ccomps * dcomps);

            auto g_z_z_zzz_z = cbuffer.data(fp_geom_11_off + 269 * ccomps * dcomps);

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

            /// set up bra offset for contr_buffer_gpxx

            const auto gp_geom_11_off = idx_geom_11_gpxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxxx_x = cbuffer.data(gp_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xxxx_y = cbuffer.data(gp_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xxxx_z = cbuffer.data(gp_geom_11_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxx_x, g_0_x_xxx_y, g_0_x_xxx_z, g_x_0_xxx_x, g_x_0_xxx_y, g_x_0_xxx_z, g_x_x_xxx_x, g_x_x_xxx_xx, g_x_x_xxx_xy, g_x_x_xxx_xz, g_x_x_xxx_y, g_x_x_xxx_z, g_x_x_xxxx_x, g_x_x_xxxx_y, g_x_x_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxxx_x[k] = -g_0_x_xxx_x[k] + g_x_0_xxx_x[k] - g_x_x_xxx_x[k] * ab_x + g_x_x_xxx_xx[k];

                g_x_x_xxxx_y[k] = -g_0_x_xxx_y[k] + g_x_0_xxx_y[k] - g_x_x_xxx_y[k] * ab_x + g_x_x_xxx_xy[k];

                g_x_x_xxxx_z[k] = -g_0_x_xxx_z[k] + g_x_0_xxx_z[k] - g_x_x_xxx_z[k] * ab_x + g_x_x_xxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxxy_x = cbuffer.data(gp_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xxxy_y = cbuffer.data(gp_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xxxy_z = cbuffer.data(gp_geom_11_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxx_x, g_x_x_xxx_xy, g_x_x_xxx_y, g_x_x_xxx_yy, g_x_x_xxx_yz, g_x_x_xxx_z, g_x_x_xxxy_x, g_x_x_xxxy_y, g_x_x_xxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxxy_x[k] = -g_x_x_xxx_x[k] * ab_y + g_x_x_xxx_xy[k];

                g_x_x_xxxy_y[k] = -g_x_x_xxx_y[k] * ab_y + g_x_x_xxx_yy[k];

                g_x_x_xxxy_z[k] = -g_x_x_xxx_z[k] * ab_y + g_x_x_xxx_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxxz_x = cbuffer.data(gp_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xxxz_y = cbuffer.data(gp_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xxxz_z = cbuffer.data(gp_geom_11_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxx_x, g_x_x_xxx_xz, g_x_x_xxx_y, g_x_x_xxx_yz, g_x_x_xxx_z, g_x_x_xxx_zz, g_x_x_xxxz_x, g_x_x_xxxz_y, g_x_x_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxxz_x[k] = -g_x_x_xxx_x[k] * ab_z + g_x_x_xxx_xz[k];

                g_x_x_xxxz_y[k] = -g_x_x_xxx_y[k] * ab_z + g_x_x_xxx_yz[k];

                g_x_x_xxxz_z[k] = -g_x_x_xxx_z[k] * ab_z + g_x_x_xxx_zz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxyy_x = cbuffer.data(gp_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_xxyy_y = cbuffer.data(gp_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_xxyy_z = cbuffer.data(gp_geom_11_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxy_x, g_x_x_xxy_xy, g_x_x_xxy_y, g_x_x_xxy_yy, g_x_x_xxy_yz, g_x_x_xxy_z, g_x_x_xxyy_x, g_x_x_xxyy_y, g_x_x_xxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxyy_x[k] = -g_x_x_xxy_x[k] * ab_y + g_x_x_xxy_xy[k];

                g_x_x_xxyy_y[k] = -g_x_x_xxy_y[k] * ab_y + g_x_x_xxy_yy[k];

                g_x_x_xxyy_z[k] = -g_x_x_xxy_z[k] * ab_y + g_x_x_xxy_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxyz_x = cbuffer.data(gp_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_xxyz_y = cbuffer.data(gp_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_xxyz_z = cbuffer.data(gp_geom_11_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxyz_x, g_x_x_xxyz_y, g_x_x_xxyz_z, g_x_x_xxz_x, g_x_x_xxz_xy, g_x_x_xxz_y, g_x_x_xxz_yy, g_x_x_xxz_yz, g_x_x_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxyz_x[k] = -g_x_x_xxz_x[k] * ab_y + g_x_x_xxz_xy[k];

                g_x_x_xxyz_y[k] = -g_x_x_xxz_y[k] * ab_y + g_x_x_xxz_yy[k];

                g_x_x_xxyz_z[k] = -g_x_x_xxz_z[k] * ab_y + g_x_x_xxz_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_x_x_xxzz_x = cbuffer.data(gp_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_xxzz_y = cbuffer.data(gp_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_xxzz_z = cbuffer.data(gp_geom_11_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xxz_x, g_x_x_xxz_xz, g_x_x_xxz_y, g_x_x_xxz_yz, g_x_x_xxz_z, g_x_x_xxz_zz, g_x_x_xxzz_x, g_x_x_xxzz_y, g_x_x_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xxzz_x[k] = -g_x_x_xxz_x[k] * ab_z + g_x_x_xxz_xz[k];

                g_x_x_xxzz_y[k] = -g_x_x_xxz_y[k] * ab_z + g_x_x_xxz_yz[k];

                g_x_x_xxzz_z[k] = -g_x_x_xxz_z[k] * ab_z + g_x_x_xxz_zz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyyy_x = cbuffer.data(gp_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_xyyy_y = cbuffer.data(gp_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_xyyy_z = cbuffer.data(gp_geom_11_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xyy_x, g_x_x_xyy_xy, g_x_x_xyy_y, g_x_x_xyy_yy, g_x_x_xyy_yz, g_x_x_xyy_z, g_x_x_xyyy_x, g_x_x_xyyy_y, g_x_x_xyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyyy_x[k] = -g_x_x_xyy_x[k] * ab_y + g_x_x_xyy_xy[k];

                g_x_x_xyyy_y[k] = -g_x_x_xyy_y[k] * ab_y + g_x_x_xyy_yy[k];

                g_x_x_xyyy_z[k] = -g_x_x_xyy_z[k] * ab_y + g_x_x_xyy_yz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyyz_x = cbuffer.data(gp_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_xyyz_y = cbuffer.data(gp_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_xyyz_z = cbuffer.data(gp_geom_11_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xyyz_x, g_x_x_xyyz_y, g_x_x_xyyz_z, g_x_x_xyz_x, g_x_x_xyz_xy, g_x_x_xyz_y, g_x_x_xyz_yy, g_x_x_xyz_yz, g_x_x_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyyz_x[k] = -g_x_x_xyz_x[k] * ab_y + g_x_x_xyz_xy[k];

                g_x_x_xyyz_y[k] = -g_x_x_xyz_y[k] * ab_y + g_x_x_xyz_yy[k];

                g_x_x_xyyz_z[k] = -g_x_x_xyz_z[k] * ab_y + g_x_x_xyz_yz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_x_x_xyzz_x = cbuffer.data(gp_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_xyzz_y = cbuffer.data(gp_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_xyzz_z = cbuffer.data(gp_geom_11_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xyzz_x, g_x_x_xyzz_y, g_x_x_xyzz_z, g_x_x_xzz_x, g_x_x_xzz_xy, g_x_x_xzz_y, g_x_x_xzz_yy, g_x_x_xzz_yz, g_x_x_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xyzz_x[k] = -g_x_x_xzz_x[k] * ab_y + g_x_x_xzz_xy[k];

                g_x_x_xyzz_y[k] = -g_x_x_xzz_y[k] * ab_y + g_x_x_xzz_yy[k];

                g_x_x_xyzz_z[k] = -g_x_x_xzz_z[k] * ab_y + g_x_x_xzz_yz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_x_x_xzzz_x = cbuffer.data(gp_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_xzzz_y = cbuffer.data(gp_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_xzzz_z = cbuffer.data(gp_geom_11_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_xzz_x, g_x_x_xzz_xz, g_x_x_xzz_y, g_x_x_xzz_yz, g_x_x_xzz_z, g_x_x_xzz_zz, g_x_x_xzzz_x, g_x_x_xzzz_y, g_x_x_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xzzz_x[k] = -g_x_x_xzz_x[k] * ab_z + g_x_x_xzz_xz[k];

                g_x_x_xzzz_y[k] = -g_x_x_xzz_y[k] * ab_z + g_x_x_xzz_yz[k];

                g_x_x_xzzz_z[k] = -g_x_x_xzz_z[k] * ab_z + g_x_x_xzz_zz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyyy_x = cbuffer.data(gp_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_yyyy_y = cbuffer.data(gp_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_yyyy_z = cbuffer.data(gp_geom_11_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yyy_x, g_x_x_yyy_xy, g_x_x_yyy_y, g_x_x_yyy_yy, g_x_x_yyy_yz, g_x_x_yyy_z, g_x_x_yyyy_x, g_x_x_yyyy_y, g_x_x_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyyy_x[k] = -g_x_x_yyy_x[k] * ab_y + g_x_x_yyy_xy[k];

                g_x_x_yyyy_y[k] = -g_x_x_yyy_y[k] * ab_y + g_x_x_yyy_yy[k];

                g_x_x_yyyy_z[k] = -g_x_x_yyy_z[k] * ab_y + g_x_x_yyy_yz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyyz_x = cbuffer.data(gp_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_yyyz_y = cbuffer.data(gp_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_yyyz_z = cbuffer.data(gp_geom_11_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yyyz_x, g_x_x_yyyz_y, g_x_x_yyyz_z, g_x_x_yyz_x, g_x_x_yyz_xy, g_x_x_yyz_y, g_x_x_yyz_yy, g_x_x_yyz_yz, g_x_x_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyyz_x[k] = -g_x_x_yyz_x[k] * ab_y + g_x_x_yyz_xy[k];

                g_x_x_yyyz_y[k] = -g_x_x_yyz_y[k] * ab_y + g_x_x_yyz_yy[k];

                g_x_x_yyyz_z[k] = -g_x_x_yyz_z[k] * ab_y + g_x_x_yyz_yz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_x_x_yyzz_x = cbuffer.data(gp_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_yyzz_y = cbuffer.data(gp_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_yyzz_z = cbuffer.data(gp_geom_11_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yyzz_x, g_x_x_yyzz_y, g_x_x_yyzz_z, g_x_x_yzz_x, g_x_x_yzz_xy, g_x_x_yzz_y, g_x_x_yzz_yy, g_x_x_yzz_yz, g_x_x_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yyzz_x[k] = -g_x_x_yzz_x[k] * ab_y + g_x_x_yzz_xy[k];

                g_x_x_yyzz_y[k] = -g_x_x_yzz_y[k] * ab_y + g_x_x_yzz_yy[k];

                g_x_x_yyzz_z[k] = -g_x_x_yzz_z[k] * ab_y + g_x_x_yzz_yz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_x_x_yzzz_x = cbuffer.data(gp_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_yzzz_y = cbuffer.data(gp_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_yzzz_z = cbuffer.data(gp_geom_11_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yzzz_x, g_x_x_yzzz_y, g_x_x_yzzz_z, g_x_x_zzz_x, g_x_x_zzz_xy, g_x_x_zzz_y, g_x_x_zzz_yy, g_x_x_zzz_yz, g_x_x_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yzzz_x[k] = -g_x_x_zzz_x[k] * ab_y + g_x_x_zzz_xy[k];

                g_x_x_yzzz_y[k] = -g_x_x_zzz_y[k] * ab_y + g_x_x_zzz_yy[k];

                g_x_x_yzzz_z[k] = -g_x_x_zzz_z[k] * ab_y + g_x_x_zzz_yz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_x_x_zzzz_x = cbuffer.data(gp_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_zzzz_y = cbuffer.data(gp_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_zzzz_z = cbuffer.data(gp_geom_11_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_zzz_x, g_x_x_zzz_xz, g_x_x_zzz_y, g_x_x_zzz_yz, g_x_x_zzz_z, g_x_x_zzz_zz, g_x_x_zzzz_x, g_x_x_zzzz_y, g_x_x_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_zzzz_x[k] = -g_x_x_zzz_x[k] * ab_z + g_x_x_zzz_xz[k];

                g_x_x_zzzz_y[k] = -g_x_x_zzz_y[k] * ab_z + g_x_x_zzz_yz[k];

                g_x_x_zzzz_z[k] = -g_x_x_zzz_z[k] * ab_z + g_x_x_zzz_zz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxxx_x = cbuffer.data(gp_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_y_xxxx_y = cbuffer.data(gp_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_y_xxxx_z = cbuffer.data(gp_geom_11_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxx_x, g_0_y_xxx_y, g_0_y_xxx_z, g_x_y_xxx_x, g_x_y_xxx_xx, g_x_y_xxx_xy, g_x_y_xxx_xz, g_x_y_xxx_y, g_x_y_xxx_z, g_x_y_xxxx_x, g_x_y_xxxx_y, g_x_y_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxxx_x[k] = -g_0_y_xxx_x[k] - g_x_y_xxx_x[k] * ab_x + g_x_y_xxx_xx[k];

                g_x_y_xxxx_y[k] = -g_0_y_xxx_y[k] - g_x_y_xxx_y[k] * ab_x + g_x_y_xxx_xy[k];

                g_x_y_xxxx_z[k] = -g_0_y_xxx_z[k] - g_x_y_xxx_z[k] * ab_x + g_x_y_xxx_xz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxxy_x = cbuffer.data(gp_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_y_xxxy_y = cbuffer.data(gp_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_y_xxxy_z = cbuffer.data(gp_geom_11_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxy_x, g_0_y_xxy_y, g_0_y_xxy_z, g_x_y_xxxy_x, g_x_y_xxxy_y, g_x_y_xxxy_z, g_x_y_xxy_x, g_x_y_xxy_xx, g_x_y_xxy_xy, g_x_y_xxy_xz, g_x_y_xxy_y, g_x_y_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxxy_x[k] = -g_0_y_xxy_x[k] - g_x_y_xxy_x[k] * ab_x + g_x_y_xxy_xx[k];

                g_x_y_xxxy_y[k] = -g_0_y_xxy_y[k] - g_x_y_xxy_y[k] * ab_x + g_x_y_xxy_xy[k];

                g_x_y_xxxy_z[k] = -g_0_y_xxy_z[k] - g_x_y_xxy_z[k] * ab_x + g_x_y_xxy_xz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxxz_x = cbuffer.data(gp_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_y_xxxz_y = cbuffer.data(gp_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_y_xxxz_z = cbuffer.data(gp_geom_11_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xxx_x, g_x_y_xxx_xz, g_x_y_xxx_y, g_x_y_xxx_yz, g_x_y_xxx_z, g_x_y_xxx_zz, g_x_y_xxxz_x, g_x_y_xxxz_y, g_x_y_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxxz_x[k] = -g_x_y_xxx_x[k] * ab_z + g_x_y_xxx_xz[k];

                g_x_y_xxxz_y[k] = -g_x_y_xxx_y[k] * ab_z + g_x_y_xxx_yz[k];

                g_x_y_xxxz_z[k] = -g_x_y_xxx_z[k] * ab_z + g_x_y_xxx_zz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxyy_x = cbuffer.data(gp_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_y_xxyy_y = cbuffer.data(gp_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_y_xxyy_z = cbuffer.data(gp_geom_11_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyy_x, g_0_y_xyy_y, g_0_y_xyy_z, g_x_y_xxyy_x, g_x_y_xxyy_y, g_x_y_xxyy_z, g_x_y_xyy_x, g_x_y_xyy_xx, g_x_y_xyy_xy, g_x_y_xyy_xz, g_x_y_xyy_y, g_x_y_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxyy_x[k] = -g_0_y_xyy_x[k] - g_x_y_xyy_x[k] * ab_x + g_x_y_xyy_xx[k];

                g_x_y_xxyy_y[k] = -g_0_y_xyy_y[k] - g_x_y_xyy_y[k] * ab_x + g_x_y_xyy_xy[k];

                g_x_y_xxyy_z[k] = -g_0_y_xyy_z[k] - g_x_y_xyy_z[k] * ab_x + g_x_y_xyy_xz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxyz_x = cbuffer.data(gp_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_y_xxyz_y = cbuffer.data(gp_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_y_xxyz_z = cbuffer.data(gp_geom_11_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xxy_x, g_x_y_xxy_xz, g_x_y_xxy_y, g_x_y_xxy_yz, g_x_y_xxy_z, g_x_y_xxy_zz, g_x_y_xxyz_x, g_x_y_xxyz_y, g_x_y_xxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxyz_x[k] = -g_x_y_xxy_x[k] * ab_z + g_x_y_xxy_xz[k];

                g_x_y_xxyz_y[k] = -g_x_y_xxy_y[k] * ab_z + g_x_y_xxy_yz[k];

                g_x_y_xxyz_z[k] = -g_x_y_xxy_z[k] * ab_z + g_x_y_xxy_zz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_x_y_xxzz_x = cbuffer.data(gp_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_y_xxzz_y = cbuffer.data(gp_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_y_xxzz_z = cbuffer.data(gp_geom_11_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xxz_x, g_x_y_xxz_xz, g_x_y_xxz_y, g_x_y_xxz_yz, g_x_y_xxz_z, g_x_y_xxz_zz, g_x_y_xxzz_x, g_x_y_xxzz_y, g_x_y_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xxzz_x[k] = -g_x_y_xxz_x[k] * ab_z + g_x_y_xxz_xz[k];

                g_x_y_xxzz_y[k] = -g_x_y_xxz_y[k] * ab_z + g_x_y_xxz_yz[k];

                g_x_y_xxzz_z[k] = -g_x_y_xxz_z[k] * ab_z + g_x_y_xxz_zz[k];
            }

            /// Set up 63-66 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyyy_x = cbuffer.data(gp_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_y_xyyy_y = cbuffer.data(gp_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_y_xyyy_z = cbuffer.data(gp_geom_11_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyy_x, g_0_y_yyy_y, g_0_y_yyy_z, g_x_y_xyyy_x, g_x_y_xyyy_y, g_x_y_xyyy_z, g_x_y_yyy_x, g_x_y_yyy_xx, g_x_y_yyy_xy, g_x_y_yyy_xz, g_x_y_yyy_y, g_x_y_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyyy_x[k] = -g_0_y_yyy_x[k] - g_x_y_yyy_x[k] * ab_x + g_x_y_yyy_xx[k];

                g_x_y_xyyy_y[k] = -g_0_y_yyy_y[k] - g_x_y_yyy_y[k] * ab_x + g_x_y_yyy_xy[k];

                g_x_y_xyyy_z[k] = -g_0_y_yyy_z[k] - g_x_y_yyy_z[k] * ab_x + g_x_y_yyy_xz[k];
            }

            /// Set up 66-69 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyyz_x = cbuffer.data(gp_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_y_xyyz_y = cbuffer.data(gp_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_y_xyyz_z = cbuffer.data(gp_geom_11_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xyy_x, g_x_y_xyy_xz, g_x_y_xyy_y, g_x_y_xyy_yz, g_x_y_xyy_z, g_x_y_xyy_zz, g_x_y_xyyz_x, g_x_y_xyyz_y, g_x_y_xyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyyz_x[k] = -g_x_y_xyy_x[k] * ab_z + g_x_y_xyy_xz[k];

                g_x_y_xyyz_y[k] = -g_x_y_xyy_y[k] * ab_z + g_x_y_xyy_yz[k];

                g_x_y_xyyz_z[k] = -g_x_y_xyy_z[k] * ab_z + g_x_y_xyy_zz[k];
            }

            /// Set up 69-72 components of targeted buffer : cbuffer.data(

            auto g_x_y_xyzz_x = cbuffer.data(gp_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_y_xyzz_y = cbuffer.data(gp_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_y_xyzz_z = cbuffer.data(gp_geom_11_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xyz_x, g_x_y_xyz_xz, g_x_y_xyz_y, g_x_y_xyz_yz, g_x_y_xyz_z, g_x_y_xyz_zz, g_x_y_xyzz_x, g_x_y_xyzz_y, g_x_y_xyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xyzz_x[k] = -g_x_y_xyz_x[k] * ab_z + g_x_y_xyz_xz[k];

                g_x_y_xyzz_y[k] = -g_x_y_xyz_y[k] * ab_z + g_x_y_xyz_yz[k];

                g_x_y_xyzz_z[k] = -g_x_y_xyz_z[k] * ab_z + g_x_y_xyz_zz[k];
            }

            /// Set up 72-75 components of targeted buffer : cbuffer.data(

            auto g_x_y_xzzz_x = cbuffer.data(gp_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_y_xzzz_y = cbuffer.data(gp_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_y_xzzz_z = cbuffer.data(gp_geom_11_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_xzz_x, g_x_y_xzz_xz, g_x_y_xzz_y, g_x_y_xzz_yz, g_x_y_xzz_z, g_x_y_xzz_zz, g_x_y_xzzz_x, g_x_y_xzzz_y, g_x_y_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xzzz_x[k] = -g_x_y_xzz_x[k] * ab_z + g_x_y_xzz_xz[k];

                g_x_y_xzzz_y[k] = -g_x_y_xzz_y[k] * ab_z + g_x_y_xzz_yz[k];

                g_x_y_xzzz_z[k] = -g_x_y_xzz_z[k] * ab_z + g_x_y_xzz_zz[k];
            }

            /// Set up 75-78 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyyy_x = cbuffer.data(gp_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_y_yyyy_y = cbuffer.data(gp_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_y_yyyy_z = cbuffer.data(gp_geom_11_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyy_x, g_x_0_yyy_y, g_x_0_yyy_z, g_x_y_yyy_x, g_x_y_yyy_xy, g_x_y_yyy_y, g_x_y_yyy_yy, g_x_y_yyy_yz, g_x_y_yyy_z, g_x_y_yyyy_x, g_x_y_yyyy_y, g_x_y_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyyy_x[k] = g_x_0_yyy_x[k] - g_x_y_yyy_x[k] * ab_y + g_x_y_yyy_xy[k];

                g_x_y_yyyy_y[k] = g_x_0_yyy_y[k] - g_x_y_yyy_y[k] * ab_y + g_x_y_yyy_yy[k];

                g_x_y_yyyy_z[k] = g_x_0_yyy_z[k] - g_x_y_yyy_z[k] * ab_y + g_x_y_yyy_yz[k];
            }

            /// Set up 78-81 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyyz_x = cbuffer.data(gp_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_y_yyyz_y = cbuffer.data(gp_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_y_yyyz_z = cbuffer.data(gp_geom_11_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yyy_x, g_x_y_yyy_xz, g_x_y_yyy_y, g_x_y_yyy_yz, g_x_y_yyy_z, g_x_y_yyy_zz, g_x_y_yyyz_x, g_x_y_yyyz_y, g_x_y_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyyz_x[k] = -g_x_y_yyy_x[k] * ab_z + g_x_y_yyy_xz[k];

                g_x_y_yyyz_y[k] = -g_x_y_yyy_y[k] * ab_z + g_x_y_yyy_yz[k];

                g_x_y_yyyz_z[k] = -g_x_y_yyy_z[k] * ab_z + g_x_y_yyy_zz[k];
            }

            /// Set up 81-84 components of targeted buffer : cbuffer.data(

            auto g_x_y_yyzz_x = cbuffer.data(gp_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_y_yyzz_y = cbuffer.data(gp_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_y_yyzz_z = cbuffer.data(gp_geom_11_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yyz_x, g_x_y_yyz_xz, g_x_y_yyz_y, g_x_y_yyz_yz, g_x_y_yyz_z, g_x_y_yyz_zz, g_x_y_yyzz_x, g_x_y_yyzz_y, g_x_y_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yyzz_x[k] = -g_x_y_yyz_x[k] * ab_z + g_x_y_yyz_xz[k];

                g_x_y_yyzz_y[k] = -g_x_y_yyz_y[k] * ab_z + g_x_y_yyz_yz[k];

                g_x_y_yyzz_z[k] = -g_x_y_yyz_z[k] * ab_z + g_x_y_yyz_zz[k];
            }

            /// Set up 84-87 components of targeted buffer : cbuffer.data(

            auto g_x_y_yzzz_x = cbuffer.data(gp_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_y_yzzz_y = cbuffer.data(gp_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_y_yzzz_z = cbuffer.data(gp_geom_11_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_yzz_x, g_x_y_yzz_xz, g_x_y_yzz_y, g_x_y_yzz_yz, g_x_y_yzz_z, g_x_y_yzz_zz, g_x_y_yzzz_x, g_x_y_yzzz_y, g_x_y_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yzzz_x[k] = -g_x_y_yzz_x[k] * ab_z + g_x_y_yzz_xz[k];

                g_x_y_yzzz_y[k] = -g_x_y_yzz_y[k] * ab_z + g_x_y_yzz_yz[k];

                g_x_y_yzzz_z[k] = -g_x_y_yzz_z[k] * ab_z + g_x_y_yzz_zz[k];
            }

            /// Set up 87-90 components of targeted buffer : cbuffer.data(

            auto g_x_y_zzzz_x = cbuffer.data(gp_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_y_zzzz_y = cbuffer.data(gp_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_y_zzzz_z = cbuffer.data(gp_geom_11_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_zzz_x, g_x_y_zzz_xz, g_x_y_zzz_y, g_x_y_zzz_yz, g_x_y_zzz_z, g_x_y_zzz_zz, g_x_y_zzzz_x, g_x_y_zzzz_y, g_x_y_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_zzzz_x[k] = -g_x_y_zzz_x[k] * ab_z + g_x_y_zzz_xz[k];

                g_x_y_zzzz_y[k] = -g_x_y_zzz_y[k] * ab_z + g_x_y_zzz_yz[k];

                g_x_y_zzzz_z[k] = -g_x_y_zzz_z[k] * ab_z + g_x_y_zzz_zz[k];
            }

            /// Set up 90-93 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxxx_x = cbuffer.data(gp_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_z_xxxx_y = cbuffer.data(gp_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_z_xxxx_z = cbuffer.data(gp_geom_11_off + 92 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxx_x, g_0_z_xxx_y, g_0_z_xxx_z, g_x_z_xxx_x, g_x_z_xxx_xx, g_x_z_xxx_xy, g_x_z_xxx_xz, g_x_z_xxx_y, g_x_z_xxx_z, g_x_z_xxxx_x, g_x_z_xxxx_y, g_x_z_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxxx_x[k] = -g_0_z_xxx_x[k] - g_x_z_xxx_x[k] * ab_x + g_x_z_xxx_xx[k];

                g_x_z_xxxx_y[k] = -g_0_z_xxx_y[k] - g_x_z_xxx_y[k] * ab_x + g_x_z_xxx_xy[k];

                g_x_z_xxxx_z[k] = -g_0_z_xxx_z[k] - g_x_z_xxx_z[k] * ab_x + g_x_z_xxx_xz[k];
            }

            /// Set up 93-96 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxxy_x = cbuffer.data(gp_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_z_xxxy_y = cbuffer.data(gp_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_z_xxxy_z = cbuffer.data(gp_geom_11_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xxx_x, g_x_z_xxx_xy, g_x_z_xxx_y, g_x_z_xxx_yy, g_x_z_xxx_yz, g_x_z_xxx_z, g_x_z_xxxy_x, g_x_z_xxxy_y, g_x_z_xxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxxy_x[k] = -g_x_z_xxx_x[k] * ab_y + g_x_z_xxx_xy[k];

                g_x_z_xxxy_y[k] = -g_x_z_xxx_y[k] * ab_y + g_x_z_xxx_yy[k];

                g_x_z_xxxy_z[k] = -g_x_z_xxx_z[k] * ab_y + g_x_z_xxx_yz[k];
            }

            /// Set up 96-99 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxxz_x = cbuffer.data(gp_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_z_xxxz_y = cbuffer.data(gp_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_z_xxxz_z = cbuffer.data(gp_geom_11_off + 98 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxz_x, g_0_z_xxz_y, g_0_z_xxz_z, g_x_z_xxxz_x, g_x_z_xxxz_y, g_x_z_xxxz_z, g_x_z_xxz_x, g_x_z_xxz_xx, g_x_z_xxz_xy, g_x_z_xxz_xz, g_x_z_xxz_y, g_x_z_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxxz_x[k] = -g_0_z_xxz_x[k] - g_x_z_xxz_x[k] * ab_x + g_x_z_xxz_xx[k];

                g_x_z_xxxz_y[k] = -g_0_z_xxz_y[k] - g_x_z_xxz_y[k] * ab_x + g_x_z_xxz_xy[k];

                g_x_z_xxxz_z[k] = -g_0_z_xxz_z[k] - g_x_z_xxz_z[k] * ab_x + g_x_z_xxz_xz[k];
            }

            /// Set up 99-102 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxyy_x = cbuffer.data(gp_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_z_xxyy_y = cbuffer.data(gp_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_z_xxyy_z = cbuffer.data(gp_geom_11_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xxy_x, g_x_z_xxy_xy, g_x_z_xxy_y, g_x_z_xxy_yy, g_x_z_xxy_yz, g_x_z_xxy_z, g_x_z_xxyy_x, g_x_z_xxyy_y, g_x_z_xxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxyy_x[k] = -g_x_z_xxy_x[k] * ab_y + g_x_z_xxy_xy[k];

                g_x_z_xxyy_y[k] = -g_x_z_xxy_y[k] * ab_y + g_x_z_xxy_yy[k];

                g_x_z_xxyy_z[k] = -g_x_z_xxy_z[k] * ab_y + g_x_z_xxy_yz[k];
            }

            /// Set up 102-105 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxyz_x = cbuffer.data(gp_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_z_xxyz_y = cbuffer.data(gp_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_z_xxyz_z = cbuffer.data(gp_geom_11_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xxyz_x, g_x_z_xxyz_y, g_x_z_xxyz_z, g_x_z_xxz_x, g_x_z_xxz_xy, g_x_z_xxz_y, g_x_z_xxz_yy, g_x_z_xxz_yz, g_x_z_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxyz_x[k] = -g_x_z_xxz_x[k] * ab_y + g_x_z_xxz_xy[k];

                g_x_z_xxyz_y[k] = -g_x_z_xxz_y[k] * ab_y + g_x_z_xxz_yy[k];

                g_x_z_xxyz_z[k] = -g_x_z_xxz_z[k] * ab_y + g_x_z_xxz_yz[k];
            }

            /// Set up 105-108 components of targeted buffer : cbuffer.data(

            auto g_x_z_xxzz_x = cbuffer.data(gp_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_z_xxzz_y = cbuffer.data(gp_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_z_xxzz_z = cbuffer.data(gp_geom_11_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzz_x, g_0_z_xzz_y, g_0_z_xzz_z, g_x_z_xxzz_x, g_x_z_xxzz_y, g_x_z_xxzz_z, g_x_z_xzz_x, g_x_z_xzz_xx, g_x_z_xzz_xy, g_x_z_xzz_xz, g_x_z_xzz_y, g_x_z_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xxzz_x[k] = -g_0_z_xzz_x[k] - g_x_z_xzz_x[k] * ab_x + g_x_z_xzz_xx[k];

                g_x_z_xxzz_y[k] = -g_0_z_xzz_y[k] - g_x_z_xzz_y[k] * ab_x + g_x_z_xzz_xy[k];

                g_x_z_xxzz_z[k] = -g_0_z_xzz_z[k] - g_x_z_xzz_z[k] * ab_x + g_x_z_xzz_xz[k];
            }

            /// Set up 108-111 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyyy_x = cbuffer.data(gp_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_z_xyyy_y = cbuffer.data(gp_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_z_xyyy_z = cbuffer.data(gp_geom_11_off + 110 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xyy_x, g_x_z_xyy_xy, g_x_z_xyy_y, g_x_z_xyy_yy, g_x_z_xyy_yz, g_x_z_xyy_z, g_x_z_xyyy_x, g_x_z_xyyy_y, g_x_z_xyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyyy_x[k] = -g_x_z_xyy_x[k] * ab_y + g_x_z_xyy_xy[k];

                g_x_z_xyyy_y[k] = -g_x_z_xyy_y[k] * ab_y + g_x_z_xyy_yy[k];

                g_x_z_xyyy_z[k] = -g_x_z_xyy_z[k] * ab_y + g_x_z_xyy_yz[k];
            }

            /// Set up 111-114 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyyz_x = cbuffer.data(gp_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_z_xyyz_y = cbuffer.data(gp_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_z_xyyz_z = cbuffer.data(gp_geom_11_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xyyz_x, g_x_z_xyyz_y, g_x_z_xyyz_z, g_x_z_xyz_x, g_x_z_xyz_xy, g_x_z_xyz_y, g_x_z_xyz_yy, g_x_z_xyz_yz, g_x_z_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyyz_x[k] = -g_x_z_xyz_x[k] * ab_y + g_x_z_xyz_xy[k];

                g_x_z_xyyz_y[k] = -g_x_z_xyz_y[k] * ab_y + g_x_z_xyz_yy[k];

                g_x_z_xyyz_z[k] = -g_x_z_xyz_z[k] * ab_y + g_x_z_xyz_yz[k];
            }

            /// Set up 114-117 components of targeted buffer : cbuffer.data(

            auto g_x_z_xyzz_x = cbuffer.data(gp_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_z_xyzz_y = cbuffer.data(gp_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_z_xyzz_z = cbuffer.data(gp_geom_11_off + 116 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_xyzz_x, g_x_z_xyzz_y, g_x_z_xyzz_z, g_x_z_xzz_x, g_x_z_xzz_xy, g_x_z_xzz_y, g_x_z_xzz_yy, g_x_z_xzz_yz, g_x_z_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xyzz_x[k] = -g_x_z_xzz_x[k] * ab_y + g_x_z_xzz_xy[k];

                g_x_z_xyzz_y[k] = -g_x_z_xzz_y[k] * ab_y + g_x_z_xzz_yy[k];

                g_x_z_xyzz_z[k] = -g_x_z_xzz_z[k] * ab_y + g_x_z_xzz_yz[k];
            }

            /// Set up 117-120 components of targeted buffer : cbuffer.data(

            auto g_x_z_xzzz_x = cbuffer.data(gp_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_z_xzzz_y = cbuffer.data(gp_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_z_xzzz_z = cbuffer.data(gp_geom_11_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzz_x, g_0_z_zzz_y, g_0_z_zzz_z, g_x_z_xzzz_x, g_x_z_xzzz_y, g_x_z_xzzz_z, g_x_z_zzz_x, g_x_z_zzz_xx, g_x_z_zzz_xy, g_x_z_zzz_xz, g_x_z_zzz_y, g_x_z_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xzzz_x[k] = -g_0_z_zzz_x[k] - g_x_z_zzz_x[k] * ab_x + g_x_z_zzz_xx[k];

                g_x_z_xzzz_y[k] = -g_0_z_zzz_y[k] - g_x_z_zzz_y[k] * ab_x + g_x_z_zzz_xy[k];

                g_x_z_xzzz_z[k] = -g_0_z_zzz_z[k] - g_x_z_zzz_z[k] * ab_x + g_x_z_zzz_xz[k];
            }

            /// Set up 120-123 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyyy_x = cbuffer.data(gp_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_z_yyyy_y = cbuffer.data(gp_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_z_yyyy_z = cbuffer.data(gp_geom_11_off + 122 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yyy_x, g_x_z_yyy_xy, g_x_z_yyy_y, g_x_z_yyy_yy, g_x_z_yyy_yz, g_x_z_yyy_z, g_x_z_yyyy_x, g_x_z_yyyy_y, g_x_z_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyyy_x[k] = -g_x_z_yyy_x[k] * ab_y + g_x_z_yyy_xy[k];

                g_x_z_yyyy_y[k] = -g_x_z_yyy_y[k] * ab_y + g_x_z_yyy_yy[k];

                g_x_z_yyyy_z[k] = -g_x_z_yyy_z[k] * ab_y + g_x_z_yyy_yz[k];
            }

            /// Set up 123-126 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyyz_x = cbuffer.data(gp_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_z_yyyz_y = cbuffer.data(gp_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_z_yyyz_z = cbuffer.data(gp_geom_11_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yyyz_x, g_x_z_yyyz_y, g_x_z_yyyz_z, g_x_z_yyz_x, g_x_z_yyz_xy, g_x_z_yyz_y, g_x_z_yyz_yy, g_x_z_yyz_yz, g_x_z_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyyz_x[k] = -g_x_z_yyz_x[k] * ab_y + g_x_z_yyz_xy[k];

                g_x_z_yyyz_y[k] = -g_x_z_yyz_y[k] * ab_y + g_x_z_yyz_yy[k];

                g_x_z_yyyz_z[k] = -g_x_z_yyz_z[k] * ab_y + g_x_z_yyz_yz[k];
            }

            /// Set up 126-129 components of targeted buffer : cbuffer.data(

            auto g_x_z_yyzz_x = cbuffer.data(gp_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_z_yyzz_y = cbuffer.data(gp_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_z_yyzz_z = cbuffer.data(gp_geom_11_off + 128 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yyzz_x, g_x_z_yyzz_y, g_x_z_yyzz_z, g_x_z_yzz_x, g_x_z_yzz_xy, g_x_z_yzz_y, g_x_z_yzz_yy, g_x_z_yzz_yz, g_x_z_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yyzz_x[k] = -g_x_z_yzz_x[k] * ab_y + g_x_z_yzz_xy[k];

                g_x_z_yyzz_y[k] = -g_x_z_yzz_y[k] * ab_y + g_x_z_yzz_yy[k];

                g_x_z_yyzz_z[k] = -g_x_z_yzz_z[k] * ab_y + g_x_z_yzz_yz[k];
            }

            /// Set up 129-132 components of targeted buffer : cbuffer.data(

            auto g_x_z_yzzz_x = cbuffer.data(gp_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_z_yzzz_y = cbuffer.data(gp_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_z_yzzz_z = cbuffer.data(gp_geom_11_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yzzz_x, g_x_z_yzzz_y, g_x_z_yzzz_z, g_x_z_zzz_x, g_x_z_zzz_xy, g_x_z_zzz_y, g_x_z_zzz_yy, g_x_z_zzz_yz, g_x_z_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yzzz_x[k] = -g_x_z_zzz_x[k] * ab_y + g_x_z_zzz_xy[k];

                g_x_z_yzzz_y[k] = -g_x_z_zzz_y[k] * ab_y + g_x_z_zzz_yy[k];

                g_x_z_yzzz_z[k] = -g_x_z_zzz_z[k] * ab_y + g_x_z_zzz_yz[k];
            }

            /// Set up 132-135 components of targeted buffer : cbuffer.data(

            auto g_x_z_zzzz_x = cbuffer.data(gp_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_z_zzzz_y = cbuffer.data(gp_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_z_zzzz_z = cbuffer.data(gp_geom_11_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzz_x, g_x_0_zzz_y, g_x_0_zzz_z, g_x_z_zzz_x, g_x_z_zzz_xz, g_x_z_zzz_y, g_x_z_zzz_yz, g_x_z_zzz_z, g_x_z_zzz_zz, g_x_z_zzzz_x, g_x_z_zzzz_y, g_x_z_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_zzzz_x[k] = g_x_0_zzz_x[k] - g_x_z_zzz_x[k] * ab_z + g_x_z_zzz_xz[k];

                g_x_z_zzzz_y[k] = g_x_0_zzz_y[k] - g_x_z_zzz_y[k] * ab_z + g_x_z_zzz_yz[k];

                g_x_z_zzzz_z[k] = g_x_0_zzz_z[k] - g_x_z_zzz_z[k] * ab_z + g_x_z_zzz_zz[k];
            }

            /// Set up 135-138 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxxx_x = cbuffer.data(gp_geom_11_off + 135 * ccomps * dcomps);

            auto g_y_x_xxxx_y = cbuffer.data(gp_geom_11_off + 136 * ccomps * dcomps);

            auto g_y_x_xxxx_z = cbuffer.data(gp_geom_11_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxx_x, g_y_0_xxx_y, g_y_0_xxx_z, g_y_x_xxx_x, g_y_x_xxx_xx, g_y_x_xxx_xy, g_y_x_xxx_xz, g_y_x_xxx_y, g_y_x_xxx_z, g_y_x_xxxx_x, g_y_x_xxxx_y, g_y_x_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxxx_x[k] = g_y_0_xxx_x[k] - g_y_x_xxx_x[k] * ab_x + g_y_x_xxx_xx[k];

                g_y_x_xxxx_y[k] = g_y_0_xxx_y[k] - g_y_x_xxx_y[k] * ab_x + g_y_x_xxx_xy[k];

                g_y_x_xxxx_z[k] = g_y_0_xxx_z[k] - g_y_x_xxx_z[k] * ab_x + g_y_x_xxx_xz[k];
            }

            /// Set up 138-141 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxxy_x = cbuffer.data(gp_geom_11_off + 138 * ccomps * dcomps);

            auto g_y_x_xxxy_y = cbuffer.data(gp_geom_11_off + 139 * ccomps * dcomps);

            auto g_y_x_xxxy_z = cbuffer.data(gp_geom_11_off + 140 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxy_x, g_y_0_xxy_y, g_y_0_xxy_z, g_y_x_xxxy_x, g_y_x_xxxy_y, g_y_x_xxxy_z, g_y_x_xxy_x, g_y_x_xxy_xx, g_y_x_xxy_xy, g_y_x_xxy_xz, g_y_x_xxy_y, g_y_x_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxxy_x[k] = g_y_0_xxy_x[k] - g_y_x_xxy_x[k] * ab_x + g_y_x_xxy_xx[k];

                g_y_x_xxxy_y[k] = g_y_0_xxy_y[k] - g_y_x_xxy_y[k] * ab_x + g_y_x_xxy_xy[k];

                g_y_x_xxxy_z[k] = g_y_0_xxy_z[k] - g_y_x_xxy_z[k] * ab_x + g_y_x_xxy_xz[k];
            }

            /// Set up 141-144 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxxz_x = cbuffer.data(gp_geom_11_off + 141 * ccomps * dcomps);

            auto g_y_x_xxxz_y = cbuffer.data(gp_geom_11_off + 142 * ccomps * dcomps);

            auto g_y_x_xxxz_z = cbuffer.data(gp_geom_11_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xxx_x, g_y_x_xxx_xz, g_y_x_xxx_y, g_y_x_xxx_yz, g_y_x_xxx_z, g_y_x_xxx_zz, g_y_x_xxxz_x, g_y_x_xxxz_y, g_y_x_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxxz_x[k] = -g_y_x_xxx_x[k] * ab_z + g_y_x_xxx_xz[k];

                g_y_x_xxxz_y[k] = -g_y_x_xxx_y[k] * ab_z + g_y_x_xxx_yz[k];

                g_y_x_xxxz_z[k] = -g_y_x_xxx_z[k] * ab_z + g_y_x_xxx_zz[k];
            }

            /// Set up 144-147 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxyy_x = cbuffer.data(gp_geom_11_off + 144 * ccomps * dcomps);

            auto g_y_x_xxyy_y = cbuffer.data(gp_geom_11_off + 145 * ccomps * dcomps);

            auto g_y_x_xxyy_z = cbuffer.data(gp_geom_11_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyy_x, g_y_0_xyy_y, g_y_0_xyy_z, g_y_x_xxyy_x, g_y_x_xxyy_y, g_y_x_xxyy_z, g_y_x_xyy_x, g_y_x_xyy_xx, g_y_x_xyy_xy, g_y_x_xyy_xz, g_y_x_xyy_y, g_y_x_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxyy_x[k] = g_y_0_xyy_x[k] - g_y_x_xyy_x[k] * ab_x + g_y_x_xyy_xx[k];

                g_y_x_xxyy_y[k] = g_y_0_xyy_y[k] - g_y_x_xyy_y[k] * ab_x + g_y_x_xyy_xy[k];

                g_y_x_xxyy_z[k] = g_y_0_xyy_z[k] - g_y_x_xyy_z[k] * ab_x + g_y_x_xyy_xz[k];
            }

            /// Set up 147-150 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxyz_x = cbuffer.data(gp_geom_11_off + 147 * ccomps * dcomps);

            auto g_y_x_xxyz_y = cbuffer.data(gp_geom_11_off + 148 * ccomps * dcomps);

            auto g_y_x_xxyz_z = cbuffer.data(gp_geom_11_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xxy_x, g_y_x_xxy_xz, g_y_x_xxy_y, g_y_x_xxy_yz, g_y_x_xxy_z, g_y_x_xxy_zz, g_y_x_xxyz_x, g_y_x_xxyz_y, g_y_x_xxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxyz_x[k] = -g_y_x_xxy_x[k] * ab_z + g_y_x_xxy_xz[k];

                g_y_x_xxyz_y[k] = -g_y_x_xxy_y[k] * ab_z + g_y_x_xxy_yz[k];

                g_y_x_xxyz_z[k] = -g_y_x_xxy_z[k] * ab_z + g_y_x_xxy_zz[k];
            }

            /// Set up 150-153 components of targeted buffer : cbuffer.data(

            auto g_y_x_xxzz_x = cbuffer.data(gp_geom_11_off + 150 * ccomps * dcomps);

            auto g_y_x_xxzz_y = cbuffer.data(gp_geom_11_off + 151 * ccomps * dcomps);

            auto g_y_x_xxzz_z = cbuffer.data(gp_geom_11_off + 152 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xxz_x, g_y_x_xxz_xz, g_y_x_xxz_y, g_y_x_xxz_yz, g_y_x_xxz_z, g_y_x_xxz_zz, g_y_x_xxzz_x, g_y_x_xxzz_y, g_y_x_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xxzz_x[k] = -g_y_x_xxz_x[k] * ab_z + g_y_x_xxz_xz[k];

                g_y_x_xxzz_y[k] = -g_y_x_xxz_y[k] * ab_z + g_y_x_xxz_yz[k];

                g_y_x_xxzz_z[k] = -g_y_x_xxz_z[k] * ab_z + g_y_x_xxz_zz[k];
            }

            /// Set up 153-156 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyyy_x = cbuffer.data(gp_geom_11_off + 153 * ccomps * dcomps);

            auto g_y_x_xyyy_y = cbuffer.data(gp_geom_11_off + 154 * ccomps * dcomps);

            auto g_y_x_xyyy_z = cbuffer.data(gp_geom_11_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyy_x, g_y_0_yyy_y, g_y_0_yyy_z, g_y_x_xyyy_x, g_y_x_xyyy_y, g_y_x_xyyy_z, g_y_x_yyy_x, g_y_x_yyy_xx, g_y_x_yyy_xy, g_y_x_yyy_xz, g_y_x_yyy_y, g_y_x_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyyy_x[k] = g_y_0_yyy_x[k] - g_y_x_yyy_x[k] * ab_x + g_y_x_yyy_xx[k];

                g_y_x_xyyy_y[k] = g_y_0_yyy_y[k] - g_y_x_yyy_y[k] * ab_x + g_y_x_yyy_xy[k];

                g_y_x_xyyy_z[k] = g_y_0_yyy_z[k] - g_y_x_yyy_z[k] * ab_x + g_y_x_yyy_xz[k];
            }

            /// Set up 156-159 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyyz_x = cbuffer.data(gp_geom_11_off + 156 * ccomps * dcomps);

            auto g_y_x_xyyz_y = cbuffer.data(gp_geom_11_off + 157 * ccomps * dcomps);

            auto g_y_x_xyyz_z = cbuffer.data(gp_geom_11_off + 158 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xyy_x, g_y_x_xyy_xz, g_y_x_xyy_y, g_y_x_xyy_yz, g_y_x_xyy_z, g_y_x_xyy_zz, g_y_x_xyyz_x, g_y_x_xyyz_y, g_y_x_xyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyyz_x[k] = -g_y_x_xyy_x[k] * ab_z + g_y_x_xyy_xz[k];

                g_y_x_xyyz_y[k] = -g_y_x_xyy_y[k] * ab_z + g_y_x_xyy_yz[k];

                g_y_x_xyyz_z[k] = -g_y_x_xyy_z[k] * ab_z + g_y_x_xyy_zz[k];
            }

            /// Set up 159-162 components of targeted buffer : cbuffer.data(

            auto g_y_x_xyzz_x = cbuffer.data(gp_geom_11_off + 159 * ccomps * dcomps);

            auto g_y_x_xyzz_y = cbuffer.data(gp_geom_11_off + 160 * ccomps * dcomps);

            auto g_y_x_xyzz_z = cbuffer.data(gp_geom_11_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xyz_x, g_y_x_xyz_xz, g_y_x_xyz_y, g_y_x_xyz_yz, g_y_x_xyz_z, g_y_x_xyz_zz, g_y_x_xyzz_x, g_y_x_xyzz_y, g_y_x_xyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xyzz_x[k] = -g_y_x_xyz_x[k] * ab_z + g_y_x_xyz_xz[k];

                g_y_x_xyzz_y[k] = -g_y_x_xyz_y[k] * ab_z + g_y_x_xyz_yz[k];

                g_y_x_xyzz_z[k] = -g_y_x_xyz_z[k] * ab_z + g_y_x_xyz_zz[k];
            }

            /// Set up 162-165 components of targeted buffer : cbuffer.data(

            auto g_y_x_xzzz_x = cbuffer.data(gp_geom_11_off + 162 * ccomps * dcomps);

            auto g_y_x_xzzz_y = cbuffer.data(gp_geom_11_off + 163 * ccomps * dcomps);

            auto g_y_x_xzzz_z = cbuffer.data(gp_geom_11_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_xzz_x, g_y_x_xzz_xz, g_y_x_xzz_y, g_y_x_xzz_yz, g_y_x_xzz_z, g_y_x_xzz_zz, g_y_x_xzzz_x, g_y_x_xzzz_y, g_y_x_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xzzz_x[k] = -g_y_x_xzz_x[k] * ab_z + g_y_x_xzz_xz[k];

                g_y_x_xzzz_y[k] = -g_y_x_xzz_y[k] * ab_z + g_y_x_xzz_yz[k];

                g_y_x_xzzz_z[k] = -g_y_x_xzz_z[k] * ab_z + g_y_x_xzz_zz[k];
            }

            /// Set up 165-168 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyyy_x = cbuffer.data(gp_geom_11_off + 165 * ccomps * dcomps);

            auto g_y_x_yyyy_y = cbuffer.data(gp_geom_11_off + 166 * ccomps * dcomps);

            auto g_y_x_yyyy_z = cbuffer.data(gp_geom_11_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyy_x, g_0_x_yyy_y, g_0_x_yyy_z, g_y_x_yyy_x, g_y_x_yyy_xy, g_y_x_yyy_y, g_y_x_yyy_yy, g_y_x_yyy_yz, g_y_x_yyy_z, g_y_x_yyyy_x, g_y_x_yyyy_y, g_y_x_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyyy_x[k] = -g_0_x_yyy_x[k] - g_y_x_yyy_x[k] * ab_y + g_y_x_yyy_xy[k];

                g_y_x_yyyy_y[k] = -g_0_x_yyy_y[k] - g_y_x_yyy_y[k] * ab_y + g_y_x_yyy_yy[k];

                g_y_x_yyyy_z[k] = -g_0_x_yyy_z[k] - g_y_x_yyy_z[k] * ab_y + g_y_x_yyy_yz[k];
            }

            /// Set up 168-171 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyyz_x = cbuffer.data(gp_geom_11_off + 168 * ccomps * dcomps);

            auto g_y_x_yyyz_y = cbuffer.data(gp_geom_11_off + 169 * ccomps * dcomps);

            auto g_y_x_yyyz_z = cbuffer.data(gp_geom_11_off + 170 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yyy_x, g_y_x_yyy_xz, g_y_x_yyy_y, g_y_x_yyy_yz, g_y_x_yyy_z, g_y_x_yyy_zz, g_y_x_yyyz_x, g_y_x_yyyz_y, g_y_x_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyyz_x[k] = -g_y_x_yyy_x[k] * ab_z + g_y_x_yyy_xz[k];

                g_y_x_yyyz_y[k] = -g_y_x_yyy_y[k] * ab_z + g_y_x_yyy_yz[k];

                g_y_x_yyyz_z[k] = -g_y_x_yyy_z[k] * ab_z + g_y_x_yyy_zz[k];
            }

            /// Set up 171-174 components of targeted buffer : cbuffer.data(

            auto g_y_x_yyzz_x = cbuffer.data(gp_geom_11_off + 171 * ccomps * dcomps);

            auto g_y_x_yyzz_y = cbuffer.data(gp_geom_11_off + 172 * ccomps * dcomps);

            auto g_y_x_yyzz_z = cbuffer.data(gp_geom_11_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yyz_x, g_y_x_yyz_xz, g_y_x_yyz_y, g_y_x_yyz_yz, g_y_x_yyz_z, g_y_x_yyz_zz, g_y_x_yyzz_x, g_y_x_yyzz_y, g_y_x_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yyzz_x[k] = -g_y_x_yyz_x[k] * ab_z + g_y_x_yyz_xz[k];

                g_y_x_yyzz_y[k] = -g_y_x_yyz_y[k] * ab_z + g_y_x_yyz_yz[k];

                g_y_x_yyzz_z[k] = -g_y_x_yyz_z[k] * ab_z + g_y_x_yyz_zz[k];
            }

            /// Set up 174-177 components of targeted buffer : cbuffer.data(

            auto g_y_x_yzzz_x = cbuffer.data(gp_geom_11_off + 174 * ccomps * dcomps);

            auto g_y_x_yzzz_y = cbuffer.data(gp_geom_11_off + 175 * ccomps * dcomps);

            auto g_y_x_yzzz_z = cbuffer.data(gp_geom_11_off + 176 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_yzz_x, g_y_x_yzz_xz, g_y_x_yzz_y, g_y_x_yzz_yz, g_y_x_yzz_z, g_y_x_yzz_zz, g_y_x_yzzz_x, g_y_x_yzzz_y, g_y_x_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yzzz_x[k] = -g_y_x_yzz_x[k] * ab_z + g_y_x_yzz_xz[k];

                g_y_x_yzzz_y[k] = -g_y_x_yzz_y[k] * ab_z + g_y_x_yzz_yz[k];

                g_y_x_yzzz_z[k] = -g_y_x_yzz_z[k] * ab_z + g_y_x_yzz_zz[k];
            }

            /// Set up 177-180 components of targeted buffer : cbuffer.data(

            auto g_y_x_zzzz_x = cbuffer.data(gp_geom_11_off + 177 * ccomps * dcomps);

            auto g_y_x_zzzz_y = cbuffer.data(gp_geom_11_off + 178 * ccomps * dcomps);

            auto g_y_x_zzzz_z = cbuffer.data(gp_geom_11_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_zzz_x, g_y_x_zzz_xz, g_y_x_zzz_y, g_y_x_zzz_yz, g_y_x_zzz_z, g_y_x_zzz_zz, g_y_x_zzzz_x, g_y_x_zzzz_y, g_y_x_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_zzzz_x[k] = -g_y_x_zzz_x[k] * ab_z + g_y_x_zzz_xz[k];

                g_y_x_zzzz_y[k] = -g_y_x_zzz_y[k] * ab_z + g_y_x_zzz_yz[k];

                g_y_x_zzzz_z[k] = -g_y_x_zzz_z[k] * ab_z + g_y_x_zzz_zz[k];
            }

            /// Set up 180-183 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxxx_x = cbuffer.data(gp_geom_11_off + 180 * ccomps * dcomps);

            auto g_y_y_xxxx_y = cbuffer.data(gp_geom_11_off + 181 * ccomps * dcomps);

            auto g_y_y_xxxx_z = cbuffer.data(gp_geom_11_off + 182 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxx_x, g_y_y_xxx_xx, g_y_y_xxx_xy, g_y_y_xxx_xz, g_y_y_xxx_y, g_y_y_xxx_z, g_y_y_xxxx_x, g_y_y_xxxx_y, g_y_y_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxxx_x[k] = -g_y_y_xxx_x[k] * ab_x + g_y_y_xxx_xx[k];

                g_y_y_xxxx_y[k] = -g_y_y_xxx_y[k] * ab_x + g_y_y_xxx_xy[k];

                g_y_y_xxxx_z[k] = -g_y_y_xxx_z[k] * ab_x + g_y_y_xxx_xz[k];
            }

            /// Set up 183-186 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxxy_x = cbuffer.data(gp_geom_11_off + 183 * ccomps * dcomps);

            auto g_y_y_xxxy_y = cbuffer.data(gp_geom_11_off + 184 * ccomps * dcomps);

            auto g_y_y_xxxy_z = cbuffer.data(gp_geom_11_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxxy_x, g_y_y_xxxy_y, g_y_y_xxxy_z, g_y_y_xxy_x, g_y_y_xxy_xx, g_y_y_xxy_xy, g_y_y_xxy_xz, g_y_y_xxy_y, g_y_y_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxxy_x[k] = -g_y_y_xxy_x[k] * ab_x + g_y_y_xxy_xx[k];

                g_y_y_xxxy_y[k] = -g_y_y_xxy_y[k] * ab_x + g_y_y_xxy_xy[k];

                g_y_y_xxxy_z[k] = -g_y_y_xxy_z[k] * ab_x + g_y_y_xxy_xz[k];
            }

            /// Set up 186-189 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxxz_x = cbuffer.data(gp_geom_11_off + 186 * ccomps * dcomps);

            auto g_y_y_xxxz_y = cbuffer.data(gp_geom_11_off + 187 * ccomps * dcomps);

            auto g_y_y_xxxz_z = cbuffer.data(gp_geom_11_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxxz_x, g_y_y_xxxz_y, g_y_y_xxxz_z, g_y_y_xxz_x, g_y_y_xxz_xx, g_y_y_xxz_xy, g_y_y_xxz_xz, g_y_y_xxz_y, g_y_y_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxxz_x[k] = -g_y_y_xxz_x[k] * ab_x + g_y_y_xxz_xx[k];

                g_y_y_xxxz_y[k] = -g_y_y_xxz_y[k] * ab_x + g_y_y_xxz_xy[k];

                g_y_y_xxxz_z[k] = -g_y_y_xxz_z[k] * ab_x + g_y_y_xxz_xz[k];
            }

            /// Set up 189-192 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxyy_x = cbuffer.data(gp_geom_11_off + 189 * ccomps * dcomps);

            auto g_y_y_xxyy_y = cbuffer.data(gp_geom_11_off + 190 * ccomps * dcomps);

            auto g_y_y_xxyy_z = cbuffer.data(gp_geom_11_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxyy_x, g_y_y_xxyy_y, g_y_y_xxyy_z, g_y_y_xyy_x, g_y_y_xyy_xx, g_y_y_xyy_xy, g_y_y_xyy_xz, g_y_y_xyy_y, g_y_y_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxyy_x[k] = -g_y_y_xyy_x[k] * ab_x + g_y_y_xyy_xx[k];

                g_y_y_xxyy_y[k] = -g_y_y_xyy_y[k] * ab_x + g_y_y_xyy_xy[k];

                g_y_y_xxyy_z[k] = -g_y_y_xyy_z[k] * ab_x + g_y_y_xyy_xz[k];
            }

            /// Set up 192-195 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxyz_x = cbuffer.data(gp_geom_11_off + 192 * ccomps * dcomps);

            auto g_y_y_xxyz_y = cbuffer.data(gp_geom_11_off + 193 * ccomps * dcomps);

            auto g_y_y_xxyz_z = cbuffer.data(gp_geom_11_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxyz_x, g_y_y_xxyz_y, g_y_y_xxyz_z, g_y_y_xyz_x, g_y_y_xyz_xx, g_y_y_xyz_xy, g_y_y_xyz_xz, g_y_y_xyz_y, g_y_y_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxyz_x[k] = -g_y_y_xyz_x[k] * ab_x + g_y_y_xyz_xx[k];

                g_y_y_xxyz_y[k] = -g_y_y_xyz_y[k] * ab_x + g_y_y_xyz_xy[k];

                g_y_y_xxyz_z[k] = -g_y_y_xyz_z[k] * ab_x + g_y_y_xyz_xz[k];
            }

            /// Set up 195-198 components of targeted buffer : cbuffer.data(

            auto g_y_y_xxzz_x = cbuffer.data(gp_geom_11_off + 195 * ccomps * dcomps);

            auto g_y_y_xxzz_y = cbuffer.data(gp_geom_11_off + 196 * ccomps * dcomps);

            auto g_y_y_xxzz_z = cbuffer.data(gp_geom_11_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xxzz_x, g_y_y_xxzz_y, g_y_y_xxzz_z, g_y_y_xzz_x, g_y_y_xzz_xx, g_y_y_xzz_xy, g_y_y_xzz_xz, g_y_y_xzz_y, g_y_y_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xxzz_x[k] = -g_y_y_xzz_x[k] * ab_x + g_y_y_xzz_xx[k];

                g_y_y_xxzz_y[k] = -g_y_y_xzz_y[k] * ab_x + g_y_y_xzz_xy[k];

                g_y_y_xxzz_z[k] = -g_y_y_xzz_z[k] * ab_x + g_y_y_xzz_xz[k];
            }

            /// Set up 198-201 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyyy_x = cbuffer.data(gp_geom_11_off + 198 * ccomps * dcomps);

            auto g_y_y_xyyy_y = cbuffer.data(gp_geom_11_off + 199 * ccomps * dcomps);

            auto g_y_y_xyyy_z = cbuffer.data(gp_geom_11_off + 200 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyyy_x, g_y_y_xyyy_y, g_y_y_xyyy_z, g_y_y_yyy_x, g_y_y_yyy_xx, g_y_y_yyy_xy, g_y_y_yyy_xz, g_y_y_yyy_y, g_y_y_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyyy_x[k] = -g_y_y_yyy_x[k] * ab_x + g_y_y_yyy_xx[k];

                g_y_y_xyyy_y[k] = -g_y_y_yyy_y[k] * ab_x + g_y_y_yyy_xy[k];

                g_y_y_xyyy_z[k] = -g_y_y_yyy_z[k] * ab_x + g_y_y_yyy_xz[k];
            }

            /// Set up 201-204 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyyz_x = cbuffer.data(gp_geom_11_off + 201 * ccomps * dcomps);

            auto g_y_y_xyyz_y = cbuffer.data(gp_geom_11_off + 202 * ccomps * dcomps);

            auto g_y_y_xyyz_z = cbuffer.data(gp_geom_11_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyyz_x, g_y_y_xyyz_y, g_y_y_xyyz_z, g_y_y_yyz_x, g_y_y_yyz_xx, g_y_y_yyz_xy, g_y_y_yyz_xz, g_y_y_yyz_y, g_y_y_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyyz_x[k] = -g_y_y_yyz_x[k] * ab_x + g_y_y_yyz_xx[k];

                g_y_y_xyyz_y[k] = -g_y_y_yyz_y[k] * ab_x + g_y_y_yyz_xy[k];

                g_y_y_xyyz_z[k] = -g_y_y_yyz_z[k] * ab_x + g_y_y_yyz_xz[k];
            }

            /// Set up 204-207 components of targeted buffer : cbuffer.data(

            auto g_y_y_xyzz_x = cbuffer.data(gp_geom_11_off + 204 * ccomps * dcomps);

            auto g_y_y_xyzz_y = cbuffer.data(gp_geom_11_off + 205 * ccomps * dcomps);

            auto g_y_y_xyzz_z = cbuffer.data(gp_geom_11_off + 206 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xyzz_x, g_y_y_xyzz_y, g_y_y_xyzz_z, g_y_y_yzz_x, g_y_y_yzz_xx, g_y_y_yzz_xy, g_y_y_yzz_xz, g_y_y_yzz_y, g_y_y_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xyzz_x[k] = -g_y_y_yzz_x[k] * ab_x + g_y_y_yzz_xx[k];

                g_y_y_xyzz_y[k] = -g_y_y_yzz_y[k] * ab_x + g_y_y_yzz_xy[k];

                g_y_y_xyzz_z[k] = -g_y_y_yzz_z[k] * ab_x + g_y_y_yzz_xz[k];
            }

            /// Set up 207-210 components of targeted buffer : cbuffer.data(

            auto g_y_y_xzzz_x = cbuffer.data(gp_geom_11_off + 207 * ccomps * dcomps);

            auto g_y_y_xzzz_y = cbuffer.data(gp_geom_11_off + 208 * ccomps * dcomps);

            auto g_y_y_xzzz_z = cbuffer.data(gp_geom_11_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xzzz_x, g_y_y_xzzz_y, g_y_y_xzzz_z, g_y_y_zzz_x, g_y_y_zzz_xx, g_y_y_zzz_xy, g_y_y_zzz_xz, g_y_y_zzz_y, g_y_y_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xzzz_x[k] = -g_y_y_zzz_x[k] * ab_x + g_y_y_zzz_xx[k];

                g_y_y_xzzz_y[k] = -g_y_y_zzz_y[k] * ab_x + g_y_y_zzz_xy[k];

                g_y_y_xzzz_z[k] = -g_y_y_zzz_z[k] * ab_x + g_y_y_zzz_xz[k];
            }

            /// Set up 210-213 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyyy_x = cbuffer.data(gp_geom_11_off + 210 * ccomps * dcomps);

            auto g_y_y_yyyy_y = cbuffer.data(gp_geom_11_off + 211 * ccomps * dcomps);

            auto g_y_y_yyyy_z = cbuffer.data(gp_geom_11_off + 212 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyy_x, g_0_y_yyy_y, g_0_y_yyy_z, g_y_0_yyy_x, g_y_0_yyy_y, g_y_0_yyy_z, g_y_y_yyy_x, g_y_y_yyy_xy, g_y_y_yyy_y, g_y_y_yyy_yy, g_y_y_yyy_yz, g_y_y_yyy_z, g_y_y_yyyy_x, g_y_y_yyyy_y, g_y_y_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyyy_x[k] = -g_0_y_yyy_x[k] + g_y_0_yyy_x[k] - g_y_y_yyy_x[k] * ab_y + g_y_y_yyy_xy[k];

                g_y_y_yyyy_y[k] = -g_0_y_yyy_y[k] + g_y_0_yyy_y[k] - g_y_y_yyy_y[k] * ab_y + g_y_y_yyy_yy[k];

                g_y_y_yyyy_z[k] = -g_0_y_yyy_z[k] + g_y_0_yyy_z[k] - g_y_y_yyy_z[k] * ab_y + g_y_y_yyy_yz[k];
            }

            /// Set up 213-216 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyyz_x = cbuffer.data(gp_geom_11_off + 213 * ccomps * dcomps);

            auto g_y_y_yyyz_y = cbuffer.data(gp_geom_11_off + 214 * ccomps * dcomps);

            auto g_y_y_yyyz_z = cbuffer.data(gp_geom_11_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yyy_x, g_y_y_yyy_xz, g_y_y_yyy_y, g_y_y_yyy_yz, g_y_y_yyy_z, g_y_y_yyy_zz, g_y_y_yyyz_x, g_y_y_yyyz_y, g_y_y_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyyz_x[k] = -g_y_y_yyy_x[k] * ab_z + g_y_y_yyy_xz[k];

                g_y_y_yyyz_y[k] = -g_y_y_yyy_y[k] * ab_z + g_y_y_yyy_yz[k];

                g_y_y_yyyz_z[k] = -g_y_y_yyy_z[k] * ab_z + g_y_y_yyy_zz[k];
            }

            /// Set up 216-219 components of targeted buffer : cbuffer.data(

            auto g_y_y_yyzz_x = cbuffer.data(gp_geom_11_off + 216 * ccomps * dcomps);

            auto g_y_y_yyzz_y = cbuffer.data(gp_geom_11_off + 217 * ccomps * dcomps);

            auto g_y_y_yyzz_z = cbuffer.data(gp_geom_11_off + 218 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yyz_x, g_y_y_yyz_xz, g_y_y_yyz_y, g_y_y_yyz_yz, g_y_y_yyz_z, g_y_y_yyz_zz, g_y_y_yyzz_x, g_y_y_yyzz_y, g_y_y_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yyzz_x[k] = -g_y_y_yyz_x[k] * ab_z + g_y_y_yyz_xz[k];

                g_y_y_yyzz_y[k] = -g_y_y_yyz_y[k] * ab_z + g_y_y_yyz_yz[k];

                g_y_y_yyzz_z[k] = -g_y_y_yyz_z[k] * ab_z + g_y_y_yyz_zz[k];
            }

            /// Set up 219-222 components of targeted buffer : cbuffer.data(

            auto g_y_y_yzzz_x = cbuffer.data(gp_geom_11_off + 219 * ccomps * dcomps);

            auto g_y_y_yzzz_y = cbuffer.data(gp_geom_11_off + 220 * ccomps * dcomps);

            auto g_y_y_yzzz_z = cbuffer.data(gp_geom_11_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_yzz_x, g_y_y_yzz_xz, g_y_y_yzz_y, g_y_y_yzz_yz, g_y_y_yzz_z, g_y_y_yzz_zz, g_y_y_yzzz_x, g_y_y_yzzz_y, g_y_y_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yzzz_x[k] = -g_y_y_yzz_x[k] * ab_z + g_y_y_yzz_xz[k];

                g_y_y_yzzz_y[k] = -g_y_y_yzz_y[k] * ab_z + g_y_y_yzz_yz[k];

                g_y_y_yzzz_z[k] = -g_y_y_yzz_z[k] * ab_z + g_y_y_yzz_zz[k];
            }

            /// Set up 222-225 components of targeted buffer : cbuffer.data(

            auto g_y_y_zzzz_x = cbuffer.data(gp_geom_11_off + 222 * ccomps * dcomps);

            auto g_y_y_zzzz_y = cbuffer.data(gp_geom_11_off + 223 * ccomps * dcomps);

            auto g_y_y_zzzz_z = cbuffer.data(gp_geom_11_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_zzz_x, g_y_y_zzz_xz, g_y_y_zzz_y, g_y_y_zzz_yz, g_y_y_zzz_z, g_y_y_zzz_zz, g_y_y_zzzz_x, g_y_y_zzzz_y, g_y_y_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_zzzz_x[k] = -g_y_y_zzz_x[k] * ab_z + g_y_y_zzz_xz[k];

                g_y_y_zzzz_y[k] = -g_y_y_zzz_y[k] * ab_z + g_y_y_zzz_yz[k];

                g_y_y_zzzz_z[k] = -g_y_y_zzz_z[k] * ab_z + g_y_y_zzz_zz[k];
            }

            /// Set up 225-228 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxxx_x = cbuffer.data(gp_geom_11_off + 225 * ccomps * dcomps);

            auto g_y_z_xxxx_y = cbuffer.data(gp_geom_11_off + 226 * ccomps * dcomps);

            auto g_y_z_xxxx_z = cbuffer.data(gp_geom_11_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxx_x, g_y_z_xxx_xx, g_y_z_xxx_xy, g_y_z_xxx_xz, g_y_z_xxx_y, g_y_z_xxx_z, g_y_z_xxxx_x, g_y_z_xxxx_y, g_y_z_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxxx_x[k] = -g_y_z_xxx_x[k] * ab_x + g_y_z_xxx_xx[k];

                g_y_z_xxxx_y[k] = -g_y_z_xxx_y[k] * ab_x + g_y_z_xxx_xy[k];

                g_y_z_xxxx_z[k] = -g_y_z_xxx_z[k] * ab_x + g_y_z_xxx_xz[k];
            }

            /// Set up 228-231 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxxy_x = cbuffer.data(gp_geom_11_off + 228 * ccomps * dcomps);

            auto g_y_z_xxxy_y = cbuffer.data(gp_geom_11_off + 229 * ccomps * dcomps);

            auto g_y_z_xxxy_z = cbuffer.data(gp_geom_11_off + 230 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxxy_x, g_y_z_xxxy_y, g_y_z_xxxy_z, g_y_z_xxy_x, g_y_z_xxy_xx, g_y_z_xxy_xy, g_y_z_xxy_xz, g_y_z_xxy_y, g_y_z_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxxy_x[k] = -g_y_z_xxy_x[k] * ab_x + g_y_z_xxy_xx[k];

                g_y_z_xxxy_y[k] = -g_y_z_xxy_y[k] * ab_x + g_y_z_xxy_xy[k];

                g_y_z_xxxy_z[k] = -g_y_z_xxy_z[k] * ab_x + g_y_z_xxy_xz[k];
            }

            /// Set up 231-234 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxxz_x = cbuffer.data(gp_geom_11_off + 231 * ccomps * dcomps);

            auto g_y_z_xxxz_y = cbuffer.data(gp_geom_11_off + 232 * ccomps * dcomps);

            auto g_y_z_xxxz_z = cbuffer.data(gp_geom_11_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxxz_x, g_y_z_xxxz_y, g_y_z_xxxz_z, g_y_z_xxz_x, g_y_z_xxz_xx, g_y_z_xxz_xy, g_y_z_xxz_xz, g_y_z_xxz_y, g_y_z_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxxz_x[k] = -g_y_z_xxz_x[k] * ab_x + g_y_z_xxz_xx[k];

                g_y_z_xxxz_y[k] = -g_y_z_xxz_y[k] * ab_x + g_y_z_xxz_xy[k];

                g_y_z_xxxz_z[k] = -g_y_z_xxz_z[k] * ab_x + g_y_z_xxz_xz[k];
            }

            /// Set up 234-237 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxyy_x = cbuffer.data(gp_geom_11_off + 234 * ccomps * dcomps);

            auto g_y_z_xxyy_y = cbuffer.data(gp_geom_11_off + 235 * ccomps * dcomps);

            auto g_y_z_xxyy_z = cbuffer.data(gp_geom_11_off + 236 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxyy_x, g_y_z_xxyy_y, g_y_z_xxyy_z, g_y_z_xyy_x, g_y_z_xyy_xx, g_y_z_xyy_xy, g_y_z_xyy_xz, g_y_z_xyy_y, g_y_z_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxyy_x[k] = -g_y_z_xyy_x[k] * ab_x + g_y_z_xyy_xx[k];

                g_y_z_xxyy_y[k] = -g_y_z_xyy_y[k] * ab_x + g_y_z_xyy_xy[k];

                g_y_z_xxyy_z[k] = -g_y_z_xyy_z[k] * ab_x + g_y_z_xyy_xz[k];
            }

            /// Set up 237-240 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxyz_x = cbuffer.data(gp_geom_11_off + 237 * ccomps * dcomps);

            auto g_y_z_xxyz_y = cbuffer.data(gp_geom_11_off + 238 * ccomps * dcomps);

            auto g_y_z_xxyz_z = cbuffer.data(gp_geom_11_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxyz_x, g_y_z_xxyz_y, g_y_z_xxyz_z, g_y_z_xyz_x, g_y_z_xyz_xx, g_y_z_xyz_xy, g_y_z_xyz_xz, g_y_z_xyz_y, g_y_z_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxyz_x[k] = -g_y_z_xyz_x[k] * ab_x + g_y_z_xyz_xx[k];

                g_y_z_xxyz_y[k] = -g_y_z_xyz_y[k] * ab_x + g_y_z_xyz_xy[k];

                g_y_z_xxyz_z[k] = -g_y_z_xyz_z[k] * ab_x + g_y_z_xyz_xz[k];
            }

            /// Set up 240-243 components of targeted buffer : cbuffer.data(

            auto g_y_z_xxzz_x = cbuffer.data(gp_geom_11_off + 240 * ccomps * dcomps);

            auto g_y_z_xxzz_y = cbuffer.data(gp_geom_11_off + 241 * ccomps * dcomps);

            auto g_y_z_xxzz_z = cbuffer.data(gp_geom_11_off + 242 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xxzz_x, g_y_z_xxzz_y, g_y_z_xxzz_z, g_y_z_xzz_x, g_y_z_xzz_xx, g_y_z_xzz_xy, g_y_z_xzz_xz, g_y_z_xzz_y, g_y_z_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xxzz_x[k] = -g_y_z_xzz_x[k] * ab_x + g_y_z_xzz_xx[k];

                g_y_z_xxzz_y[k] = -g_y_z_xzz_y[k] * ab_x + g_y_z_xzz_xy[k];

                g_y_z_xxzz_z[k] = -g_y_z_xzz_z[k] * ab_x + g_y_z_xzz_xz[k];
            }

            /// Set up 243-246 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyyy_x = cbuffer.data(gp_geom_11_off + 243 * ccomps * dcomps);

            auto g_y_z_xyyy_y = cbuffer.data(gp_geom_11_off + 244 * ccomps * dcomps);

            auto g_y_z_xyyy_z = cbuffer.data(gp_geom_11_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyyy_x, g_y_z_xyyy_y, g_y_z_xyyy_z, g_y_z_yyy_x, g_y_z_yyy_xx, g_y_z_yyy_xy, g_y_z_yyy_xz, g_y_z_yyy_y, g_y_z_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyyy_x[k] = -g_y_z_yyy_x[k] * ab_x + g_y_z_yyy_xx[k];

                g_y_z_xyyy_y[k] = -g_y_z_yyy_y[k] * ab_x + g_y_z_yyy_xy[k];

                g_y_z_xyyy_z[k] = -g_y_z_yyy_z[k] * ab_x + g_y_z_yyy_xz[k];
            }

            /// Set up 246-249 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyyz_x = cbuffer.data(gp_geom_11_off + 246 * ccomps * dcomps);

            auto g_y_z_xyyz_y = cbuffer.data(gp_geom_11_off + 247 * ccomps * dcomps);

            auto g_y_z_xyyz_z = cbuffer.data(gp_geom_11_off + 248 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyyz_x, g_y_z_xyyz_y, g_y_z_xyyz_z, g_y_z_yyz_x, g_y_z_yyz_xx, g_y_z_yyz_xy, g_y_z_yyz_xz, g_y_z_yyz_y, g_y_z_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyyz_x[k] = -g_y_z_yyz_x[k] * ab_x + g_y_z_yyz_xx[k];

                g_y_z_xyyz_y[k] = -g_y_z_yyz_y[k] * ab_x + g_y_z_yyz_xy[k];

                g_y_z_xyyz_z[k] = -g_y_z_yyz_z[k] * ab_x + g_y_z_yyz_xz[k];
            }

            /// Set up 249-252 components of targeted buffer : cbuffer.data(

            auto g_y_z_xyzz_x = cbuffer.data(gp_geom_11_off + 249 * ccomps * dcomps);

            auto g_y_z_xyzz_y = cbuffer.data(gp_geom_11_off + 250 * ccomps * dcomps);

            auto g_y_z_xyzz_z = cbuffer.data(gp_geom_11_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xyzz_x, g_y_z_xyzz_y, g_y_z_xyzz_z, g_y_z_yzz_x, g_y_z_yzz_xx, g_y_z_yzz_xy, g_y_z_yzz_xz, g_y_z_yzz_y, g_y_z_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xyzz_x[k] = -g_y_z_yzz_x[k] * ab_x + g_y_z_yzz_xx[k];

                g_y_z_xyzz_y[k] = -g_y_z_yzz_y[k] * ab_x + g_y_z_yzz_xy[k];

                g_y_z_xyzz_z[k] = -g_y_z_yzz_z[k] * ab_x + g_y_z_yzz_xz[k];
            }

            /// Set up 252-255 components of targeted buffer : cbuffer.data(

            auto g_y_z_xzzz_x = cbuffer.data(gp_geom_11_off + 252 * ccomps * dcomps);

            auto g_y_z_xzzz_y = cbuffer.data(gp_geom_11_off + 253 * ccomps * dcomps);

            auto g_y_z_xzzz_z = cbuffer.data(gp_geom_11_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xzzz_x, g_y_z_xzzz_y, g_y_z_xzzz_z, g_y_z_zzz_x, g_y_z_zzz_xx, g_y_z_zzz_xy, g_y_z_zzz_xz, g_y_z_zzz_y, g_y_z_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xzzz_x[k] = -g_y_z_zzz_x[k] * ab_x + g_y_z_zzz_xx[k];

                g_y_z_xzzz_y[k] = -g_y_z_zzz_y[k] * ab_x + g_y_z_zzz_xy[k];

                g_y_z_xzzz_z[k] = -g_y_z_zzz_z[k] * ab_x + g_y_z_zzz_xz[k];
            }

            /// Set up 255-258 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyyy_x = cbuffer.data(gp_geom_11_off + 255 * ccomps * dcomps);

            auto g_y_z_yyyy_y = cbuffer.data(gp_geom_11_off + 256 * ccomps * dcomps);

            auto g_y_z_yyyy_z = cbuffer.data(gp_geom_11_off + 257 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyy_x, g_0_z_yyy_y, g_0_z_yyy_z, g_y_z_yyy_x, g_y_z_yyy_xy, g_y_z_yyy_y, g_y_z_yyy_yy, g_y_z_yyy_yz, g_y_z_yyy_z, g_y_z_yyyy_x, g_y_z_yyyy_y, g_y_z_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyyy_x[k] = -g_0_z_yyy_x[k] - g_y_z_yyy_x[k] * ab_y + g_y_z_yyy_xy[k];

                g_y_z_yyyy_y[k] = -g_0_z_yyy_y[k] - g_y_z_yyy_y[k] * ab_y + g_y_z_yyy_yy[k];

                g_y_z_yyyy_z[k] = -g_0_z_yyy_z[k] - g_y_z_yyy_z[k] * ab_y + g_y_z_yyy_yz[k];
            }

            /// Set up 258-261 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyyz_x = cbuffer.data(gp_geom_11_off + 258 * ccomps * dcomps);

            auto g_y_z_yyyz_y = cbuffer.data(gp_geom_11_off + 259 * ccomps * dcomps);

            auto g_y_z_yyyz_z = cbuffer.data(gp_geom_11_off + 260 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyz_x, g_0_z_yyz_y, g_0_z_yyz_z, g_y_z_yyyz_x, g_y_z_yyyz_y, g_y_z_yyyz_z, g_y_z_yyz_x, g_y_z_yyz_xy, g_y_z_yyz_y, g_y_z_yyz_yy, g_y_z_yyz_yz, g_y_z_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyyz_x[k] = -g_0_z_yyz_x[k] - g_y_z_yyz_x[k] * ab_y + g_y_z_yyz_xy[k];

                g_y_z_yyyz_y[k] = -g_0_z_yyz_y[k] - g_y_z_yyz_y[k] * ab_y + g_y_z_yyz_yy[k];

                g_y_z_yyyz_z[k] = -g_0_z_yyz_z[k] - g_y_z_yyz_z[k] * ab_y + g_y_z_yyz_yz[k];
            }

            /// Set up 261-264 components of targeted buffer : cbuffer.data(

            auto g_y_z_yyzz_x = cbuffer.data(gp_geom_11_off + 261 * ccomps * dcomps);

            auto g_y_z_yyzz_y = cbuffer.data(gp_geom_11_off + 262 * ccomps * dcomps);

            auto g_y_z_yyzz_z = cbuffer.data(gp_geom_11_off + 263 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzz_x, g_0_z_yzz_y, g_0_z_yzz_z, g_y_z_yyzz_x, g_y_z_yyzz_y, g_y_z_yyzz_z, g_y_z_yzz_x, g_y_z_yzz_xy, g_y_z_yzz_y, g_y_z_yzz_yy, g_y_z_yzz_yz, g_y_z_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yyzz_x[k] = -g_0_z_yzz_x[k] - g_y_z_yzz_x[k] * ab_y + g_y_z_yzz_xy[k];

                g_y_z_yyzz_y[k] = -g_0_z_yzz_y[k] - g_y_z_yzz_y[k] * ab_y + g_y_z_yzz_yy[k];

                g_y_z_yyzz_z[k] = -g_0_z_yzz_z[k] - g_y_z_yzz_z[k] * ab_y + g_y_z_yzz_yz[k];
            }

            /// Set up 264-267 components of targeted buffer : cbuffer.data(

            auto g_y_z_yzzz_x = cbuffer.data(gp_geom_11_off + 264 * ccomps * dcomps);

            auto g_y_z_yzzz_y = cbuffer.data(gp_geom_11_off + 265 * ccomps * dcomps);

            auto g_y_z_yzzz_z = cbuffer.data(gp_geom_11_off + 266 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzz_x, g_0_z_zzz_y, g_0_z_zzz_z, g_y_z_yzzz_x, g_y_z_yzzz_y, g_y_z_yzzz_z, g_y_z_zzz_x, g_y_z_zzz_xy, g_y_z_zzz_y, g_y_z_zzz_yy, g_y_z_zzz_yz, g_y_z_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yzzz_x[k] = -g_0_z_zzz_x[k] - g_y_z_zzz_x[k] * ab_y + g_y_z_zzz_xy[k];

                g_y_z_yzzz_y[k] = -g_0_z_zzz_y[k] - g_y_z_zzz_y[k] * ab_y + g_y_z_zzz_yy[k];

                g_y_z_yzzz_z[k] = -g_0_z_zzz_z[k] - g_y_z_zzz_z[k] * ab_y + g_y_z_zzz_yz[k];
            }

            /// Set up 267-270 components of targeted buffer : cbuffer.data(

            auto g_y_z_zzzz_x = cbuffer.data(gp_geom_11_off + 267 * ccomps * dcomps);

            auto g_y_z_zzzz_y = cbuffer.data(gp_geom_11_off + 268 * ccomps * dcomps);

            auto g_y_z_zzzz_z = cbuffer.data(gp_geom_11_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzz_x, g_y_0_zzz_y, g_y_0_zzz_z, g_y_z_zzz_x, g_y_z_zzz_xz, g_y_z_zzz_y, g_y_z_zzz_yz, g_y_z_zzz_z, g_y_z_zzz_zz, g_y_z_zzzz_x, g_y_z_zzzz_y, g_y_z_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_zzzz_x[k] = g_y_0_zzz_x[k] - g_y_z_zzz_x[k] * ab_z + g_y_z_zzz_xz[k];

                g_y_z_zzzz_y[k] = g_y_0_zzz_y[k] - g_y_z_zzz_y[k] * ab_z + g_y_z_zzz_yz[k];

                g_y_z_zzzz_z[k] = g_y_0_zzz_z[k] - g_y_z_zzz_z[k] * ab_z + g_y_z_zzz_zz[k];
            }

            /// Set up 270-273 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxxx_x = cbuffer.data(gp_geom_11_off + 270 * ccomps * dcomps);

            auto g_z_x_xxxx_y = cbuffer.data(gp_geom_11_off + 271 * ccomps * dcomps);

            auto g_z_x_xxxx_z = cbuffer.data(gp_geom_11_off + 272 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxx_x, g_z_0_xxx_y, g_z_0_xxx_z, g_z_x_xxx_x, g_z_x_xxx_xx, g_z_x_xxx_xy, g_z_x_xxx_xz, g_z_x_xxx_y, g_z_x_xxx_z, g_z_x_xxxx_x, g_z_x_xxxx_y, g_z_x_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxxx_x[k] = g_z_0_xxx_x[k] - g_z_x_xxx_x[k] * ab_x + g_z_x_xxx_xx[k];

                g_z_x_xxxx_y[k] = g_z_0_xxx_y[k] - g_z_x_xxx_y[k] * ab_x + g_z_x_xxx_xy[k];

                g_z_x_xxxx_z[k] = g_z_0_xxx_z[k] - g_z_x_xxx_z[k] * ab_x + g_z_x_xxx_xz[k];
            }

            /// Set up 273-276 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxxy_x = cbuffer.data(gp_geom_11_off + 273 * ccomps * dcomps);

            auto g_z_x_xxxy_y = cbuffer.data(gp_geom_11_off + 274 * ccomps * dcomps);

            auto g_z_x_xxxy_z = cbuffer.data(gp_geom_11_off + 275 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xxx_x, g_z_x_xxx_xy, g_z_x_xxx_y, g_z_x_xxx_yy, g_z_x_xxx_yz, g_z_x_xxx_z, g_z_x_xxxy_x, g_z_x_xxxy_y, g_z_x_xxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxxy_x[k] = -g_z_x_xxx_x[k] * ab_y + g_z_x_xxx_xy[k];

                g_z_x_xxxy_y[k] = -g_z_x_xxx_y[k] * ab_y + g_z_x_xxx_yy[k];

                g_z_x_xxxy_z[k] = -g_z_x_xxx_z[k] * ab_y + g_z_x_xxx_yz[k];
            }

            /// Set up 276-279 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxxz_x = cbuffer.data(gp_geom_11_off + 276 * ccomps * dcomps);

            auto g_z_x_xxxz_y = cbuffer.data(gp_geom_11_off + 277 * ccomps * dcomps);

            auto g_z_x_xxxz_z = cbuffer.data(gp_geom_11_off + 278 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxz_x, g_z_0_xxz_y, g_z_0_xxz_z, g_z_x_xxxz_x, g_z_x_xxxz_y, g_z_x_xxxz_z, g_z_x_xxz_x, g_z_x_xxz_xx, g_z_x_xxz_xy, g_z_x_xxz_xz, g_z_x_xxz_y, g_z_x_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxxz_x[k] = g_z_0_xxz_x[k] - g_z_x_xxz_x[k] * ab_x + g_z_x_xxz_xx[k];

                g_z_x_xxxz_y[k] = g_z_0_xxz_y[k] - g_z_x_xxz_y[k] * ab_x + g_z_x_xxz_xy[k];

                g_z_x_xxxz_z[k] = g_z_0_xxz_z[k] - g_z_x_xxz_z[k] * ab_x + g_z_x_xxz_xz[k];
            }

            /// Set up 279-282 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxyy_x = cbuffer.data(gp_geom_11_off + 279 * ccomps * dcomps);

            auto g_z_x_xxyy_y = cbuffer.data(gp_geom_11_off + 280 * ccomps * dcomps);

            auto g_z_x_xxyy_z = cbuffer.data(gp_geom_11_off + 281 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xxy_x, g_z_x_xxy_xy, g_z_x_xxy_y, g_z_x_xxy_yy, g_z_x_xxy_yz, g_z_x_xxy_z, g_z_x_xxyy_x, g_z_x_xxyy_y, g_z_x_xxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxyy_x[k] = -g_z_x_xxy_x[k] * ab_y + g_z_x_xxy_xy[k];

                g_z_x_xxyy_y[k] = -g_z_x_xxy_y[k] * ab_y + g_z_x_xxy_yy[k];

                g_z_x_xxyy_z[k] = -g_z_x_xxy_z[k] * ab_y + g_z_x_xxy_yz[k];
            }

            /// Set up 282-285 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxyz_x = cbuffer.data(gp_geom_11_off + 282 * ccomps * dcomps);

            auto g_z_x_xxyz_y = cbuffer.data(gp_geom_11_off + 283 * ccomps * dcomps);

            auto g_z_x_xxyz_z = cbuffer.data(gp_geom_11_off + 284 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xxyz_x, g_z_x_xxyz_y, g_z_x_xxyz_z, g_z_x_xxz_x, g_z_x_xxz_xy, g_z_x_xxz_y, g_z_x_xxz_yy, g_z_x_xxz_yz, g_z_x_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxyz_x[k] = -g_z_x_xxz_x[k] * ab_y + g_z_x_xxz_xy[k];

                g_z_x_xxyz_y[k] = -g_z_x_xxz_y[k] * ab_y + g_z_x_xxz_yy[k];

                g_z_x_xxyz_z[k] = -g_z_x_xxz_z[k] * ab_y + g_z_x_xxz_yz[k];
            }

            /// Set up 285-288 components of targeted buffer : cbuffer.data(

            auto g_z_x_xxzz_x = cbuffer.data(gp_geom_11_off + 285 * ccomps * dcomps);

            auto g_z_x_xxzz_y = cbuffer.data(gp_geom_11_off + 286 * ccomps * dcomps);

            auto g_z_x_xxzz_z = cbuffer.data(gp_geom_11_off + 287 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xzz_x, g_z_0_xzz_y, g_z_0_xzz_z, g_z_x_xxzz_x, g_z_x_xxzz_y, g_z_x_xxzz_z, g_z_x_xzz_x, g_z_x_xzz_xx, g_z_x_xzz_xy, g_z_x_xzz_xz, g_z_x_xzz_y, g_z_x_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xxzz_x[k] = g_z_0_xzz_x[k] - g_z_x_xzz_x[k] * ab_x + g_z_x_xzz_xx[k];

                g_z_x_xxzz_y[k] = g_z_0_xzz_y[k] - g_z_x_xzz_y[k] * ab_x + g_z_x_xzz_xy[k];

                g_z_x_xxzz_z[k] = g_z_0_xzz_z[k] - g_z_x_xzz_z[k] * ab_x + g_z_x_xzz_xz[k];
            }

            /// Set up 288-291 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyyy_x = cbuffer.data(gp_geom_11_off + 288 * ccomps * dcomps);

            auto g_z_x_xyyy_y = cbuffer.data(gp_geom_11_off + 289 * ccomps * dcomps);

            auto g_z_x_xyyy_z = cbuffer.data(gp_geom_11_off + 290 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xyy_x, g_z_x_xyy_xy, g_z_x_xyy_y, g_z_x_xyy_yy, g_z_x_xyy_yz, g_z_x_xyy_z, g_z_x_xyyy_x, g_z_x_xyyy_y, g_z_x_xyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyyy_x[k] = -g_z_x_xyy_x[k] * ab_y + g_z_x_xyy_xy[k];

                g_z_x_xyyy_y[k] = -g_z_x_xyy_y[k] * ab_y + g_z_x_xyy_yy[k];

                g_z_x_xyyy_z[k] = -g_z_x_xyy_z[k] * ab_y + g_z_x_xyy_yz[k];
            }

            /// Set up 291-294 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyyz_x = cbuffer.data(gp_geom_11_off + 291 * ccomps * dcomps);

            auto g_z_x_xyyz_y = cbuffer.data(gp_geom_11_off + 292 * ccomps * dcomps);

            auto g_z_x_xyyz_z = cbuffer.data(gp_geom_11_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xyyz_x, g_z_x_xyyz_y, g_z_x_xyyz_z, g_z_x_xyz_x, g_z_x_xyz_xy, g_z_x_xyz_y, g_z_x_xyz_yy, g_z_x_xyz_yz, g_z_x_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyyz_x[k] = -g_z_x_xyz_x[k] * ab_y + g_z_x_xyz_xy[k];

                g_z_x_xyyz_y[k] = -g_z_x_xyz_y[k] * ab_y + g_z_x_xyz_yy[k];

                g_z_x_xyyz_z[k] = -g_z_x_xyz_z[k] * ab_y + g_z_x_xyz_yz[k];
            }

            /// Set up 294-297 components of targeted buffer : cbuffer.data(

            auto g_z_x_xyzz_x = cbuffer.data(gp_geom_11_off + 294 * ccomps * dcomps);

            auto g_z_x_xyzz_y = cbuffer.data(gp_geom_11_off + 295 * ccomps * dcomps);

            auto g_z_x_xyzz_z = cbuffer.data(gp_geom_11_off + 296 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_xyzz_x, g_z_x_xyzz_y, g_z_x_xyzz_z, g_z_x_xzz_x, g_z_x_xzz_xy, g_z_x_xzz_y, g_z_x_xzz_yy, g_z_x_xzz_yz, g_z_x_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xyzz_x[k] = -g_z_x_xzz_x[k] * ab_y + g_z_x_xzz_xy[k];

                g_z_x_xyzz_y[k] = -g_z_x_xzz_y[k] * ab_y + g_z_x_xzz_yy[k];

                g_z_x_xyzz_z[k] = -g_z_x_xzz_z[k] * ab_y + g_z_x_xzz_yz[k];
            }

            /// Set up 297-300 components of targeted buffer : cbuffer.data(

            auto g_z_x_xzzz_x = cbuffer.data(gp_geom_11_off + 297 * ccomps * dcomps);

            auto g_z_x_xzzz_y = cbuffer.data(gp_geom_11_off + 298 * ccomps * dcomps);

            auto g_z_x_xzzz_z = cbuffer.data(gp_geom_11_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzz_x, g_z_0_zzz_y, g_z_0_zzz_z, g_z_x_xzzz_x, g_z_x_xzzz_y, g_z_x_xzzz_z, g_z_x_zzz_x, g_z_x_zzz_xx, g_z_x_zzz_xy, g_z_x_zzz_xz, g_z_x_zzz_y, g_z_x_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xzzz_x[k] = g_z_0_zzz_x[k] - g_z_x_zzz_x[k] * ab_x + g_z_x_zzz_xx[k];

                g_z_x_xzzz_y[k] = g_z_0_zzz_y[k] - g_z_x_zzz_y[k] * ab_x + g_z_x_zzz_xy[k];

                g_z_x_xzzz_z[k] = g_z_0_zzz_z[k] - g_z_x_zzz_z[k] * ab_x + g_z_x_zzz_xz[k];
            }

            /// Set up 300-303 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyyy_x = cbuffer.data(gp_geom_11_off + 300 * ccomps * dcomps);

            auto g_z_x_yyyy_y = cbuffer.data(gp_geom_11_off + 301 * ccomps * dcomps);

            auto g_z_x_yyyy_z = cbuffer.data(gp_geom_11_off + 302 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yyy_x, g_z_x_yyy_xy, g_z_x_yyy_y, g_z_x_yyy_yy, g_z_x_yyy_yz, g_z_x_yyy_z, g_z_x_yyyy_x, g_z_x_yyyy_y, g_z_x_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyyy_x[k] = -g_z_x_yyy_x[k] * ab_y + g_z_x_yyy_xy[k];

                g_z_x_yyyy_y[k] = -g_z_x_yyy_y[k] * ab_y + g_z_x_yyy_yy[k];

                g_z_x_yyyy_z[k] = -g_z_x_yyy_z[k] * ab_y + g_z_x_yyy_yz[k];
            }

            /// Set up 303-306 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyyz_x = cbuffer.data(gp_geom_11_off + 303 * ccomps * dcomps);

            auto g_z_x_yyyz_y = cbuffer.data(gp_geom_11_off + 304 * ccomps * dcomps);

            auto g_z_x_yyyz_z = cbuffer.data(gp_geom_11_off + 305 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yyyz_x, g_z_x_yyyz_y, g_z_x_yyyz_z, g_z_x_yyz_x, g_z_x_yyz_xy, g_z_x_yyz_y, g_z_x_yyz_yy, g_z_x_yyz_yz, g_z_x_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyyz_x[k] = -g_z_x_yyz_x[k] * ab_y + g_z_x_yyz_xy[k];

                g_z_x_yyyz_y[k] = -g_z_x_yyz_y[k] * ab_y + g_z_x_yyz_yy[k];

                g_z_x_yyyz_z[k] = -g_z_x_yyz_z[k] * ab_y + g_z_x_yyz_yz[k];
            }

            /// Set up 306-309 components of targeted buffer : cbuffer.data(

            auto g_z_x_yyzz_x = cbuffer.data(gp_geom_11_off + 306 * ccomps * dcomps);

            auto g_z_x_yyzz_y = cbuffer.data(gp_geom_11_off + 307 * ccomps * dcomps);

            auto g_z_x_yyzz_z = cbuffer.data(gp_geom_11_off + 308 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yyzz_x, g_z_x_yyzz_y, g_z_x_yyzz_z, g_z_x_yzz_x, g_z_x_yzz_xy, g_z_x_yzz_y, g_z_x_yzz_yy, g_z_x_yzz_yz, g_z_x_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yyzz_x[k] = -g_z_x_yzz_x[k] * ab_y + g_z_x_yzz_xy[k];

                g_z_x_yyzz_y[k] = -g_z_x_yzz_y[k] * ab_y + g_z_x_yzz_yy[k];

                g_z_x_yyzz_z[k] = -g_z_x_yzz_z[k] * ab_y + g_z_x_yzz_yz[k];
            }

            /// Set up 309-312 components of targeted buffer : cbuffer.data(

            auto g_z_x_yzzz_x = cbuffer.data(gp_geom_11_off + 309 * ccomps * dcomps);

            auto g_z_x_yzzz_y = cbuffer.data(gp_geom_11_off + 310 * ccomps * dcomps);

            auto g_z_x_yzzz_z = cbuffer.data(gp_geom_11_off + 311 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yzzz_x, g_z_x_yzzz_y, g_z_x_yzzz_z, g_z_x_zzz_x, g_z_x_zzz_xy, g_z_x_zzz_y, g_z_x_zzz_yy, g_z_x_zzz_yz, g_z_x_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yzzz_x[k] = -g_z_x_zzz_x[k] * ab_y + g_z_x_zzz_xy[k];

                g_z_x_yzzz_y[k] = -g_z_x_zzz_y[k] * ab_y + g_z_x_zzz_yy[k];

                g_z_x_yzzz_z[k] = -g_z_x_zzz_z[k] * ab_y + g_z_x_zzz_yz[k];
            }

            /// Set up 312-315 components of targeted buffer : cbuffer.data(

            auto g_z_x_zzzz_x = cbuffer.data(gp_geom_11_off + 312 * ccomps * dcomps);

            auto g_z_x_zzzz_y = cbuffer.data(gp_geom_11_off + 313 * ccomps * dcomps);

            auto g_z_x_zzzz_z = cbuffer.data(gp_geom_11_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzz_x, g_0_x_zzz_y, g_0_x_zzz_z, g_z_x_zzz_x, g_z_x_zzz_xz, g_z_x_zzz_y, g_z_x_zzz_yz, g_z_x_zzz_z, g_z_x_zzz_zz, g_z_x_zzzz_x, g_z_x_zzzz_y, g_z_x_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_zzzz_x[k] = -g_0_x_zzz_x[k] - g_z_x_zzz_x[k] * ab_z + g_z_x_zzz_xz[k];

                g_z_x_zzzz_y[k] = -g_0_x_zzz_y[k] - g_z_x_zzz_y[k] * ab_z + g_z_x_zzz_yz[k];

                g_z_x_zzzz_z[k] = -g_0_x_zzz_z[k] - g_z_x_zzz_z[k] * ab_z + g_z_x_zzz_zz[k];
            }

            /// Set up 315-318 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxxx_x = cbuffer.data(gp_geom_11_off + 315 * ccomps * dcomps);

            auto g_z_y_xxxx_y = cbuffer.data(gp_geom_11_off + 316 * ccomps * dcomps);

            auto g_z_y_xxxx_z = cbuffer.data(gp_geom_11_off + 317 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxx_x, g_z_y_xxx_xx, g_z_y_xxx_xy, g_z_y_xxx_xz, g_z_y_xxx_y, g_z_y_xxx_z, g_z_y_xxxx_x, g_z_y_xxxx_y, g_z_y_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxxx_x[k] = -g_z_y_xxx_x[k] * ab_x + g_z_y_xxx_xx[k];

                g_z_y_xxxx_y[k] = -g_z_y_xxx_y[k] * ab_x + g_z_y_xxx_xy[k];

                g_z_y_xxxx_z[k] = -g_z_y_xxx_z[k] * ab_x + g_z_y_xxx_xz[k];
            }

            /// Set up 318-321 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxxy_x = cbuffer.data(gp_geom_11_off + 318 * ccomps * dcomps);

            auto g_z_y_xxxy_y = cbuffer.data(gp_geom_11_off + 319 * ccomps * dcomps);

            auto g_z_y_xxxy_z = cbuffer.data(gp_geom_11_off + 320 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxxy_x, g_z_y_xxxy_y, g_z_y_xxxy_z, g_z_y_xxy_x, g_z_y_xxy_xx, g_z_y_xxy_xy, g_z_y_xxy_xz, g_z_y_xxy_y, g_z_y_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxxy_x[k] = -g_z_y_xxy_x[k] * ab_x + g_z_y_xxy_xx[k];

                g_z_y_xxxy_y[k] = -g_z_y_xxy_y[k] * ab_x + g_z_y_xxy_xy[k];

                g_z_y_xxxy_z[k] = -g_z_y_xxy_z[k] * ab_x + g_z_y_xxy_xz[k];
            }

            /// Set up 321-324 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxxz_x = cbuffer.data(gp_geom_11_off + 321 * ccomps * dcomps);

            auto g_z_y_xxxz_y = cbuffer.data(gp_geom_11_off + 322 * ccomps * dcomps);

            auto g_z_y_xxxz_z = cbuffer.data(gp_geom_11_off + 323 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxxz_x, g_z_y_xxxz_y, g_z_y_xxxz_z, g_z_y_xxz_x, g_z_y_xxz_xx, g_z_y_xxz_xy, g_z_y_xxz_xz, g_z_y_xxz_y, g_z_y_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxxz_x[k] = -g_z_y_xxz_x[k] * ab_x + g_z_y_xxz_xx[k];

                g_z_y_xxxz_y[k] = -g_z_y_xxz_y[k] * ab_x + g_z_y_xxz_xy[k];

                g_z_y_xxxz_z[k] = -g_z_y_xxz_z[k] * ab_x + g_z_y_xxz_xz[k];
            }

            /// Set up 324-327 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxyy_x = cbuffer.data(gp_geom_11_off + 324 * ccomps * dcomps);

            auto g_z_y_xxyy_y = cbuffer.data(gp_geom_11_off + 325 * ccomps * dcomps);

            auto g_z_y_xxyy_z = cbuffer.data(gp_geom_11_off + 326 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxyy_x, g_z_y_xxyy_y, g_z_y_xxyy_z, g_z_y_xyy_x, g_z_y_xyy_xx, g_z_y_xyy_xy, g_z_y_xyy_xz, g_z_y_xyy_y, g_z_y_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxyy_x[k] = -g_z_y_xyy_x[k] * ab_x + g_z_y_xyy_xx[k];

                g_z_y_xxyy_y[k] = -g_z_y_xyy_y[k] * ab_x + g_z_y_xyy_xy[k];

                g_z_y_xxyy_z[k] = -g_z_y_xyy_z[k] * ab_x + g_z_y_xyy_xz[k];
            }

            /// Set up 327-330 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxyz_x = cbuffer.data(gp_geom_11_off + 327 * ccomps * dcomps);

            auto g_z_y_xxyz_y = cbuffer.data(gp_geom_11_off + 328 * ccomps * dcomps);

            auto g_z_y_xxyz_z = cbuffer.data(gp_geom_11_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxyz_x, g_z_y_xxyz_y, g_z_y_xxyz_z, g_z_y_xyz_x, g_z_y_xyz_xx, g_z_y_xyz_xy, g_z_y_xyz_xz, g_z_y_xyz_y, g_z_y_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxyz_x[k] = -g_z_y_xyz_x[k] * ab_x + g_z_y_xyz_xx[k];

                g_z_y_xxyz_y[k] = -g_z_y_xyz_y[k] * ab_x + g_z_y_xyz_xy[k];

                g_z_y_xxyz_z[k] = -g_z_y_xyz_z[k] * ab_x + g_z_y_xyz_xz[k];
            }

            /// Set up 330-333 components of targeted buffer : cbuffer.data(

            auto g_z_y_xxzz_x = cbuffer.data(gp_geom_11_off + 330 * ccomps * dcomps);

            auto g_z_y_xxzz_y = cbuffer.data(gp_geom_11_off + 331 * ccomps * dcomps);

            auto g_z_y_xxzz_z = cbuffer.data(gp_geom_11_off + 332 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xxzz_x, g_z_y_xxzz_y, g_z_y_xxzz_z, g_z_y_xzz_x, g_z_y_xzz_xx, g_z_y_xzz_xy, g_z_y_xzz_xz, g_z_y_xzz_y, g_z_y_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xxzz_x[k] = -g_z_y_xzz_x[k] * ab_x + g_z_y_xzz_xx[k];

                g_z_y_xxzz_y[k] = -g_z_y_xzz_y[k] * ab_x + g_z_y_xzz_xy[k];

                g_z_y_xxzz_z[k] = -g_z_y_xzz_z[k] * ab_x + g_z_y_xzz_xz[k];
            }

            /// Set up 333-336 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyyy_x = cbuffer.data(gp_geom_11_off + 333 * ccomps * dcomps);

            auto g_z_y_xyyy_y = cbuffer.data(gp_geom_11_off + 334 * ccomps * dcomps);

            auto g_z_y_xyyy_z = cbuffer.data(gp_geom_11_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyyy_x, g_z_y_xyyy_y, g_z_y_xyyy_z, g_z_y_yyy_x, g_z_y_yyy_xx, g_z_y_yyy_xy, g_z_y_yyy_xz, g_z_y_yyy_y, g_z_y_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyyy_x[k] = -g_z_y_yyy_x[k] * ab_x + g_z_y_yyy_xx[k];

                g_z_y_xyyy_y[k] = -g_z_y_yyy_y[k] * ab_x + g_z_y_yyy_xy[k];

                g_z_y_xyyy_z[k] = -g_z_y_yyy_z[k] * ab_x + g_z_y_yyy_xz[k];
            }

            /// Set up 336-339 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyyz_x = cbuffer.data(gp_geom_11_off + 336 * ccomps * dcomps);

            auto g_z_y_xyyz_y = cbuffer.data(gp_geom_11_off + 337 * ccomps * dcomps);

            auto g_z_y_xyyz_z = cbuffer.data(gp_geom_11_off + 338 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyyz_x, g_z_y_xyyz_y, g_z_y_xyyz_z, g_z_y_yyz_x, g_z_y_yyz_xx, g_z_y_yyz_xy, g_z_y_yyz_xz, g_z_y_yyz_y, g_z_y_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyyz_x[k] = -g_z_y_yyz_x[k] * ab_x + g_z_y_yyz_xx[k];

                g_z_y_xyyz_y[k] = -g_z_y_yyz_y[k] * ab_x + g_z_y_yyz_xy[k];

                g_z_y_xyyz_z[k] = -g_z_y_yyz_z[k] * ab_x + g_z_y_yyz_xz[k];
            }

            /// Set up 339-342 components of targeted buffer : cbuffer.data(

            auto g_z_y_xyzz_x = cbuffer.data(gp_geom_11_off + 339 * ccomps * dcomps);

            auto g_z_y_xyzz_y = cbuffer.data(gp_geom_11_off + 340 * ccomps * dcomps);

            auto g_z_y_xyzz_z = cbuffer.data(gp_geom_11_off + 341 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xyzz_x, g_z_y_xyzz_y, g_z_y_xyzz_z, g_z_y_yzz_x, g_z_y_yzz_xx, g_z_y_yzz_xy, g_z_y_yzz_xz, g_z_y_yzz_y, g_z_y_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xyzz_x[k] = -g_z_y_yzz_x[k] * ab_x + g_z_y_yzz_xx[k];

                g_z_y_xyzz_y[k] = -g_z_y_yzz_y[k] * ab_x + g_z_y_yzz_xy[k];

                g_z_y_xyzz_z[k] = -g_z_y_yzz_z[k] * ab_x + g_z_y_yzz_xz[k];
            }

            /// Set up 342-345 components of targeted buffer : cbuffer.data(

            auto g_z_y_xzzz_x = cbuffer.data(gp_geom_11_off + 342 * ccomps * dcomps);

            auto g_z_y_xzzz_y = cbuffer.data(gp_geom_11_off + 343 * ccomps * dcomps);

            auto g_z_y_xzzz_z = cbuffer.data(gp_geom_11_off + 344 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xzzz_x, g_z_y_xzzz_y, g_z_y_xzzz_z, g_z_y_zzz_x, g_z_y_zzz_xx, g_z_y_zzz_xy, g_z_y_zzz_xz, g_z_y_zzz_y, g_z_y_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xzzz_x[k] = -g_z_y_zzz_x[k] * ab_x + g_z_y_zzz_xx[k];

                g_z_y_xzzz_y[k] = -g_z_y_zzz_y[k] * ab_x + g_z_y_zzz_xy[k];

                g_z_y_xzzz_z[k] = -g_z_y_zzz_z[k] * ab_x + g_z_y_zzz_xz[k];
            }

            /// Set up 345-348 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyyy_x = cbuffer.data(gp_geom_11_off + 345 * ccomps * dcomps);

            auto g_z_y_yyyy_y = cbuffer.data(gp_geom_11_off + 346 * ccomps * dcomps);

            auto g_z_y_yyyy_z = cbuffer.data(gp_geom_11_off + 347 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyy_x, g_z_0_yyy_y, g_z_0_yyy_z, g_z_y_yyy_x, g_z_y_yyy_xy, g_z_y_yyy_y, g_z_y_yyy_yy, g_z_y_yyy_yz, g_z_y_yyy_z, g_z_y_yyyy_x, g_z_y_yyyy_y, g_z_y_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyyy_x[k] = g_z_0_yyy_x[k] - g_z_y_yyy_x[k] * ab_y + g_z_y_yyy_xy[k];

                g_z_y_yyyy_y[k] = g_z_0_yyy_y[k] - g_z_y_yyy_y[k] * ab_y + g_z_y_yyy_yy[k];

                g_z_y_yyyy_z[k] = g_z_0_yyy_z[k] - g_z_y_yyy_z[k] * ab_y + g_z_y_yyy_yz[k];
            }

            /// Set up 348-351 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyyz_x = cbuffer.data(gp_geom_11_off + 348 * ccomps * dcomps);

            auto g_z_y_yyyz_y = cbuffer.data(gp_geom_11_off + 349 * ccomps * dcomps);

            auto g_z_y_yyyz_z = cbuffer.data(gp_geom_11_off + 350 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyz_x, g_z_0_yyz_y, g_z_0_yyz_z, g_z_y_yyyz_x, g_z_y_yyyz_y, g_z_y_yyyz_z, g_z_y_yyz_x, g_z_y_yyz_xy, g_z_y_yyz_y, g_z_y_yyz_yy, g_z_y_yyz_yz, g_z_y_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyyz_x[k] = g_z_0_yyz_x[k] - g_z_y_yyz_x[k] * ab_y + g_z_y_yyz_xy[k];

                g_z_y_yyyz_y[k] = g_z_0_yyz_y[k] - g_z_y_yyz_y[k] * ab_y + g_z_y_yyz_yy[k];

                g_z_y_yyyz_z[k] = g_z_0_yyz_z[k] - g_z_y_yyz_z[k] * ab_y + g_z_y_yyz_yz[k];
            }

            /// Set up 351-354 components of targeted buffer : cbuffer.data(

            auto g_z_y_yyzz_x = cbuffer.data(gp_geom_11_off + 351 * ccomps * dcomps);

            auto g_z_y_yyzz_y = cbuffer.data(gp_geom_11_off + 352 * ccomps * dcomps);

            auto g_z_y_yyzz_z = cbuffer.data(gp_geom_11_off + 353 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yzz_x, g_z_0_yzz_y, g_z_0_yzz_z, g_z_y_yyzz_x, g_z_y_yyzz_y, g_z_y_yyzz_z, g_z_y_yzz_x, g_z_y_yzz_xy, g_z_y_yzz_y, g_z_y_yzz_yy, g_z_y_yzz_yz, g_z_y_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yyzz_x[k] = g_z_0_yzz_x[k] - g_z_y_yzz_x[k] * ab_y + g_z_y_yzz_xy[k];

                g_z_y_yyzz_y[k] = g_z_0_yzz_y[k] - g_z_y_yzz_y[k] * ab_y + g_z_y_yzz_yy[k];

                g_z_y_yyzz_z[k] = g_z_0_yzz_z[k] - g_z_y_yzz_z[k] * ab_y + g_z_y_yzz_yz[k];
            }

            /// Set up 354-357 components of targeted buffer : cbuffer.data(

            auto g_z_y_yzzz_x = cbuffer.data(gp_geom_11_off + 354 * ccomps * dcomps);

            auto g_z_y_yzzz_y = cbuffer.data(gp_geom_11_off + 355 * ccomps * dcomps);

            auto g_z_y_yzzz_z = cbuffer.data(gp_geom_11_off + 356 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzz_x, g_z_0_zzz_y, g_z_0_zzz_z, g_z_y_yzzz_x, g_z_y_yzzz_y, g_z_y_yzzz_z, g_z_y_zzz_x, g_z_y_zzz_xy, g_z_y_zzz_y, g_z_y_zzz_yy, g_z_y_zzz_yz, g_z_y_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yzzz_x[k] = g_z_0_zzz_x[k] - g_z_y_zzz_x[k] * ab_y + g_z_y_zzz_xy[k];

                g_z_y_yzzz_y[k] = g_z_0_zzz_y[k] - g_z_y_zzz_y[k] * ab_y + g_z_y_zzz_yy[k];

                g_z_y_yzzz_z[k] = g_z_0_zzz_z[k] - g_z_y_zzz_z[k] * ab_y + g_z_y_zzz_yz[k];
            }

            /// Set up 357-360 components of targeted buffer : cbuffer.data(

            auto g_z_y_zzzz_x = cbuffer.data(gp_geom_11_off + 357 * ccomps * dcomps);

            auto g_z_y_zzzz_y = cbuffer.data(gp_geom_11_off + 358 * ccomps * dcomps);

            auto g_z_y_zzzz_z = cbuffer.data(gp_geom_11_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzz_x, g_0_y_zzz_y, g_0_y_zzz_z, g_z_y_zzz_x, g_z_y_zzz_xz, g_z_y_zzz_y, g_z_y_zzz_yz, g_z_y_zzz_z, g_z_y_zzz_zz, g_z_y_zzzz_x, g_z_y_zzzz_y, g_z_y_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_zzzz_x[k] = -g_0_y_zzz_x[k] - g_z_y_zzz_x[k] * ab_z + g_z_y_zzz_xz[k];

                g_z_y_zzzz_y[k] = -g_0_y_zzz_y[k] - g_z_y_zzz_y[k] * ab_z + g_z_y_zzz_yz[k];

                g_z_y_zzzz_z[k] = -g_0_y_zzz_z[k] - g_z_y_zzz_z[k] * ab_z + g_z_y_zzz_zz[k];
            }

            /// Set up 360-363 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxxx_x = cbuffer.data(gp_geom_11_off + 360 * ccomps * dcomps);

            auto g_z_z_xxxx_y = cbuffer.data(gp_geom_11_off + 361 * ccomps * dcomps);

            auto g_z_z_xxxx_z = cbuffer.data(gp_geom_11_off + 362 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxx_x, g_z_z_xxx_xx, g_z_z_xxx_xy, g_z_z_xxx_xz, g_z_z_xxx_y, g_z_z_xxx_z, g_z_z_xxxx_x, g_z_z_xxxx_y, g_z_z_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxxx_x[k] = -g_z_z_xxx_x[k] * ab_x + g_z_z_xxx_xx[k];

                g_z_z_xxxx_y[k] = -g_z_z_xxx_y[k] * ab_x + g_z_z_xxx_xy[k];

                g_z_z_xxxx_z[k] = -g_z_z_xxx_z[k] * ab_x + g_z_z_xxx_xz[k];
            }

            /// Set up 363-366 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxxy_x = cbuffer.data(gp_geom_11_off + 363 * ccomps * dcomps);

            auto g_z_z_xxxy_y = cbuffer.data(gp_geom_11_off + 364 * ccomps * dcomps);

            auto g_z_z_xxxy_z = cbuffer.data(gp_geom_11_off + 365 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxxy_x, g_z_z_xxxy_y, g_z_z_xxxy_z, g_z_z_xxy_x, g_z_z_xxy_xx, g_z_z_xxy_xy, g_z_z_xxy_xz, g_z_z_xxy_y, g_z_z_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxxy_x[k] = -g_z_z_xxy_x[k] * ab_x + g_z_z_xxy_xx[k];

                g_z_z_xxxy_y[k] = -g_z_z_xxy_y[k] * ab_x + g_z_z_xxy_xy[k];

                g_z_z_xxxy_z[k] = -g_z_z_xxy_z[k] * ab_x + g_z_z_xxy_xz[k];
            }

            /// Set up 366-369 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxxz_x = cbuffer.data(gp_geom_11_off + 366 * ccomps * dcomps);

            auto g_z_z_xxxz_y = cbuffer.data(gp_geom_11_off + 367 * ccomps * dcomps);

            auto g_z_z_xxxz_z = cbuffer.data(gp_geom_11_off + 368 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxxz_x, g_z_z_xxxz_y, g_z_z_xxxz_z, g_z_z_xxz_x, g_z_z_xxz_xx, g_z_z_xxz_xy, g_z_z_xxz_xz, g_z_z_xxz_y, g_z_z_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxxz_x[k] = -g_z_z_xxz_x[k] * ab_x + g_z_z_xxz_xx[k];

                g_z_z_xxxz_y[k] = -g_z_z_xxz_y[k] * ab_x + g_z_z_xxz_xy[k];

                g_z_z_xxxz_z[k] = -g_z_z_xxz_z[k] * ab_x + g_z_z_xxz_xz[k];
            }

            /// Set up 369-372 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxyy_x = cbuffer.data(gp_geom_11_off + 369 * ccomps * dcomps);

            auto g_z_z_xxyy_y = cbuffer.data(gp_geom_11_off + 370 * ccomps * dcomps);

            auto g_z_z_xxyy_z = cbuffer.data(gp_geom_11_off + 371 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxyy_x, g_z_z_xxyy_y, g_z_z_xxyy_z, g_z_z_xyy_x, g_z_z_xyy_xx, g_z_z_xyy_xy, g_z_z_xyy_xz, g_z_z_xyy_y, g_z_z_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxyy_x[k] = -g_z_z_xyy_x[k] * ab_x + g_z_z_xyy_xx[k];

                g_z_z_xxyy_y[k] = -g_z_z_xyy_y[k] * ab_x + g_z_z_xyy_xy[k];

                g_z_z_xxyy_z[k] = -g_z_z_xyy_z[k] * ab_x + g_z_z_xyy_xz[k];
            }

            /// Set up 372-375 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxyz_x = cbuffer.data(gp_geom_11_off + 372 * ccomps * dcomps);

            auto g_z_z_xxyz_y = cbuffer.data(gp_geom_11_off + 373 * ccomps * dcomps);

            auto g_z_z_xxyz_z = cbuffer.data(gp_geom_11_off + 374 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxyz_x, g_z_z_xxyz_y, g_z_z_xxyz_z, g_z_z_xyz_x, g_z_z_xyz_xx, g_z_z_xyz_xy, g_z_z_xyz_xz, g_z_z_xyz_y, g_z_z_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxyz_x[k] = -g_z_z_xyz_x[k] * ab_x + g_z_z_xyz_xx[k];

                g_z_z_xxyz_y[k] = -g_z_z_xyz_y[k] * ab_x + g_z_z_xyz_xy[k];

                g_z_z_xxyz_z[k] = -g_z_z_xyz_z[k] * ab_x + g_z_z_xyz_xz[k];
            }

            /// Set up 375-378 components of targeted buffer : cbuffer.data(

            auto g_z_z_xxzz_x = cbuffer.data(gp_geom_11_off + 375 * ccomps * dcomps);

            auto g_z_z_xxzz_y = cbuffer.data(gp_geom_11_off + 376 * ccomps * dcomps);

            auto g_z_z_xxzz_z = cbuffer.data(gp_geom_11_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xxzz_x, g_z_z_xxzz_y, g_z_z_xxzz_z, g_z_z_xzz_x, g_z_z_xzz_xx, g_z_z_xzz_xy, g_z_z_xzz_xz, g_z_z_xzz_y, g_z_z_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xxzz_x[k] = -g_z_z_xzz_x[k] * ab_x + g_z_z_xzz_xx[k];

                g_z_z_xxzz_y[k] = -g_z_z_xzz_y[k] * ab_x + g_z_z_xzz_xy[k];

                g_z_z_xxzz_z[k] = -g_z_z_xzz_z[k] * ab_x + g_z_z_xzz_xz[k];
            }

            /// Set up 378-381 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyyy_x = cbuffer.data(gp_geom_11_off + 378 * ccomps * dcomps);

            auto g_z_z_xyyy_y = cbuffer.data(gp_geom_11_off + 379 * ccomps * dcomps);

            auto g_z_z_xyyy_z = cbuffer.data(gp_geom_11_off + 380 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyyy_x, g_z_z_xyyy_y, g_z_z_xyyy_z, g_z_z_yyy_x, g_z_z_yyy_xx, g_z_z_yyy_xy, g_z_z_yyy_xz, g_z_z_yyy_y, g_z_z_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyyy_x[k] = -g_z_z_yyy_x[k] * ab_x + g_z_z_yyy_xx[k];

                g_z_z_xyyy_y[k] = -g_z_z_yyy_y[k] * ab_x + g_z_z_yyy_xy[k];

                g_z_z_xyyy_z[k] = -g_z_z_yyy_z[k] * ab_x + g_z_z_yyy_xz[k];
            }

            /// Set up 381-384 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyyz_x = cbuffer.data(gp_geom_11_off + 381 * ccomps * dcomps);

            auto g_z_z_xyyz_y = cbuffer.data(gp_geom_11_off + 382 * ccomps * dcomps);

            auto g_z_z_xyyz_z = cbuffer.data(gp_geom_11_off + 383 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyyz_x, g_z_z_xyyz_y, g_z_z_xyyz_z, g_z_z_yyz_x, g_z_z_yyz_xx, g_z_z_yyz_xy, g_z_z_yyz_xz, g_z_z_yyz_y, g_z_z_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyyz_x[k] = -g_z_z_yyz_x[k] * ab_x + g_z_z_yyz_xx[k];

                g_z_z_xyyz_y[k] = -g_z_z_yyz_y[k] * ab_x + g_z_z_yyz_xy[k];

                g_z_z_xyyz_z[k] = -g_z_z_yyz_z[k] * ab_x + g_z_z_yyz_xz[k];
            }

            /// Set up 384-387 components of targeted buffer : cbuffer.data(

            auto g_z_z_xyzz_x = cbuffer.data(gp_geom_11_off + 384 * ccomps * dcomps);

            auto g_z_z_xyzz_y = cbuffer.data(gp_geom_11_off + 385 * ccomps * dcomps);

            auto g_z_z_xyzz_z = cbuffer.data(gp_geom_11_off + 386 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xyzz_x, g_z_z_xyzz_y, g_z_z_xyzz_z, g_z_z_yzz_x, g_z_z_yzz_xx, g_z_z_yzz_xy, g_z_z_yzz_xz, g_z_z_yzz_y, g_z_z_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xyzz_x[k] = -g_z_z_yzz_x[k] * ab_x + g_z_z_yzz_xx[k];

                g_z_z_xyzz_y[k] = -g_z_z_yzz_y[k] * ab_x + g_z_z_yzz_xy[k];

                g_z_z_xyzz_z[k] = -g_z_z_yzz_z[k] * ab_x + g_z_z_yzz_xz[k];
            }

            /// Set up 387-390 components of targeted buffer : cbuffer.data(

            auto g_z_z_xzzz_x = cbuffer.data(gp_geom_11_off + 387 * ccomps * dcomps);

            auto g_z_z_xzzz_y = cbuffer.data(gp_geom_11_off + 388 * ccomps * dcomps);

            auto g_z_z_xzzz_z = cbuffer.data(gp_geom_11_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xzzz_x, g_z_z_xzzz_y, g_z_z_xzzz_z, g_z_z_zzz_x, g_z_z_zzz_xx, g_z_z_zzz_xy, g_z_z_zzz_xz, g_z_z_zzz_y, g_z_z_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xzzz_x[k] = -g_z_z_zzz_x[k] * ab_x + g_z_z_zzz_xx[k];

                g_z_z_xzzz_y[k] = -g_z_z_zzz_y[k] * ab_x + g_z_z_zzz_xy[k];

                g_z_z_xzzz_z[k] = -g_z_z_zzz_z[k] * ab_x + g_z_z_zzz_xz[k];
            }

            /// Set up 390-393 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyyy_x = cbuffer.data(gp_geom_11_off + 390 * ccomps * dcomps);

            auto g_z_z_yyyy_y = cbuffer.data(gp_geom_11_off + 391 * ccomps * dcomps);

            auto g_z_z_yyyy_z = cbuffer.data(gp_geom_11_off + 392 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yyy_x, g_z_z_yyy_xy, g_z_z_yyy_y, g_z_z_yyy_yy, g_z_z_yyy_yz, g_z_z_yyy_z, g_z_z_yyyy_x, g_z_z_yyyy_y, g_z_z_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyyy_x[k] = -g_z_z_yyy_x[k] * ab_y + g_z_z_yyy_xy[k];

                g_z_z_yyyy_y[k] = -g_z_z_yyy_y[k] * ab_y + g_z_z_yyy_yy[k];

                g_z_z_yyyy_z[k] = -g_z_z_yyy_z[k] * ab_y + g_z_z_yyy_yz[k];
            }

            /// Set up 393-396 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyyz_x = cbuffer.data(gp_geom_11_off + 393 * ccomps * dcomps);

            auto g_z_z_yyyz_y = cbuffer.data(gp_geom_11_off + 394 * ccomps * dcomps);

            auto g_z_z_yyyz_z = cbuffer.data(gp_geom_11_off + 395 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yyyz_x, g_z_z_yyyz_y, g_z_z_yyyz_z, g_z_z_yyz_x, g_z_z_yyz_xy, g_z_z_yyz_y, g_z_z_yyz_yy, g_z_z_yyz_yz, g_z_z_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyyz_x[k] = -g_z_z_yyz_x[k] * ab_y + g_z_z_yyz_xy[k];

                g_z_z_yyyz_y[k] = -g_z_z_yyz_y[k] * ab_y + g_z_z_yyz_yy[k];

                g_z_z_yyyz_z[k] = -g_z_z_yyz_z[k] * ab_y + g_z_z_yyz_yz[k];
            }

            /// Set up 396-399 components of targeted buffer : cbuffer.data(

            auto g_z_z_yyzz_x = cbuffer.data(gp_geom_11_off + 396 * ccomps * dcomps);

            auto g_z_z_yyzz_y = cbuffer.data(gp_geom_11_off + 397 * ccomps * dcomps);

            auto g_z_z_yyzz_z = cbuffer.data(gp_geom_11_off + 398 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yyzz_x, g_z_z_yyzz_y, g_z_z_yyzz_z, g_z_z_yzz_x, g_z_z_yzz_xy, g_z_z_yzz_y, g_z_z_yzz_yy, g_z_z_yzz_yz, g_z_z_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yyzz_x[k] = -g_z_z_yzz_x[k] * ab_y + g_z_z_yzz_xy[k];

                g_z_z_yyzz_y[k] = -g_z_z_yzz_y[k] * ab_y + g_z_z_yzz_yy[k];

                g_z_z_yyzz_z[k] = -g_z_z_yzz_z[k] * ab_y + g_z_z_yzz_yz[k];
            }

            /// Set up 399-402 components of targeted buffer : cbuffer.data(

            auto g_z_z_yzzz_x = cbuffer.data(gp_geom_11_off + 399 * ccomps * dcomps);

            auto g_z_z_yzzz_y = cbuffer.data(gp_geom_11_off + 400 * ccomps * dcomps);

            auto g_z_z_yzzz_z = cbuffer.data(gp_geom_11_off + 401 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yzzz_x, g_z_z_yzzz_y, g_z_z_yzzz_z, g_z_z_zzz_x, g_z_z_zzz_xy, g_z_z_zzz_y, g_z_z_zzz_yy, g_z_z_zzz_yz, g_z_z_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yzzz_x[k] = -g_z_z_zzz_x[k] * ab_y + g_z_z_zzz_xy[k];

                g_z_z_yzzz_y[k] = -g_z_z_zzz_y[k] * ab_y + g_z_z_zzz_yy[k];

                g_z_z_yzzz_z[k] = -g_z_z_zzz_z[k] * ab_y + g_z_z_zzz_yz[k];
            }

            /// Set up 402-405 components of targeted buffer : cbuffer.data(

            auto g_z_z_zzzz_x = cbuffer.data(gp_geom_11_off + 402 * ccomps * dcomps);

            auto g_z_z_zzzz_y = cbuffer.data(gp_geom_11_off + 403 * ccomps * dcomps);

            auto g_z_z_zzzz_z = cbuffer.data(gp_geom_11_off + 404 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzz_x, g_0_z_zzz_y, g_0_z_zzz_z, g_z_0_zzz_x, g_z_0_zzz_y, g_z_0_zzz_z, g_z_z_zzz_x, g_z_z_zzz_xz, g_z_z_zzz_y, g_z_z_zzz_yz, g_z_z_zzz_z, g_z_z_zzz_zz, g_z_z_zzzz_x, g_z_z_zzzz_y, g_z_z_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_zzzz_x[k] = -g_0_z_zzz_x[k] + g_z_0_zzz_x[k] - g_z_z_zzz_x[k] * ab_z + g_z_z_zzz_xz[k];

                g_z_z_zzzz_y[k] = -g_0_z_zzz_y[k] + g_z_0_zzz_y[k] - g_z_z_zzz_y[k] * ab_z + g_z_z_zzz_yz[k];

                g_z_z_zzzz_z[k] = -g_0_z_zzz_z[k] + g_z_0_zzz_z[k] - g_z_z_zzz_z[k] * ab_z + g_z_z_zzz_zz[k];
            }
        }
    }
}

} // erirec namespace

