#include "ElectronRepulsionGeom1010ContrRecGSXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom1010_hrr_electron_repulsion_gsxx(CSimdArray<double>& cbuffer,
                                              const size_t idx_geom_1010_gsxx,
                                              const size_t idx_geom_0010_fsxx,
                                              const size_t idx_geom_1010_fsxx,
                                              const size_t idx_geom_1010_fpxx,
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
            /// Set up components of auxilary buffer : FSSS

            const auto fs_geom_0010_off = idx_geom_0010_fsxx + i * dcomps + j;

            auto g_0_0_x_0_xxx_0 = cbuffer.data(fs_geom_0010_off + 0 * ccomps * dcomps);

            auto g_0_0_x_0_xxy_0 = cbuffer.data(fs_geom_0010_off + 1 * ccomps * dcomps);

            auto g_0_0_x_0_xxz_0 = cbuffer.data(fs_geom_0010_off + 2 * ccomps * dcomps);

            auto g_0_0_x_0_xyy_0 = cbuffer.data(fs_geom_0010_off + 3 * ccomps * dcomps);

            auto g_0_0_x_0_xyz_0 = cbuffer.data(fs_geom_0010_off + 4 * ccomps * dcomps);

            auto g_0_0_x_0_xzz_0 = cbuffer.data(fs_geom_0010_off + 5 * ccomps * dcomps);

            auto g_0_0_x_0_yyy_0 = cbuffer.data(fs_geom_0010_off + 6 * ccomps * dcomps);

            auto g_0_0_x_0_yyz_0 = cbuffer.data(fs_geom_0010_off + 7 * ccomps * dcomps);

            auto g_0_0_x_0_yzz_0 = cbuffer.data(fs_geom_0010_off + 8 * ccomps * dcomps);

            auto g_0_0_x_0_zzz_0 = cbuffer.data(fs_geom_0010_off + 9 * ccomps * dcomps);

            auto g_0_0_y_0_xxx_0 = cbuffer.data(fs_geom_0010_off + 10 * ccomps * dcomps);

            auto g_0_0_y_0_xxy_0 = cbuffer.data(fs_geom_0010_off + 11 * ccomps * dcomps);

            auto g_0_0_y_0_xxz_0 = cbuffer.data(fs_geom_0010_off + 12 * ccomps * dcomps);

            auto g_0_0_y_0_xyy_0 = cbuffer.data(fs_geom_0010_off + 13 * ccomps * dcomps);

            auto g_0_0_y_0_xyz_0 = cbuffer.data(fs_geom_0010_off + 14 * ccomps * dcomps);

            auto g_0_0_y_0_xzz_0 = cbuffer.data(fs_geom_0010_off + 15 * ccomps * dcomps);

            auto g_0_0_y_0_yyy_0 = cbuffer.data(fs_geom_0010_off + 16 * ccomps * dcomps);

            auto g_0_0_y_0_yyz_0 = cbuffer.data(fs_geom_0010_off + 17 * ccomps * dcomps);

            auto g_0_0_y_0_yzz_0 = cbuffer.data(fs_geom_0010_off + 18 * ccomps * dcomps);

            auto g_0_0_y_0_zzz_0 = cbuffer.data(fs_geom_0010_off + 19 * ccomps * dcomps);

            auto g_0_0_z_0_xxx_0 = cbuffer.data(fs_geom_0010_off + 20 * ccomps * dcomps);

            auto g_0_0_z_0_xxy_0 = cbuffer.data(fs_geom_0010_off + 21 * ccomps * dcomps);

            auto g_0_0_z_0_xxz_0 = cbuffer.data(fs_geom_0010_off + 22 * ccomps * dcomps);

            auto g_0_0_z_0_xyy_0 = cbuffer.data(fs_geom_0010_off + 23 * ccomps * dcomps);

            auto g_0_0_z_0_xyz_0 = cbuffer.data(fs_geom_0010_off + 24 * ccomps * dcomps);

            auto g_0_0_z_0_xzz_0 = cbuffer.data(fs_geom_0010_off + 25 * ccomps * dcomps);

            auto g_0_0_z_0_yyy_0 = cbuffer.data(fs_geom_0010_off + 26 * ccomps * dcomps);

            auto g_0_0_z_0_yyz_0 = cbuffer.data(fs_geom_0010_off + 27 * ccomps * dcomps);

            auto g_0_0_z_0_yzz_0 = cbuffer.data(fs_geom_0010_off + 28 * ccomps * dcomps);

            auto g_0_0_z_0_zzz_0 = cbuffer.data(fs_geom_0010_off + 29 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FSSS

            const auto fs_geom_1010_off = idx_geom_1010_fsxx + i * dcomps + j;

            auto g_x_0_x_0_xxx_0 = cbuffer.data(fs_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_0 = cbuffer.data(fs_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_0 = cbuffer.data(fs_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_0 = cbuffer.data(fs_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_0 = cbuffer.data(fs_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_0 = cbuffer.data(fs_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_0 = cbuffer.data(fs_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_0 = cbuffer.data(fs_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_0 = cbuffer.data(fs_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_0 = cbuffer.data(fs_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_0 = cbuffer.data(fs_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_0 = cbuffer.data(fs_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_0 = cbuffer.data(fs_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_0 = cbuffer.data(fs_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_0 = cbuffer.data(fs_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_0 = cbuffer.data(fs_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_0 = cbuffer.data(fs_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_0 = cbuffer.data(fs_geom_1010_off + 17 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_0 = cbuffer.data(fs_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_0 = cbuffer.data(fs_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_0 = cbuffer.data(fs_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_0 = cbuffer.data(fs_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_0 = cbuffer.data(fs_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_0 = cbuffer.data(fs_geom_1010_off + 23 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_0 = cbuffer.data(fs_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_0 = cbuffer.data(fs_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_0 = cbuffer.data(fs_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_0 = cbuffer.data(fs_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_0 = cbuffer.data(fs_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_0 = cbuffer.data(fs_geom_1010_off + 29 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_0 = cbuffer.data(fs_geom_1010_off + 30 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_0 = cbuffer.data(fs_geom_1010_off + 31 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_0 = cbuffer.data(fs_geom_1010_off + 32 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_0 = cbuffer.data(fs_geom_1010_off + 33 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_0 = cbuffer.data(fs_geom_1010_off + 34 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_0 = cbuffer.data(fs_geom_1010_off + 35 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_0 = cbuffer.data(fs_geom_1010_off + 36 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_0 = cbuffer.data(fs_geom_1010_off + 37 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_0 = cbuffer.data(fs_geom_1010_off + 38 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_0 = cbuffer.data(fs_geom_1010_off + 39 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_0 = cbuffer.data(fs_geom_1010_off + 40 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_0 = cbuffer.data(fs_geom_1010_off + 41 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_0 = cbuffer.data(fs_geom_1010_off + 42 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_0 = cbuffer.data(fs_geom_1010_off + 43 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_0 = cbuffer.data(fs_geom_1010_off + 44 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_0 = cbuffer.data(fs_geom_1010_off + 45 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_0 = cbuffer.data(fs_geom_1010_off + 46 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_0 = cbuffer.data(fs_geom_1010_off + 47 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_0 = cbuffer.data(fs_geom_1010_off + 48 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_0 = cbuffer.data(fs_geom_1010_off + 49 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_0 = cbuffer.data(fs_geom_1010_off + 50 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_0 = cbuffer.data(fs_geom_1010_off + 51 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_0 = cbuffer.data(fs_geom_1010_off + 52 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_0 = cbuffer.data(fs_geom_1010_off + 53 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_0 = cbuffer.data(fs_geom_1010_off + 54 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_0 = cbuffer.data(fs_geom_1010_off + 55 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_0 = cbuffer.data(fs_geom_1010_off + 56 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_0 = cbuffer.data(fs_geom_1010_off + 57 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_0 = cbuffer.data(fs_geom_1010_off + 58 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_0 = cbuffer.data(fs_geom_1010_off + 59 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_0 = cbuffer.data(fs_geom_1010_off + 60 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_0 = cbuffer.data(fs_geom_1010_off + 61 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_0 = cbuffer.data(fs_geom_1010_off + 62 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_0 = cbuffer.data(fs_geom_1010_off + 63 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_0 = cbuffer.data(fs_geom_1010_off + 64 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_0 = cbuffer.data(fs_geom_1010_off + 65 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_0 = cbuffer.data(fs_geom_1010_off + 66 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_0 = cbuffer.data(fs_geom_1010_off + 67 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_0 = cbuffer.data(fs_geom_1010_off + 68 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_0 = cbuffer.data(fs_geom_1010_off + 69 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_0 = cbuffer.data(fs_geom_1010_off + 70 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_0 = cbuffer.data(fs_geom_1010_off + 71 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_0 = cbuffer.data(fs_geom_1010_off + 72 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_0 = cbuffer.data(fs_geom_1010_off + 73 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_0 = cbuffer.data(fs_geom_1010_off + 74 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_0 = cbuffer.data(fs_geom_1010_off + 75 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_0 = cbuffer.data(fs_geom_1010_off + 76 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_0 = cbuffer.data(fs_geom_1010_off + 77 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_0 = cbuffer.data(fs_geom_1010_off + 78 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_0 = cbuffer.data(fs_geom_1010_off + 79 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_0 = cbuffer.data(fs_geom_1010_off + 80 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_0 = cbuffer.data(fs_geom_1010_off + 81 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_0 = cbuffer.data(fs_geom_1010_off + 82 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_0 = cbuffer.data(fs_geom_1010_off + 83 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_0 = cbuffer.data(fs_geom_1010_off + 84 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_0 = cbuffer.data(fs_geom_1010_off + 85 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_0 = cbuffer.data(fs_geom_1010_off + 86 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_0 = cbuffer.data(fs_geom_1010_off + 87 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_0 = cbuffer.data(fs_geom_1010_off + 88 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_0 = cbuffer.data(fs_geom_1010_off + 89 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FPSS

            const auto fp_geom_1010_off = idx_geom_1010_fpxx + i * dcomps + j;

            auto g_x_0_x_0_xxx_x = cbuffer.data(fp_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_y = cbuffer.data(fp_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_z = cbuffer.data(fp_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_x = cbuffer.data(fp_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_y = cbuffer.data(fp_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_z = cbuffer.data(fp_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_x = cbuffer.data(fp_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_y = cbuffer.data(fp_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_z = cbuffer.data(fp_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_x = cbuffer.data(fp_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_y = cbuffer.data(fp_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_z = cbuffer.data(fp_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_x = cbuffer.data(fp_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_y = cbuffer.data(fp_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_z = cbuffer.data(fp_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_x = cbuffer.data(fp_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_y = cbuffer.data(fp_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_z = cbuffer.data(fp_geom_1010_off + 17 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_x = cbuffer.data(fp_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_y = cbuffer.data(fp_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_z = cbuffer.data(fp_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_x = cbuffer.data(fp_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_y = cbuffer.data(fp_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_z = cbuffer.data(fp_geom_1010_off + 23 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_x = cbuffer.data(fp_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_y = cbuffer.data(fp_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_z = cbuffer.data(fp_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_x = cbuffer.data(fp_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_y = cbuffer.data(fp_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_z = cbuffer.data(fp_geom_1010_off + 29 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_x = cbuffer.data(fp_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_y = cbuffer.data(fp_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_z = cbuffer.data(fp_geom_1010_off + 32 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_x = cbuffer.data(fp_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_y = cbuffer.data(fp_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_z = cbuffer.data(fp_geom_1010_off + 35 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_x = cbuffer.data(fp_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_y = cbuffer.data(fp_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_z = cbuffer.data(fp_geom_1010_off + 38 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_x = cbuffer.data(fp_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_y = cbuffer.data(fp_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_z = cbuffer.data(fp_geom_1010_off + 41 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_x = cbuffer.data(fp_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_y = cbuffer.data(fp_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_z = cbuffer.data(fp_geom_1010_off + 44 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_x = cbuffer.data(fp_geom_1010_off + 45 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_y = cbuffer.data(fp_geom_1010_off + 46 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_z = cbuffer.data(fp_geom_1010_off + 47 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_x = cbuffer.data(fp_geom_1010_off + 48 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_y = cbuffer.data(fp_geom_1010_off + 49 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_z = cbuffer.data(fp_geom_1010_off + 50 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_x = cbuffer.data(fp_geom_1010_off + 51 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_y = cbuffer.data(fp_geom_1010_off + 52 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_z = cbuffer.data(fp_geom_1010_off + 53 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_x = cbuffer.data(fp_geom_1010_off + 54 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_y = cbuffer.data(fp_geom_1010_off + 55 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_z = cbuffer.data(fp_geom_1010_off + 56 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_x = cbuffer.data(fp_geom_1010_off + 57 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_y = cbuffer.data(fp_geom_1010_off + 58 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_z = cbuffer.data(fp_geom_1010_off + 59 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_x = cbuffer.data(fp_geom_1010_off + 60 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_y = cbuffer.data(fp_geom_1010_off + 61 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_z = cbuffer.data(fp_geom_1010_off + 62 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_x = cbuffer.data(fp_geom_1010_off + 63 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_y = cbuffer.data(fp_geom_1010_off + 64 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_z = cbuffer.data(fp_geom_1010_off + 65 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_x = cbuffer.data(fp_geom_1010_off + 66 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_y = cbuffer.data(fp_geom_1010_off + 67 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_z = cbuffer.data(fp_geom_1010_off + 68 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_x = cbuffer.data(fp_geom_1010_off + 69 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_y = cbuffer.data(fp_geom_1010_off + 70 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_z = cbuffer.data(fp_geom_1010_off + 71 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_x = cbuffer.data(fp_geom_1010_off + 72 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_y = cbuffer.data(fp_geom_1010_off + 73 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_z = cbuffer.data(fp_geom_1010_off + 74 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_x = cbuffer.data(fp_geom_1010_off + 75 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_y = cbuffer.data(fp_geom_1010_off + 76 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_z = cbuffer.data(fp_geom_1010_off + 77 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_x = cbuffer.data(fp_geom_1010_off + 78 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_y = cbuffer.data(fp_geom_1010_off + 79 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_z = cbuffer.data(fp_geom_1010_off + 80 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_x = cbuffer.data(fp_geom_1010_off + 81 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_y = cbuffer.data(fp_geom_1010_off + 82 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_z = cbuffer.data(fp_geom_1010_off + 83 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_x = cbuffer.data(fp_geom_1010_off + 84 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_y = cbuffer.data(fp_geom_1010_off + 85 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_z = cbuffer.data(fp_geom_1010_off + 86 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_x = cbuffer.data(fp_geom_1010_off + 87 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_y = cbuffer.data(fp_geom_1010_off + 88 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_z = cbuffer.data(fp_geom_1010_off + 89 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_x = cbuffer.data(fp_geom_1010_off + 90 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_y = cbuffer.data(fp_geom_1010_off + 91 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_z = cbuffer.data(fp_geom_1010_off + 92 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_x = cbuffer.data(fp_geom_1010_off + 93 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_y = cbuffer.data(fp_geom_1010_off + 94 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_z = cbuffer.data(fp_geom_1010_off + 95 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_x = cbuffer.data(fp_geom_1010_off + 96 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_y = cbuffer.data(fp_geom_1010_off + 97 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_z = cbuffer.data(fp_geom_1010_off + 98 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_x = cbuffer.data(fp_geom_1010_off + 99 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_y = cbuffer.data(fp_geom_1010_off + 100 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_z = cbuffer.data(fp_geom_1010_off + 101 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_x = cbuffer.data(fp_geom_1010_off + 102 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_y = cbuffer.data(fp_geom_1010_off + 103 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_z = cbuffer.data(fp_geom_1010_off + 104 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_x = cbuffer.data(fp_geom_1010_off + 105 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_y = cbuffer.data(fp_geom_1010_off + 106 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_z = cbuffer.data(fp_geom_1010_off + 107 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_x = cbuffer.data(fp_geom_1010_off + 108 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_y = cbuffer.data(fp_geom_1010_off + 109 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_z = cbuffer.data(fp_geom_1010_off + 110 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_x = cbuffer.data(fp_geom_1010_off + 111 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_y = cbuffer.data(fp_geom_1010_off + 112 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_z = cbuffer.data(fp_geom_1010_off + 113 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_x = cbuffer.data(fp_geom_1010_off + 114 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_y = cbuffer.data(fp_geom_1010_off + 115 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_z = cbuffer.data(fp_geom_1010_off + 116 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_x = cbuffer.data(fp_geom_1010_off + 117 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_y = cbuffer.data(fp_geom_1010_off + 118 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_z = cbuffer.data(fp_geom_1010_off + 119 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_x = cbuffer.data(fp_geom_1010_off + 120 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_y = cbuffer.data(fp_geom_1010_off + 121 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_z = cbuffer.data(fp_geom_1010_off + 122 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_x = cbuffer.data(fp_geom_1010_off + 123 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_y = cbuffer.data(fp_geom_1010_off + 124 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_z = cbuffer.data(fp_geom_1010_off + 125 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_x = cbuffer.data(fp_geom_1010_off + 126 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_y = cbuffer.data(fp_geom_1010_off + 127 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_z = cbuffer.data(fp_geom_1010_off + 128 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_x = cbuffer.data(fp_geom_1010_off + 129 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_y = cbuffer.data(fp_geom_1010_off + 130 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_z = cbuffer.data(fp_geom_1010_off + 131 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_x = cbuffer.data(fp_geom_1010_off + 132 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_y = cbuffer.data(fp_geom_1010_off + 133 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_z = cbuffer.data(fp_geom_1010_off + 134 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_x = cbuffer.data(fp_geom_1010_off + 135 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_y = cbuffer.data(fp_geom_1010_off + 136 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_z = cbuffer.data(fp_geom_1010_off + 137 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_x = cbuffer.data(fp_geom_1010_off + 138 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_y = cbuffer.data(fp_geom_1010_off + 139 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_z = cbuffer.data(fp_geom_1010_off + 140 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_x = cbuffer.data(fp_geom_1010_off + 141 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_y = cbuffer.data(fp_geom_1010_off + 142 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_z = cbuffer.data(fp_geom_1010_off + 143 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_x = cbuffer.data(fp_geom_1010_off + 144 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_y = cbuffer.data(fp_geom_1010_off + 145 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_z = cbuffer.data(fp_geom_1010_off + 146 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_x = cbuffer.data(fp_geom_1010_off + 147 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_y = cbuffer.data(fp_geom_1010_off + 148 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_z = cbuffer.data(fp_geom_1010_off + 149 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_x = cbuffer.data(fp_geom_1010_off + 150 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_y = cbuffer.data(fp_geom_1010_off + 151 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_z = cbuffer.data(fp_geom_1010_off + 152 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_x = cbuffer.data(fp_geom_1010_off + 153 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_y = cbuffer.data(fp_geom_1010_off + 154 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_z = cbuffer.data(fp_geom_1010_off + 155 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_x = cbuffer.data(fp_geom_1010_off + 156 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_y = cbuffer.data(fp_geom_1010_off + 157 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_z = cbuffer.data(fp_geom_1010_off + 158 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_x = cbuffer.data(fp_geom_1010_off + 159 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_y = cbuffer.data(fp_geom_1010_off + 160 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_z = cbuffer.data(fp_geom_1010_off + 161 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_x = cbuffer.data(fp_geom_1010_off + 162 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_y = cbuffer.data(fp_geom_1010_off + 163 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_z = cbuffer.data(fp_geom_1010_off + 164 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_x = cbuffer.data(fp_geom_1010_off + 165 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_y = cbuffer.data(fp_geom_1010_off + 166 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_z = cbuffer.data(fp_geom_1010_off + 167 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_x = cbuffer.data(fp_geom_1010_off + 168 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_y = cbuffer.data(fp_geom_1010_off + 169 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_z = cbuffer.data(fp_geom_1010_off + 170 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_x = cbuffer.data(fp_geom_1010_off + 171 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_y = cbuffer.data(fp_geom_1010_off + 172 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_z = cbuffer.data(fp_geom_1010_off + 173 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_x = cbuffer.data(fp_geom_1010_off + 174 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_y = cbuffer.data(fp_geom_1010_off + 175 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_z = cbuffer.data(fp_geom_1010_off + 176 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_x = cbuffer.data(fp_geom_1010_off + 177 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_y = cbuffer.data(fp_geom_1010_off + 178 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_z = cbuffer.data(fp_geom_1010_off + 179 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_x = cbuffer.data(fp_geom_1010_off + 180 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_y = cbuffer.data(fp_geom_1010_off + 181 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_z = cbuffer.data(fp_geom_1010_off + 182 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_x = cbuffer.data(fp_geom_1010_off + 183 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_y = cbuffer.data(fp_geom_1010_off + 184 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_z = cbuffer.data(fp_geom_1010_off + 185 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_x = cbuffer.data(fp_geom_1010_off + 186 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_y = cbuffer.data(fp_geom_1010_off + 187 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_z = cbuffer.data(fp_geom_1010_off + 188 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_x = cbuffer.data(fp_geom_1010_off + 189 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_y = cbuffer.data(fp_geom_1010_off + 190 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_z = cbuffer.data(fp_geom_1010_off + 191 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_x = cbuffer.data(fp_geom_1010_off + 192 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_y = cbuffer.data(fp_geom_1010_off + 193 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_z = cbuffer.data(fp_geom_1010_off + 194 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_x = cbuffer.data(fp_geom_1010_off + 195 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_y = cbuffer.data(fp_geom_1010_off + 196 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_z = cbuffer.data(fp_geom_1010_off + 197 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_x = cbuffer.data(fp_geom_1010_off + 198 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_y = cbuffer.data(fp_geom_1010_off + 199 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_z = cbuffer.data(fp_geom_1010_off + 200 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_x = cbuffer.data(fp_geom_1010_off + 201 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_y = cbuffer.data(fp_geom_1010_off + 202 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_z = cbuffer.data(fp_geom_1010_off + 203 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_x = cbuffer.data(fp_geom_1010_off + 204 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_y = cbuffer.data(fp_geom_1010_off + 205 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_z = cbuffer.data(fp_geom_1010_off + 206 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_x = cbuffer.data(fp_geom_1010_off + 207 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_y = cbuffer.data(fp_geom_1010_off + 208 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_z = cbuffer.data(fp_geom_1010_off + 209 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_x = cbuffer.data(fp_geom_1010_off + 210 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_y = cbuffer.data(fp_geom_1010_off + 211 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_z = cbuffer.data(fp_geom_1010_off + 212 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_x = cbuffer.data(fp_geom_1010_off + 213 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_y = cbuffer.data(fp_geom_1010_off + 214 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_z = cbuffer.data(fp_geom_1010_off + 215 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_x = cbuffer.data(fp_geom_1010_off + 216 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_y = cbuffer.data(fp_geom_1010_off + 217 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_z = cbuffer.data(fp_geom_1010_off + 218 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_x = cbuffer.data(fp_geom_1010_off + 219 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_y = cbuffer.data(fp_geom_1010_off + 220 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_z = cbuffer.data(fp_geom_1010_off + 221 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_x = cbuffer.data(fp_geom_1010_off + 222 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_y = cbuffer.data(fp_geom_1010_off + 223 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_z = cbuffer.data(fp_geom_1010_off + 224 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_x = cbuffer.data(fp_geom_1010_off + 225 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_y = cbuffer.data(fp_geom_1010_off + 226 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_z = cbuffer.data(fp_geom_1010_off + 227 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_x = cbuffer.data(fp_geom_1010_off + 228 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_y = cbuffer.data(fp_geom_1010_off + 229 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_z = cbuffer.data(fp_geom_1010_off + 230 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_x = cbuffer.data(fp_geom_1010_off + 231 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_y = cbuffer.data(fp_geom_1010_off + 232 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_z = cbuffer.data(fp_geom_1010_off + 233 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_x = cbuffer.data(fp_geom_1010_off + 234 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_y = cbuffer.data(fp_geom_1010_off + 235 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_z = cbuffer.data(fp_geom_1010_off + 236 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_x = cbuffer.data(fp_geom_1010_off + 237 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_y = cbuffer.data(fp_geom_1010_off + 238 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_z = cbuffer.data(fp_geom_1010_off + 239 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_x = cbuffer.data(fp_geom_1010_off + 240 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_y = cbuffer.data(fp_geom_1010_off + 241 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_z = cbuffer.data(fp_geom_1010_off + 242 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_x = cbuffer.data(fp_geom_1010_off + 243 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_y = cbuffer.data(fp_geom_1010_off + 244 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_z = cbuffer.data(fp_geom_1010_off + 245 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_x = cbuffer.data(fp_geom_1010_off + 246 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_y = cbuffer.data(fp_geom_1010_off + 247 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_z = cbuffer.data(fp_geom_1010_off + 248 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_x = cbuffer.data(fp_geom_1010_off + 249 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_y = cbuffer.data(fp_geom_1010_off + 250 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_z = cbuffer.data(fp_geom_1010_off + 251 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_x = cbuffer.data(fp_geom_1010_off + 252 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_y = cbuffer.data(fp_geom_1010_off + 253 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_z = cbuffer.data(fp_geom_1010_off + 254 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_x = cbuffer.data(fp_geom_1010_off + 255 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_y = cbuffer.data(fp_geom_1010_off + 256 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_z = cbuffer.data(fp_geom_1010_off + 257 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_x = cbuffer.data(fp_geom_1010_off + 258 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_y = cbuffer.data(fp_geom_1010_off + 259 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_z = cbuffer.data(fp_geom_1010_off + 260 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_x = cbuffer.data(fp_geom_1010_off + 261 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_y = cbuffer.data(fp_geom_1010_off + 262 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_z = cbuffer.data(fp_geom_1010_off + 263 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_x = cbuffer.data(fp_geom_1010_off + 264 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_y = cbuffer.data(fp_geom_1010_off + 265 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_z = cbuffer.data(fp_geom_1010_off + 266 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_x = cbuffer.data(fp_geom_1010_off + 267 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_y = cbuffer.data(fp_geom_1010_off + 268 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_z = cbuffer.data(fp_geom_1010_off + 269 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_gsxx

            const auto gs_geom_1010_off = idx_geom_1010_gsxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxx_0 = cbuffer.data(gs_geom_1010_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_xxx_0, g_x_0_x_0_xxx_0, g_x_0_x_0_xxx_x, g_x_0_x_0_xxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxx_0[k] = -g_0_0_x_0_xxx_0[k] - g_x_0_x_0_xxx_0[k] * ab_x + g_x_0_x_0_xxx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxy_0 = cbuffer.data(gs_geom_1010_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxx_0, g_x_0_x_0_xxx_y, g_x_0_x_0_xxxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxy_0[k] = -g_x_0_x_0_xxx_0[k] * ab_y + g_x_0_x_0_xxx_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxz_0 = cbuffer.data(gs_geom_1010_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxx_0, g_x_0_x_0_xxx_z, g_x_0_x_0_xxxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxz_0[k] = -g_x_0_x_0_xxx_0[k] * ab_z + g_x_0_x_0_xxx_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxyy_0 = cbuffer.data(gs_geom_1010_off + 3 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxy_0, g_x_0_x_0_xxy_y, g_x_0_x_0_xxyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxyy_0[k] = -g_x_0_x_0_xxy_0[k] * ab_y + g_x_0_x_0_xxy_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxyz_0 = cbuffer.data(gs_geom_1010_off + 4 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxyz_0, g_x_0_x_0_xxz_0, g_x_0_x_0_xxz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxyz_0[k] = -g_x_0_x_0_xxz_0[k] * ab_y + g_x_0_x_0_xxz_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxzz_0 = cbuffer.data(gs_geom_1010_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxz_0, g_x_0_x_0_xxz_z, g_x_0_x_0_xxzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxzz_0[k] = -g_x_0_x_0_xxz_0[k] * ab_z + g_x_0_x_0_xxz_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyyy_0 = cbuffer.data(gs_geom_1010_off + 6 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyy_0, g_x_0_x_0_xyy_y, g_x_0_x_0_xyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyyy_0[k] = -g_x_0_x_0_xyy_0[k] * ab_y + g_x_0_x_0_xyy_y[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyyz_0 = cbuffer.data(gs_geom_1010_off + 7 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyyz_0, g_x_0_x_0_xyz_0, g_x_0_x_0_xyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyyz_0[k] = -g_x_0_x_0_xyz_0[k] * ab_y + g_x_0_x_0_xyz_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyzz_0 = cbuffer.data(gs_geom_1010_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyzz_0, g_x_0_x_0_xzz_0, g_x_0_x_0_xzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyzz_0[k] = -g_x_0_x_0_xzz_0[k] * ab_y + g_x_0_x_0_xzz_y[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xzzz_0 = cbuffer.data(gs_geom_1010_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xzz_0, g_x_0_x_0_xzz_z, g_x_0_x_0_xzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xzzz_0[k] = -g_x_0_x_0_xzz_0[k] * ab_z + g_x_0_x_0_xzz_z[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyyy_0 = cbuffer.data(gs_geom_1010_off + 10 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyy_0, g_x_0_x_0_yyy_y, g_x_0_x_0_yyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyyy_0[k] = -g_x_0_x_0_yyy_0[k] * ab_y + g_x_0_x_0_yyy_y[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyyz_0 = cbuffer.data(gs_geom_1010_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyyz_0, g_x_0_x_0_yyz_0, g_x_0_x_0_yyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyyz_0[k] = -g_x_0_x_0_yyz_0[k] * ab_y + g_x_0_x_0_yyz_y[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyzz_0 = cbuffer.data(gs_geom_1010_off + 12 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyzz_0, g_x_0_x_0_yzz_0, g_x_0_x_0_yzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyzz_0[k] = -g_x_0_x_0_yzz_0[k] * ab_y + g_x_0_x_0_yzz_y[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yzzz_0 = cbuffer.data(gs_geom_1010_off + 13 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yzzz_0, g_x_0_x_0_zzz_0, g_x_0_x_0_zzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yzzz_0[k] = -g_x_0_x_0_zzz_0[k] * ab_y + g_x_0_x_0_zzz_y[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_zzzz_0 = cbuffer.data(gs_geom_1010_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_zzz_0, g_x_0_x_0_zzz_z, g_x_0_x_0_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_zzzz_0[k] = -g_x_0_x_0_zzz_0[k] * ab_z + g_x_0_x_0_zzz_z[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxx_0 = cbuffer.data(gs_geom_1010_off + 15 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_xxx_0, g_x_0_y_0_xxx_0, g_x_0_y_0_xxx_x, g_x_0_y_0_xxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxx_0[k] = -g_0_0_y_0_xxx_0[k] - g_x_0_y_0_xxx_0[k] * ab_x + g_x_0_y_0_xxx_x[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxy_0 = cbuffer.data(gs_geom_1010_off + 16 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxx_0, g_x_0_y_0_xxx_y, g_x_0_y_0_xxxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxy_0[k] = -g_x_0_y_0_xxx_0[k] * ab_y + g_x_0_y_0_xxx_y[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxz_0 = cbuffer.data(gs_geom_1010_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxx_0, g_x_0_y_0_xxx_z, g_x_0_y_0_xxxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxz_0[k] = -g_x_0_y_0_xxx_0[k] * ab_z + g_x_0_y_0_xxx_z[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxyy_0 = cbuffer.data(gs_geom_1010_off + 18 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxy_0, g_x_0_y_0_xxy_y, g_x_0_y_0_xxyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxyy_0[k] = -g_x_0_y_0_xxy_0[k] * ab_y + g_x_0_y_0_xxy_y[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxyz_0 = cbuffer.data(gs_geom_1010_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxyz_0, g_x_0_y_0_xxz_0, g_x_0_y_0_xxz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxyz_0[k] = -g_x_0_y_0_xxz_0[k] * ab_y + g_x_0_y_0_xxz_y[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxzz_0 = cbuffer.data(gs_geom_1010_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxz_0, g_x_0_y_0_xxz_z, g_x_0_y_0_xxzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxzz_0[k] = -g_x_0_y_0_xxz_0[k] * ab_z + g_x_0_y_0_xxz_z[k];
            }

            /// Set up 21-22 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyyy_0 = cbuffer.data(gs_geom_1010_off + 21 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyy_0, g_x_0_y_0_xyy_y, g_x_0_y_0_xyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyyy_0[k] = -g_x_0_y_0_xyy_0[k] * ab_y + g_x_0_y_0_xyy_y[k];
            }

            /// Set up 22-23 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyyz_0 = cbuffer.data(gs_geom_1010_off + 22 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyyz_0, g_x_0_y_0_xyz_0, g_x_0_y_0_xyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyyz_0[k] = -g_x_0_y_0_xyz_0[k] * ab_y + g_x_0_y_0_xyz_y[k];
            }

            /// Set up 23-24 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyzz_0 = cbuffer.data(gs_geom_1010_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyzz_0, g_x_0_y_0_xzz_0, g_x_0_y_0_xzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyzz_0[k] = -g_x_0_y_0_xzz_0[k] * ab_y + g_x_0_y_0_xzz_y[k];
            }

            /// Set up 24-25 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xzzz_0 = cbuffer.data(gs_geom_1010_off + 24 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xzz_0, g_x_0_y_0_xzz_z, g_x_0_y_0_xzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xzzz_0[k] = -g_x_0_y_0_xzz_0[k] * ab_z + g_x_0_y_0_xzz_z[k];
            }

            /// Set up 25-26 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyyy_0 = cbuffer.data(gs_geom_1010_off + 25 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyy_0, g_x_0_y_0_yyy_y, g_x_0_y_0_yyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyyy_0[k] = -g_x_0_y_0_yyy_0[k] * ab_y + g_x_0_y_0_yyy_y[k];
            }

            /// Set up 26-27 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyyz_0 = cbuffer.data(gs_geom_1010_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyyz_0, g_x_0_y_0_yyz_0, g_x_0_y_0_yyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyyz_0[k] = -g_x_0_y_0_yyz_0[k] * ab_y + g_x_0_y_0_yyz_y[k];
            }

            /// Set up 27-28 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyzz_0 = cbuffer.data(gs_geom_1010_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyzz_0, g_x_0_y_0_yzz_0, g_x_0_y_0_yzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyzz_0[k] = -g_x_0_y_0_yzz_0[k] * ab_y + g_x_0_y_0_yzz_y[k];
            }

            /// Set up 28-29 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yzzz_0 = cbuffer.data(gs_geom_1010_off + 28 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yzzz_0, g_x_0_y_0_zzz_0, g_x_0_y_0_zzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yzzz_0[k] = -g_x_0_y_0_zzz_0[k] * ab_y + g_x_0_y_0_zzz_y[k];
            }

            /// Set up 29-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_zzzz_0 = cbuffer.data(gs_geom_1010_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_zzz_0, g_x_0_y_0_zzz_z, g_x_0_y_0_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_zzzz_0[k] = -g_x_0_y_0_zzz_0[k] * ab_z + g_x_0_y_0_zzz_z[k];
            }

            /// Set up 30-31 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxx_0 = cbuffer.data(gs_geom_1010_off + 30 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_xxx_0, g_x_0_z_0_xxx_0, g_x_0_z_0_xxx_x, g_x_0_z_0_xxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxx_0[k] = -g_0_0_z_0_xxx_0[k] - g_x_0_z_0_xxx_0[k] * ab_x + g_x_0_z_0_xxx_x[k];
            }

            /// Set up 31-32 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxy_0 = cbuffer.data(gs_geom_1010_off + 31 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxx_0, g_x_0_z_0_xxx_y, g_x_0_z_0_xxxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxy_0[k] = -g_x_0_z_0_xxx_0[k] * ab_y + g_x_0_z_0_xxx_y[k];
            }

            /// Set up 32-33 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxz_0 = cbuffer.data(gs_geom_1010_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxx_0, g_x_0_z_0_xxx_z, g_x_0_z_0_xxxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxz_0[k] = -g_x_0_z_0_xxx_0[k] * ab_z + g_x_0_z_0_xxx_z[k];
            }

            /// Set up 33-34 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxyy_0 = cbuffer.data(gs_geom_1010_off + 33 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxy_0, g_x_0_z_0_xxy_y, g_x_0_z_0_xxyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxyy_0[k] = -g_x_0_z_0_xxy_0[k] * ab_y + g_x_0_z_0_xxy_y[k];
            }

            /// Set up 34-35 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxyz_0 = cbuffer.data(gs_geom_1010_off + 34 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxyz_0, g_x_0_z_0_xxz_0, g_x_0_z_0_xxz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxyz_0[k] = -g_x_0_z_0_xxz_0[k] * ab_y + g_x_0_z_0_xxz_y[k];
            }

            /// Set up 35-36 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxzz_0 = cbuffer.data(gs_geom_1010_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxz_0, g_x_0_z_0_xxz_z, g_x_0_z_0_xxzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxzz_0[k] = -g_x_0_z_0_xxz_0[k] * ab_z + g_x_0_z_0_xxz_z[k];
            }

            /// Set up 36-37 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyyy_0 = cbuffer.data(gs_geom_1010_off + 36 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyy_0, g_x_0_z_0_xyy_y, g_x_0_z_0_xyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyyy_0[k] = -g_x_0_z_0_xyy_0[k] * ab_y + g_x_0_z_0_xyy_y[k];
            }

            /// Set up 37-38 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyyz_0 = cbuffer.data(gs_geom_1010_off + 37 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyyz_0, g_x_0_z_0_xyz_0, g_x_0_z_0_xyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyyz_0[k] = -g_x_0_z_0_xyz_0[k] * ab_y + g_x_0_z_0_xyz_y[k];
            }

            /// Set up 38-39 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyzz_0 = cbuffer.data(gs_geom_1010_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyzz_0, g_x_0_z_0_xzz_0, g_x_0_z_0_xzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyzz_0[k] = -g_x_0_z_0_xzz_0[k] * ab_y + g_x_0_z_0_xzz_y[k];
            }

            /// Set up 39-40 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xzzz_0 = cbuffer.data(gs_geom_1010_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xzz_0, g_x_0_z_0_xzz_z, g_x_0_z_0_xzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xzzz_0[k] = -g_x_0_z_0_xzz_0[k] * ab_z + g_x_0_z_0_xzz_z[k];
            }

            /// Set up 40-41 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyyy_0 = cbuffer.data(gs_geom_1010_off + 40 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyy_0, g_x_0_z_0_yyy_y, g_x_0_z_0_yyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyyy_0[k] = -g_x_0_z_0_yyy_0[k] * ab_y + g_x_0_z_0_yyy_y[k];
            }

            /// Set up 41-42 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyyz_0 = cbuffer.data(gs_geom_1010_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyyz_0, g_x_0_z_0_yyz_0, g_x_0_z_0_yyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyyz_0[k] = -g_x_0_z_0_yyz_0[k] * ab_y + g_x_0_z_0_yyz_y[k];
            }

            /// Set up 42-43 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyzz_0 = cbuffer.data(gs_geom_1010_off + 42 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyzz_0, g_x_0_z_0_yzz_0, g_x_0_z_0_yzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyzz_0[k] = -g_x_0_z_0_yzz_0[k] * ab_y + g_x_0_z_0_yzz_y[k];
            }

            /// Set up 43-44 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yzzz_0 = cbuffer.data(gs_geom_1010_off + 43 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yzzz_0, g_x_0_z_0_zzz_0, g_x_0_z_0_zzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yzzz_0[k] = -g_x_0_z_0_zzz_0[k] * ab_y + g_x_0_z_0_zzz_y[k];
            }

            /// Set up 44-45 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_zzzz_0 = cbuffer.data(gs_geom_1010_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_zzz_0, g_x_0_z_0_zzz_z, g_x_0_z_0_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_zzzz_0[k] = -g_x_0_z_0_zzz_0[k] * ab_z + g_x_0_z_0_zzz_z[k];
            }

            /// Set up 45-46 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxx_0 = cbuffer.data(gs_geom_1010_off + 45 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxx_0, g_y_0_x_0_xxx_x, g_y_0_x_0_xxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxx_0[k] = -g_y_0_x_0_xxx_0[k] * ab_x + g_y_0_x_0_xxx_x[k];
            }

            /// Set up 46-47 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxy_0 = cbuffer.data(gs_geom_1010_off + 46 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxy_0, g_y_0_x_0_xxy_0, g_y_0_x_0_xxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxy_0[k] = -g_y_0_x_0_xxy_0[k] * ab_x + g_y_0_x_0_xxy_x[k];
            }

            /// Set up 47-48 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxz_0 = cbuffer.data(gs_geom_1010_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxz_0, g_y_0_x_0_xxz_0, g_y_0_x_0_xxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxz_0[k] = -g_y_0_x_0_xxz_0[k] * ab_x + g_y_0_x_0_xxz_x[k];
            }

            /// Set up 48-49 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxyy_0 = cbuffer.data(gs_geom_1010_off + 48 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxyy_0, g_y_0_x_0_xyy_0, g_y_0_x_0_xyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxyy_0[k] = -g_y_0_x_0_xyy_0[k] * ab_x + g_y_0_x_0_xyy_x[k];
            }

            /// Set up 49-50 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxyz_0 = cbuffer.data(gs_geom_1010_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxyz_0, g_y_0_x_0_xyz_0, g_y_0_x_0_xyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxyz_0[k] = -g_y_0_x_0_xyz_0[k] * ab_x + g_y_0_x_0_xyz_x[k];
            }

            /// Set up 50-51 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxzz_0 = cbuffer.data(gs_geom_1010_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxzz_0, g_y_0_x_0_xzz_0, g_y_0_x_0_xzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxzz_0[k] = -g_y_0_x_0_xzz_0[k] * ab_x + g_y_0_x_0_xzz_x[k];
            }

            /// Set up 51-52 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyyy_0 = cbuffer.data(gs_geom_1010_off + 51 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyyy_0, g_y_0_x_0_yyy_0, g_y_0_x_0_yyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyyy_0[k] = -g_y_0_x_0_yyy_0[k] * ab_x + g_y_0_x_0_yyy_x[k];
            }

            /// Set up 52-53 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyyz_0 = cbuffer.data(gs_geom_1010_off + 52 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyyz_0, g_y_0_x_0_yyz_0, g_y_0_x_0_yyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyyz_0[k] = -g_y_0_x_0_yyz_0[k] * ab_x + g_y_0_x_0_yyz_x[k];
            }

            /// Set up 53-54 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyzz_0 = cbuffer.data(gs_geom_1010_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyzz_0, g_y_0_x_0_yzz_0, g_y_0_x_0_yzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyzz_0[k] = -g_y_0_x_0_yzz_0[k] * ab_x + g_y_0_x_0_yzz_x[k];
            }

            /// Set up 54-55 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xzzz_0 = cbuffer.data(gs_geom_1010_off + 54 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xzzz_0, g_y_0_x_0_zzz_0, g_y_0_x_0_zzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xzzz_0[k] = -g_y_0_x_0_zzz_0[k] * ab_x + g_y_0_x_0_zzz_x[k];
            }

            /// Set up 55-56 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyyy_0 = cbuffer.data(gs_geom_1010_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_yyy_0, g_y_0_x_0_yyy_0, g_y_0_x_0_yyy_y, g_y_0_x_0_yyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyyy_0[k] = -g_0_0_x_0_yyy_0[k] - g_y_0_x_0_yyy_0[k] * ab_y + g_y_0_x_0_yyy_y[k];
            }

            /// Set up 56-57 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyyz_0 = cbuffer.data(gs_geom_1010_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yyy_0, g_y_0_x_0_yyy_z, g_y_0_x_0_yyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyyz_0[k] = -g_y_0_x_0_yyy_0[k] * ab_z + g_y_0_x_0_yyy_z[k];
            }

            /// Set up 57-58 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyzz_0 = cbuffer.data(gs_geom_1010_off + 57 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yyz_0, g_y_0_x_0_yyz_z, g_y_0_x_0_yyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyzz_0[k] = -g_y_0_x_0_yyz_0[k] * ab_z + g_y_0_x_0_yyz_z[k];
            }

            /// Set up 58-59 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yzzz_0 = cbuffer.data(gs_geom_1010_off + 58 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yzz_0, g_y_0_x_0_yzz_z, g_y_0_x_0_yzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yzzz_0[k] = -g_y_0_x_0_yzz_0[k] * ab_z + g_y_0_x_0_yzz_z[k];
            }

            /// Set up 59-60 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_zzzz_0 = cbuffer.data(gs_geom_1010_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_zzz_0, g_y_0_x_0_zzz_z, g_y_0_x_0_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_zzzz_0[k] = -g_y_0_x_0_zzz_0[k] * ab_z + g_y_0_x_0_zzz_z[k];
            }

            /// Set up 60-61 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxx_0 = cbuffer.data(gs_geom_1010_off + 60 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxx_0, g_y_0_y_0_xxx_x, g_y_0_y_0_xxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxx_0[k] = -g_y_0_y_0_xxx_0[k] * ab_x + g_y_0_y_0_xxx_x[k];
            }

            /// Set up 61-62 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxy_0 = cbuffer.data(gs_geom_1010_off + 61 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxy_0, g_y_0_y_0_xxy_0, g_y_0_y_0_xxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxy_0[k] = -g_y_0_y_0_xxy_0[k] * ab_x + g_y_0_y_0_xxy_x[k];
            }

            /// Set up 62-63 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxz_0 = cbuffer.data(gs_geom_1010_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxz_0, g_y_0_y_0_xxz_0, g_y_0_y_0_xxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxz_0[k] = -g_y_0_y_0_xxz_0[k] * ab_x + g_y_0_y_0_xxz_x[k];
            }

            /// Set up 63-64 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxyy_0 = cbuffer.data(gs_geom_1010_off + 63 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxyy_0, g_y_0_y_0_xyy_0, g_y_0_y_0_xyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxyy_0[k] = -g_y_0_y_0_xyy_0[k] * ab_x + g_y_0_y_0_xyy_x[k];
            }

            /// Set up 64-65 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxyz_0 = cbuffer.data(gs_geom_1010_off + 64 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxyz_0, g_y_0_y_0_xyz_0, g_y_0_y_0_xyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxyz_0[k] = -g_y_0_y_0_xyz_0[k] * ab_x + g_y_0_y_0_xyz_x[k];
            }

            /// Set up 65-66 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxzz_0 = cbuffer.data(gs_geom_1010_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxzz_0, g_y_0_y_0_xzz_0, g_y_0_y_0_xzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxzz_0[k] = -g_y_0_y_0_xzz_0[k] * ab_x + g_y_0_y_0_xzz_x[k];
            }

            /// Set up 66-67 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyyy_0 = cbuffer.data(gs_geom_1010_off + 66 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyyy_0, g_y_0_y_0_yyy_0, g_y_0_y_0_yyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyyy_0[k] = -g_y_0_y_0_yyy_0[k] * ab_x + g_y_0_y_0_yyy_x[k];
            }

            /// Set up 67-68 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyyz_0 = cbuffer.data(gs_geom_1010_off + 67 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyyz_0, g_y_0_y_0_yyz_0, g_y_0_y_0_yyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyyz_0[k] = -g_y_0_y_0_yyz_0[k] * ab_x + g_y_0_y_0_yyz_x[k];
            }

            /// Set up 68-69 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyzz_0 = cbuffer.data(gs_geom_1010_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyzz_0, g_y_0_y_0_yzz_0, g_y_0_y_0_yzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyzz_0[k] = -g_y_0_y_0_yzz_0[k] * ab_x + g_y_0_y_0_yzz_x[k];
            }

            /// Set up 69-70 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xzzz_0 = cbuffer.data(gs_geom_1010_off + 69 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xzzz_0, g_y_0_y_0_zzz_0, g_y_0_y_0_zzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xzzz_0[k] = -g_y_0_y_0_zzz_0[k] * ab_x + g_y_0_y_0_zzz_x[k];
            }

            /// Set up 70-71 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyyy_0 = cbuffer.data(gs_geom_1010_off + 70 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_yyy_0, g_y_0_y_0_yyy_0, g_y_0_y_0_yyy_y, g_y_0_y_0_yyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyyy_0[k] = -g_0_0_y_0_yyy_0[k] - g_y_0_y_0_yyy_0[k] * ab_y + g_y_0_y_0_yyy_y[k];
            }

            /// Set up 71-72 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyyz_0 = cbuffer.data(gs_geom_1010_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yyy_0, g_y_0_y_0_yyy_z, g_y_0_y_0_yyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyyz_0[k] = -g_y_0_y_0_yyy_0[k] * ab_z + g_y_0_y_0_yyy_z[k];
            }

            /// Set up 72-73 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyzz_0 = cbuffer.data(gs_geom_1010_off + 72 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yyz_0, g_y_0_y_0_yyz_z, g_y_0_y_0_yyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyzz_0[k] = -g_y_0_y_0_yyz_0[k] * ab_z + g_y_0_y_0_yyz_z[k];
            }

            /// Set up 73-74 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yzzz_0 = cbuffer.data(gs_geom_1010_off + 73 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yzz_0, g_y_0_y_0_yzz_z, g_y_0_y_0_yzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yzzz_0[k] = -g_y_0_y_0_yzz_0[k] * ab_z + g_y_0_y_0_yzz_z[k];
            }

            /// Set up 74-75 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_zzzz_0 = cbuffer.data(gs_geom_1010_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_zzz_0, g_y_0_y_0_zzz_z, g_y_0_y_0_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_zzzz_0[k] = -g_y_0_y_0_zzz_0[k] * ab_z + g_y_0_y_0_zzz_z[k];
            }

            /// Set up 75-76 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxx_0 = cbuffer.data(gs_geom_1010_off + 75 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxx_0, g_y_0_z_0_xxx_x, g_y_0_z_0_xxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxx_0[k] = -g_y_0_z_0_xxx_0[k] * ab_x + g_y_0_z_0_xxx_x[k];
            }

            /// Set up 76-77 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxy_0 = cbuffer.data(gs_geom_1010_off + 76 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxy_0, g_y_0_z_0_xxy_0, g_y_0_z_0_xxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxy_0[k] = -g_y_0_z_0_xxy_0[k] * ab_x + g_y_0_z_0_xxy_x[k];
            }

            /// Set up 77-78 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxz_0 = cbuffer.data(gs_geom_1010_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxz_0, g_y_0_z_0_xxz_0, g_y_0_z_0_xxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxz_0[k] = -g_y_0_z_0_xxz_0[k] * ab_x + g_y_0_z_0_xxz_x[k];
            }

            /// Set up 78-79 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxyy_0 = cbuffer.data(gs_geom_1010_off + 78 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxyy_0, g_y_0_z_0_xyy_0, g_y_0_z_0_xyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxyy_0[k] = -g_y_0_z_0_xyy_0[k] * ab_x + g_y_0_z_0_xyy_x[k];
            }

            /// Set up 79-80 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxyz_0 = cbuffer.data(gs_geom_1010_off + 79 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxyz_0, g_y_0_z_0_xyz_0, g_y_0_z_0_xyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxyz_0[k] = -g_y_0_z_0_xyz_0[k] * ab_x + g_y_0_z_0_xyz_x[k];
            }

            /// Set up 80-81 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxzz_0 = cbuffer.data(gs_geom_1010_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxzz_0, g_y_0_z_0_xzz_0, g_y_0_z_0_xzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxzz_0[k] = -g_y_0_z_0_xzz_0[k] * ab_x + g_y_0_z_0_xzz_x[k];
            }

            /// Set up 81-82 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyyy_0 = cbuffer.data(gs_geom_1010_off + 81 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyyy_0, g_y_0_z_0_yyy_0, g_y_0_z_0_yyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyyy_0[k] = -g_y_0_z_0_yyy_0[k] * ab_x + g_y_0_z_0_yyy_x[k];
            }

            /// Set up 82-83 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyyz_0 = cbuffer.data(gs_geom_1010_off + 82 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyyz_0, g_y_0_z_0_yyz_0, g_y_0_z_0_yyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyyz_0[k] = -g_y_0_z_0_yyz_0[k] * ab_x + g_y_0_z_0_yyz_x[k];
            }

            /// Set up 83-84 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyzz_0 = cbuffer.data(gs_geom_1010_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyzz_0, g_y_0_z_0_yzz_0, g_y_0_z_0_yzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyzz_0[k] = -g_y_0_z_0_yzz_0[k] * ab_x + g_y_0_z_0_yzz_x[k];
            }

            /// Set up 84-85 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xzzz_0 = cbuffer.data(gs_geom_1010_off + 84 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xzzz_0, g_y_0_z_0_zzz_0, g_y_0_z_0_zzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xzzz_0[k] = -g_y_0_z_0_zzz_0[k] * ab_x + g_y_0_z_0_zzz_x[k];
            }

            /// Set up 85-86 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyyy_0 = cbuffer.data(gs_geom_1010_off + 85 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_yyy_0, g_y_0_z_0_yyy_0, g_y_0_z_0_yyy_y, g_y_0_z_0_yyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyyy_0[k] = -g_0_0_z_0_yyy_0[k] - g_y_0_z_0_yyy_0[k] * ab_y + g_y_0_z_0_yyy_y[k];
            }

            /// Set up 86-87 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyyz_0 = cbuffer.data(gs_geom_1010_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yyy_0, g_y_0_z_0_yyy_z, g_y_0_z_0_yyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyyz_0[k] = -g_y_0_z_0_yyy_0[k] * ab_z + g_y_0_z_0_yyy_z[k];
            }

            /// Set up 87-88 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyzz_0 = cbuffer.data(gs_geom_1010_off + 87 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yyz_0, g_y_0_z_0_yyz_z, g_y_0_z_0_yyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyzz_0[k] = -g_y_0_z_0_yyz_0[k] * ab_z + g_y_0_z_0_yyz_z[k];
            }

            /// Set up 88-89 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yzzz_0 = cbuffer.data(gs_geom_1010_off + 88 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yzz_0, g_y_0_z_0_yzz_z, g_y_0_z_0_yzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yzzz_0[k] = -g_y_0_z_0_yzz_0[k] * ab_z + g_y_0_z_0_yzz_z[k];
            }

            /// Set up 89-90 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_zzzz_0 = cbuffer.data(gs_geom_1010_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_zzz_0, g_y_0_z_0_zzz_z, g_y_0_z_0_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_zzzz_0[k] = -g_y_0_z_0_zzz_0[k] * ab_z + g_y_0_z_0_zzz_z[k];
            }

            /// Set up 90-91 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxx_0 = cbuffer.data(gs_geom_1010_off + 90 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxx_0, g_z_0_x_0_xxx_x, g_z_0_x_0_xxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxx_0[k] = -g_z_0_x_0_xxx_0[k] * ab_x + g_z_0_x_0_xxx_x[k];
            }

            /// Set up 91-92 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxy_0 = cbuffer.data(gs_geom_1010_off + 91 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxy_0, g_z_0_x_0_xxy_0, g_z_0_x_0_xxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxy_0[k] = -g_z_0_x_0_xxy_0[k] * ab_x + g_z_0_x_0_xxy_x[k];
            }

            /// Set up 92-93 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxz_0 = cbuffer.data(gs_geom_1010_off + 92 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxz_0, g_z_0_x_0_xxz_0, g_z_0_x_0_xxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxz_0[k] = -g_z_0_x_0_xxz_0[k] * ab_x + g_z_0_x_0_xxz_x[k];
            }

            /// Set up 93-94 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxyy_0 = cbuffer.data(gs_geom_1010_off + 93 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxyy_0, g_z_0_x_0_xyy_0, g_z_0_x_0_xyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxyy_0[k] = -g_z_0_x_0_xyy_0[k] * ab_x + g_z_0_x_0_xyy_x[k];
            }

            /// Set up 94-95 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxyz_0 = cbuffer.data(gs_geom_1010_off + 94 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxyz_0, g_z_0_x_0_xyz_0, g_z_0_x_0_xyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxyz_0[k] = -g_z_0_x_0_xyz_0[k] * ab_x + g_z_0_x_0_xyz_x[k];
            }

            /// Set up 95-96 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxzz_0 = cbuffer.data(gs_geom_1010_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxzz_0, g_z_0_x_0_xzz_0, g_z_0_x_0_xzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxzz_0[k] = -g_z_0_x_0_xzz_0[k] * ab_x + g_z_0_x_0_xzz_x[k];
            }

            /// Set up 96-97 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyyy_0 = cbuffer.data(gs_geom_1010_off + 96 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyyy_0, g_z_0_x_0_yyy_0, g_z_0_x_0_yyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyyy_0[k] = -g_z_0_x_0_yyy_0[k] * ab_x + g_z_0_x_0_yyy_x[k];
            }

            /// Set up 97-98 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyyz_0 = cbuffer.data(gs_geom_1010_off + 97 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyyz_0, g_z_0_x_0_yyz_0, g_z_0_x_0_yyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyyz_0[k] = -g_z_0_x_0_yyz_0[k] * ab_x + g_z_0_x_0_yyz_x[k];
            }

            /// Set up 98-99 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyzz_0 = cbuffer.data(gs_geom_1010_off + 98 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyzz_0, g_z_0_x_0_yzz_0, g_z_0_x_0_yzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyzz_0[k] = -g_z_0_x_0_yzz_0[k] * ab_x + g_z_0_x_0_yzz_x[k];
            }

            /// Set up 99-100 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xzzz_0 = cbuffer.data(gs_geom_1010_off + 99 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xzzz_0, g_z_0_x_0_zzz_0, g_z_0_x_0_zzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xzzz_0[k] = -g_z_0_x_0_zzz_0[k] * ab_x + g_z_0_x_0_zzz_x[k];
            }

            /// Set up 100-101 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyyy_0 = cbuffer.data(gs_geom_1010_off + 100 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyy_0, g_z_0_x_0_yyy_y, g_z_0_x_0_yyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyyy_0[k] = -g_z_0_x_0_yyy_0[k] * ab_y + g_z_0_x_0_yyy_y[k];
            }

            /// Set up 101-102 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyyz_0 = cbuffer.data(gs_geom_1010_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyyz_0, g_z_0_x_0_yyz_0, g_z_0_x_0_yyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyyz_0[k] = -g_z_0_x_0_yyz_0[k] * ab_y + g_z_0_x_0_yyz_y[k];
            }

            /// Set up 102-103 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyzz_0 = cbuffer.data(gs_geom_1010_off + 102 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyzz_0, g_z_0_x_0_yzz_0, g_z_0_x_0_yzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyzz_0[k] = -g_z_0_x_0_yzz_0[k] * ab_y + g_z_0_x_0_yzz_y[k];
            }

            /// Set up 103-104 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yzzz_0 = cbuffer.data(gs_geom_1010_off + 103 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yzzz_0, g_z_0_x_0_zzz_0, g_z_0_x_0_zzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yzzz_0[k] = -g_z_0_x_0_zzz_0[k] * ab_y + g_z_0_x_0_zzz_y[k];
            }

            /// Set up 104-105 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_zzzz_0 = cbuffer.data(gs_geom_1010_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_zzz_0, g_z_0_x_0_zzz_0, g_z_0_x_0_zzz_z, g_z_0_x_0_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_zzzz_0[k] = -g_0_0_x_0_zzz_0[k] - g_z_0_x_0_zzz_0[k] * ab_z + g_z_0_x_0_zzz_z[k];
            }

            /// Set up 105-106 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxx_0 = cbuffer.data(gs_geom_1010_off + 105 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxx_0, g_z_0_y_0_xxx_x, g_z_0_y_0_xxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxx_0[k] = -g_z_0_y_0_xxx_0[k] * ab_x + g_z_0_y_0_xxx_x[k];
            }

            /// Set up 106-107 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxy_0 = cbuffer.data(gs_geom_1010_off + 106 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxy_0, g_z_0_y_0_xxy_0, g_z_0_y_0_xxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxy_0[k] = -g_z_0_y_0_xxy_0[k] * ab_x + g_z_0_y_0_xxy_x[k];
            }

            /// Set up 107-108 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxz_0 = cbuffer.data(gs_geom_1010_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxz_0, g_z_0_y_0_xxz_0, g_z_0_y_0_xxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxz_0[k] = -g_z_0_y_0_xxz_0[k] * ab_x + g_z_0_y_0_xxz_x[k];
            }

            /// Set up 108-109 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxyy_0 = cbuffer.data(gs_geom_1010_off + 108 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxyy_0, g_z_0_y_0_xyy_0, g_z_0_y_0_xyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxyy_0[k] = -g_z_0_y_0_xyy_0[k] * ab_x + g_z_0_y_0_xyy_x[k];
            }

            /// Set up 109-110 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxyz_0 = cbuffer.data(gs_geom_1010_off + 109 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxyz_0, g_z_0_y_0_xyz_0, g_z_0_y_0_xyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxyz_0[k] = -g_z_0_y_0_xyz_0[k] * ab_x + g_z_0_y_0_xyz_x[k];
            }

            /// Set up 110-111 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxzz_0 = cbuffer.data(gs_geom_1010_off + 110 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxzz_0, g_z_0_y_0_xzz_0, g_z_0_y_0_xzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxzz_0[k] = -g_z_0_y_0_xzz_0[k] * ab_x + g_z_0_y_0_xzz_x[k];
            }

            /// Set up 111-112 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyyy_0 = cbuffer.data(gs_geom_1010_off + 111 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyyy_0, g_z_0_y_0_yyy_0, g_z_0_y_0_yyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyyy_0[k] = -g_z_0_y_0_yyy_0[k] * ab_x + g_z_0_y_0_yyy_x[k];
            }

            /// Set up 112-113 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyyz_0 = cbuffer.data(gs_geom_1010_off + 112 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyyz_0, g_z_0_y_0_yyz_0, g_z_0_y_0_yyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyyz_0[k] = -g_z_0_y_0_yyz_0[k] * ab_x + g_z_0_y_0_yyz_x[k];
            }

            /// Set up 113-114 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyzz_0 = cbuffer.data(gs_geom_1010_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyzz_0, g_z_0_y_0_yzz_0, g_z_0_y_0_yzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyzz_0[k] = -g_z_0_y_0_yzz_0[k] * ab_x + g_z_0_y_0_yzz_x[k];
            }

            /// Set up 114-115 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xzzz_0 = cbuffer.data(gs_geom_1010_off + 114 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xzzz_0, g_z_0_y_0_zzz_0, g_z_0_y_0_zzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xzzz_0[k] = -g_z_0_y_0_zzz_0[k] * ab_x + g_z_0_y_0_zzz_x[k];
            }

            /// Set up 115-116 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyyy_0 = cbuffer.data(gs_geom_1010_off + 115 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyy_0, g_z_0_y_0_yyy_y, g_z_0_y_0_yyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyyy_0[k] = -g_z_0_y_0_yyy_0[k] * ab_y + g_z_0_y_0_yyy_y[k];
            }

            /// Set up 116-117 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyyz_0 = cbuffer.data(gs_geom_1010_off + 116 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyyz_0, g_z_0_y_0_yyz_0, g_z_0_y_0_yyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyyz_0[k] = -g_z_0_y_0_yyz_0[k] * ab_y + g_z_0_y_0_yyz_y[k];
            }

            /// Set up 117-118 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyzz_0 = cbuffer.data(gs_geom_1010_off + 117 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyzz_0, g_z_0_y_0_yzz_0, g_z_0_y_0_yzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyzz_0[k] = -g_z_0_y_0_yzz_0[k] * ab_y + g_z_0_y_0_yzz_y[k];
            }

            /// Set up 118-119 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yzzz_0 = cbuffer.data(gs_geom_1010_off + 118 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yzzz_0, g_z_0_y_0_zzz_0, g_z_0_y_0_zzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yzzz_0[k] = -g_z_0_y_0_zzz_0[k] * ab_y + g_z_0_y_0_zzz_y[k];
            }

            /// Set up 119-120 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_zzzz_0 = cbuffer.data(gs_geom_1010_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_zzz_0, g_z_0_y_0_zzz_0, g_z_0_y_0_zzz_z, g_z_0_y_0_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_zzzz_0[k] = -g_0_0_y_0_zzz_0[k] - g_z_0_y_0_zzz_0[k] * ab_z + g_z_0_y_0_zzz_z[k];
            }

            /// Set up 120-121 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxx_0 = cbuffer.data(gs_geom_1010_off + 120 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxx_0, g_z_0_z_0_xxx_x, g_z_0_z_0_xxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxx_0[k] = -g_z_0_z_0_xxx_0[k] * ab_x + g_z_0_z_0_xxx_x[k];
            }

            /// Set up 121-122 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxy_0 = cbuffer.data(gs_geom_1010_off + 121 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxy_0, g_z_0_z_0_xxy_0, g_z_0_z_0_xxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxy_0[k] = -g_z_0_z_0_xxy_0[k] * ab_x + g_z_0_z_0_xxy_x[k];
            }

            /// Set up 122-123 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxz_0 = cbuffer.data(gs_geom_1010_off + 122 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxz_0, g_z_0_z_0_xxz_0, g_z_0_z_0_xxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxz_0[k] = -g_z_0_z_0_xxz_0[k] * ab_x + g_z_0_z_0_xxz_x[k];
            }

            /// Set up 123-124 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxyy_0 = cbuffer.data(gs_geom_1010_off + 123 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxyy_0, g_z_0_z_0_xyy_0, g_z_0_z_0_xyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxyy_0[k] = -g_z_0_z_0_xyy_0[k] * ab_x + g_z_0_z_0_xyy_x[k];
            }

            /// Set up 124-125 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxyz_0 = cbuffer.data(gs_geom_1010_off + 124 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxyz_0, g_z_0_z_0_xyz_0, g_z_0_z_0_xyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxyz_0[k] = -g_z_0_z_0_xyz_0[k] * ab_x + g_z_0_z_0_xyz_x[k];
            }

            /// Set up 125-126 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxzz_0 = cbuffer.data(gs_geom_1010_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxzz_0, g_z_0_z_0_xzz_0, g_z_0_z_0_xzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxzz_0[k] = -g_z_0_z_0_xzz_0[k] * ab_x + g_z_0_z_0_xzz_x[k];
            }

            /// Set up 126-127 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyyy_0 = cbuffer.data(gs_geom_1010_off + 126 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyyy_0, g_z_0_z_0_yyy_0, g_z_0_z_0_yyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyyy_0[k] = -g_z_0_z_0_yyy_0[k] * ab_x + g_z_0_z_0_yyy_x[k];
            }

            /// Set up 127-128 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyyz_0 = cbuffer.data(gs_geom_1010_off + 127 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyyz_0, g_z_0_z_0_yyz_0, g_z_0_z_0_yyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyyz_0[k] = -g_z_0_z_0_yyz_0[k] * ab_x + g_z_0_z_0_yyz_x[k];
            }

            /// Set up 128-129 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyzz_0 = cbuffer.data(gs_geom_1010_off + 128 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyzz_0, g_z_0_z_0_yzz_0, g_z_0_z_0_yzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyzz_0[k] = -g_z_0_z_0_yzz_0[k] * ab_x + g_z_0_z_0_yzz_x[k];
            }

            /// Set up 129-130 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xzzz_0 = cbuffer.data(gs_geom_1010_off + 129 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xzzz_0, g_z_0_z_0_zzz_0, g_z_0_z_0_zzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xzzz_0[k] = -g_z_0_z_0_zzz_0[k] * ab_x + g_z_0_z_0_zzz_x[k];
            }

            /// Set up 130-131 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyyy_0 = cbuffer.data(gs_geom_1010_off + 130 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyy_0, g_z_0_z_0_yyy_y, g_z_0_z_0_yyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyyy_0[k] = -g_z_0_z_0_yyy_0[k] * ab_y + g_z_0_z_0_yyy_y[k];
            }

            /// Set up 131-132 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyyz_0 = cbuffer.data(gs_geom_1010_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyyz_0, g_z_0_z_0_yyz_0, g_z_0_z_0_yyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyyz_0[k] = -g_z_0_z_0_yyz_0[k] * ab_y + g_z_0_z_0_yyz_y[k];
            }

            /// Set up 132-133 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyzz_0 = cbuffer.data(gs_geom_1010_off + 132 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyzz_0, g_z_0_z_0_yzz_0, g_z_0_z_0_yzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyzz_0[k] = -g_z_0_z_0_yzz_0[k] * ab_y + g_z_0_z_0_yzz_y[k];
            }

            /// Set up 133-134 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yzzz_0 = cbuffer.data(gs_geom_1010_off + 133 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yzzz_0, g_z_0_z_0_zzz_0, g_z_0_z_0_zzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yzzz_0[k] = -g_z_0_z_0_zzz_0[k] * ab_y + g_z_0_z_0_zzz_y[k];
            }

            /// Set up 134-135 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_zzzz_0 = cbuffer.data(gs_geom_1010_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_zzz_0, g_z_0_z_0_zzz_0, g_z_0_z_0_zzz_z, g_z_0_z_0_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_zzzz_0[k] = -g_0_0_z_0_zzz_0[k] - g_z_0_z_0_zzz_0[k] * ab_z + g_z_0_z_0_zzz_z[k];
            }
        }
    }
}

} // erirec namespace

