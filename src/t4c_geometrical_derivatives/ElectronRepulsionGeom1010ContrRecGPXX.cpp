#include "ElectronRepulsionGeom1010ContrRecGPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom1010_hrr_electron_repulsion_gpxx(CSimdArray<double>& cbuffer,
                                              const size_t idx_geom_1010_gpxx,
                                              const size_t idx_geom_0010_fpxx,
                                              const size_t idx_geom_1010_fpxx,
                                              const size_t idx_geom_1010_fdxx,
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

            const auto fp_geom_0010_off = idx_geom_0010_fpxx + i * dcomps + j;

            auto g_0_0_x_0_xxx_x = cbuffer.data(fp_geom_0010_off + 0 * ccomps * dcomps);

            auto g_0_0_x_0_xxx_y = cbuffer.data(fp_geom_0010_off + 1 * ccomps * dcomps);

            auto g_0_0_x_0_xxx_z = cbuffer.data(fp_geom_0010_off + 2 * ccomps * dcomps);

            auto g_0_0_x_0_xxy_x = cbuffer.data(fp_geom_0010_off + 3 * ccomps * dcomps);

            auto g_0_0_x_0_xxy_y = cbuffer.data(fp_geom_0010_off + 4 * ccomps * dcomps);

            auto g_0_0_x_0_xxy_z = cbuffer.data(fp_geom_0010_off + 5 * ccomps * dcomps);

            auto g_0_0_x_0_xxz_x = cbuffer.data(fp_geom_0010_off + 6 * ccomps * dcomps);

            auto g_0_0_x_0_xxz_y = cbuffer.data(fp_geom_0010_off + 7 * ccomps * dcomps);

            auto g_0_0_x_0_xxz_z = cbuffer.data(fp_geom_0010_off + 8 * ccomps * dcomps);

            auto g_0_0_x_0_xyy_x = cbuffer.data(fp_geom_0010_off + 9 * ccomps * dcomps);

            auto g_0_0_x_0_xyy_y = cbuffer.data(fp_geom_0010_off + 10 * ccomps * dcomps);

            auto g_0_0_x_0_xyy_z = cbuffer.data(fp_geom_0010_off + 11 * ccomps * dcomps);

            auto g_0_0_x_0_xyz_x = cbuffer.data(fp_geom_0010_off + 12 * ccomps * dcomps);

            auto g_0_0_x_0_xyz_y = cbuffer.data(fp_geom_0010_off + 13 * ccomps * dcomps);

            auto g_0_0_x_0_xyz_z = cbuffer.data(fp_geom_0010_off + 14 * ccomps * dcomps);

            auto g_0_0_x_0_xzz_x = cbuffer.data(fp_geom_0010_off + 15 * ccomps * dcomps);

            auto g_0_0_x_0_xzz_y = cbuffer.data(fp_geom_0010_off + 16 * ccomps * dcomps);

            auto g_0_0_x_0_xzz_z = cbuffer.data(fp_geom_0010_off + 17 * ccomps * dcomps);

            auto g_0_0_x_0_yyy_x = cbuffer.data(fp_geom_0010_off + 18 * ccomps * dcomps);

            auto g_0_0_x_0_yyy_y = cbuffer.data(fp_geom_0010_off + 19 * ccomps * dcomps);

            auto g_0_0_x_0_yyy_z = cbuffer.data(fp_geom_0010_off + 20 * ccomps * dcomps);

            auto g_0_0_x_0_yyz_x = cbuffer.data(fp_geom_0010_off + 21 * ccomps * dcomps);

            auto g_0_0_x_0_yyz_y = cbuffer.data(fp_geom_0010_off + 22 * ccomps * dcomps);

            auto g_0_0_x_0_yyz_z = cbuffer.data(fp_geom_0010_off + 23 * ccomps * dcomps);

            auto g_0_0_x_0_yzz_x = cbuffer.data(fp_geom_0010_off + 24 * ccomps * dcomps);

            auto g_0_0_x_0_yzz_y = cbuffer.data(fp_geom_0010_off + 25 * ccomps * dcomps);

            auto g_0_0_x_0_yzz_z = cbuffer.data(fp_geom_0010_off + 26 * ccomps * dcomps);

            auto g_0_0_x_0_zzz_x = cbuffer.data(fp_geom_0010_off + 27 * ccomps * dcomps);

            auto g_0_0_x_0_zzz_y = cbuffer.data(fp_geom_0010_off + 28 * ccomps * dcomps);

            auto g_0_0_x_0_zzz_z = cbuffer.data(fp_geom_0010_off + 29 * ccomps * dcomps);

            auto g_0_0_y_0_xxx_x = cbuffer.data(fp_geom_0010_off + 30 * ccomps * dcomps);

            auto g_0_0_y_0_xxx_y = cbuffer.data(fp_geom_0010_off + 31 * ccomps * dcomps);

            auto g_0_0_y_0_xxx_z = cbuffer.data(fp_geom_0010_off + 32 * ccomps * dcomps);

            auto g_0_0_y_0_xxy_x = cbuffer.data(fp_geom_0010_off + 33 * ccomps * dcomps);

            auto g_0_0_y_0_xxy_y = cbuffer.data(fp_geom_0010_off + 34 * ccomps * dcomps);

            auto g_0_0_y_0_xxy_z = cbuffer.data(fp_geom_0010_off + 35 * ccomps * dcomps);

            auto g_0_0_y_0_xxz_x = cbuffer.data(fp_geom_0010_off + 36 * ccomps * dcomps);

            auto g_0_0_y_0_xxz_y = cbuffer.data(fp_geom_0010_off + 37 * ccomps * dcomps);

            auto g_0_0_y_0_xxz_z = cbuffer.data(fp_geom_0010_off + 38 * ccomps * dcomps);

            auto g_0_0_y_0_xyy_x = cbuffer.data(fp_geom_0010_off + 39 * ccomps * dcomps);

            auto g_0_0_y_0_xyy_y = cbuffer.data(fp_geom_0010_off + 40 * ccomps * dcomps);

            auto g_0_0_y_0_xyy_z = cbuffer.data(fp_geom_0010_off + 41 * ccomps * dcomps);

            auto g_0_0_y_0_xyz_x = cbuffer.data(fp_geom_0010_off + 42 * ccomps * dcomps);

            auto g_0_0_y_0_xyz_y = cbuffer.data(fp_geom_0010_off + 43 * ccomps * dcomps);

            auto g_0_0_y_0_xyz_z = cbuffer.data(fp_geom_0010_off + 44 * ccomps * dcomps);

            auto g_0_0_y_0_xzz_x = cbuffer.data(fp_geom_0010_off + 45 * ccomps * dcomps);

            auto g_0_0_y_0_xzz_y = cbuffer.data(fp_geom_0010_off + 46 * ccomps * dcomps);

            auto g_0_0_y_0_xzz_z = cbuffer.data(fp_geom_0010_off + 47 * ccomps * dcomps);

            auto g_0_0_y_0_yyy_x = cbuffer.data(fp_geom_0010_off + 48 * ccomps * dcomps);

            auto g_0_0_y_0_yyy_y = cbuffer.data(fp_geom_0010_off + 49 * ccomps * dcomps);

            auto g_0_0_y_0_yyy_z = cbuffer.data(fp_geom_0010_off + 50 * ccomps * dcomps);

            auto g_0_0_y_0_yyz_x = cbuffer.data(fp_geom_0010_off + 51 * ccomps * dcomps);

            auto g_0_0_y_0_yyz_y = cbuffer.data(fp_geom_0010_off + 52 * ccomps * dcomps);

            auto g_0_0_y_0_yyz_z = cbuffer.data(fp_geom_0010_off + 53 * ccomps * dcomps);

            auto g_0_0_y_0_yzz_x = cbuffer.data(fp_geom_0010_off + 54 * ccomps * dcomps);

            auto g_0_0_y_0_yzz_y = cbuffer.data(fp_geom_0010_off + 55 * ccomps * dcomps);

            auto g_0_0_y_0_yzz_z = cbuffer.data(fp_geom_0010_off + 56 * ccomps * dcomps);

            auto g_0_0_y_0_zzz_x = cbuffer.data(fp_geom_0010_off + 57 * ccomps * dcomps);

            auto g_0_0_y_0_zzz_y = cbuffer.data(fp_geom_0010_off + 58 * ccomps * dcomps);

            auto g_0_0_y_0_zzz_z = cbuffer.data(fp_geom_0010_off + 59 * ccomps * dcomps);

            auto g_0_0_z_0_xxx_x = cbuffer.data(fp_geom_0010_off + 60 * ccomps * dcomps);

            auto g_0_0_z_0_xxx_y = cbuffer.data(fp_geom_0010_off + 61 * ccomps * dcomps);

            auto g_0_0_z_0_xxx_z = cbuffer.data(fp_geom_0010_off + 62 * ccomps * dcomps);

            auto g_0_0_z_0_xxy_x = cbuffer.data(fp_geom_0010_off + 63 * ccomps * dcomps);

            auto g_0_0_z_0_xxy_y = cbuffer.data(fp_geom_0010_off + 64 * ccomps * dcomps);

            auto g_0_0_z_0_xxy_z = cbuffer.data(fp_geom_0010_off + 65 * ccomps * dcomps);

            auto g_0_0_z_0_xxz_x = cbuffer.data(fp_geom_0010_off + 66 * ccomps * dcomps);

            auto g_0_0_z_0_xxz_y = cbuffer.data(fp_geom_0010_off + 67 * ccomps * dcomps);

            auto g_0_0_z_0_xxz_z = cbuffer.data(fp_geom_0010_off + 68 * ccomps * dcomps);

            auto g_0_0_z_0_xyy_x = cbuffer.data(fp_geom_0010_off + 69 * ccomps * dcomps);

            auto g_0_0_z_0_xyy_y = cbuffer.data(fp_geom_0010_off + 70 * ccomps * dcomps);

            auto g_0_0_z_0_xyy_z = cbuffer.data(fp_geom_0010_off + 71 * ccomps * dcomps);

            auto g_0_0_z_0_xyz_x = cbuffer.data(fp_geom_0010_off + 72 * ccomps * dcomps);

            auto g_0_0_z_0_xyz_y = cbuffer.data(fp_geom_0010_off + 73 * ccomps * dcomps);

            auto g_0_0_z_0_xyz_z = cbuffer.data(fp_geom_0010_off + 74 * ccomps * dcomps);

            auto g_0_0_z_0_xzz_x = cbuffer.data(fp_geom_0010_off + 75 * ccomps * dcomps);

            auto g_0_0_z_0_xzz_y = cbuffer.data(fp_geom_0010_off + 76 * ccomps * dcomps);

            auto g_0_0_z_0_xzz_z = cbuffer.data(fp_geom_0010_off + 77 * ccomps * dcomps);

            auto g_0_0_z_0_yyy_x = cbuffer.data(fp_geom_0010_off + 78 * ccomps * dcomps);

            auto g_0_0_z_0_yyy_y = cbuffer.data(fp_geom_0010_off + 79 * ccomps * dcomps);

            auto g_0_0_z_0_yyy_z = cbuffer.data(fp_geom_0010_off + 80 * ccomps * dcomps);

            auto g_0_0_z_0_yyz_x = cbuffer.data(fp_geom_0010_off + 81 * ccomps * dcomps);

            auto g_0_0_z_0_yyz_y = cbuffer.data(fp_geom_0010_off + 82 * ccomps * dcomps);

            auto g_0_0_z_0_yyz_z = cbuffer.data(fp_geom_0010_off + 83 * ccomps * dcomps);

            auto g_0_0_z_0_yzz_x = cbuffer.data(fp_geom_0010_off + 84 * ccomps * dcomps);

            auto g_0_0_z_0_yzz_y = cbuffer.data(fp_geom_0010_off + 85 * ccomps * dcomps);

            auto g_0_0_z_0_yzz_z = cbuffer.data(fp_geom_0010_off + 86 * ccomps * dcomps);

            auto g_0_0_z_0_zzz_x = cbuffer.data(fp_geom_0010_off + 87 * ccomps * dcomps);

            auto g_0_0_z_0_zzz_y = cbuffer.data(fp_geom_0010_off + 88 * ccomps * dcomps);

            auto g_0_0_z_0_zzz_z = cbuffer.data(fp_geom_0010_off + 89 * ccomps * dcomps);

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

            /// Set up components of auxilary buffer : FDSS

            const auto fd_geom_1010_off = idx_geom_1010_fdxx + i * dcomps + j;

            auto g_x_0_x_0_xxx_xx = cbuffer.data(fd_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xy = cbuffer.data(fd_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_xz = cbuffer.data(fd_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_yy = cbuffer.data(fd_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_yz = cbuffer.data(fd_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_xxx_zz = cbuffer.data(fd_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xx = cbuffer.data(fd_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xy = cbuffer.data(fd_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_xz = cbuffer.data(fd_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_yy = cbuffer.data(fd_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_yz = cbuffer.data(fd_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_xxy_zz = cbuffer.data(fd_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xx = cbuffer.data(fd_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xy = cbuffer.data(fd_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_xz = cbuffer.data(fd_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_yy = cbuffer.data(fd_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_yz = cbuffer.data(fd_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_x_0_xxz_zz = cbuffer.data(fd_geom_1010_off + 17 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xx = cbuffer.data(fd_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xy = cbuffer.data(fd_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_xz = cbuffer.data(fd_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_yy = cbuffer.data(fd_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_yz = cbuffer.data(fd_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_x_0_xyy_zz = cbuffer.data(fd_geom_1010_off + 23 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xx = cbuffer.data(fd_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xy = cbuffer.data(fd_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_xz = cbuffer.data(fd_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_yy = cbuffer.data(fd_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_yz = cbuffer.data(fd_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_x_0_xyz_zz = cbuffer.data(fd_geom_1010_off + 29 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xx = cbuffer.data(fd_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xy = cbuffer.data(fd_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_xz = cbuffer.data(fd_geom_1010_off + 32 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_yy = cbuffer.data(fd_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_yz = cbuffer.data(fd_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_x_0_xzz_zz = cbuffer.data(fd_geom_1010_off + 35 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xx = cbuffer.data(fd_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xy = cbuffer.data(fd_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_xz = cbuffer.data(fd_geom_1010_off + 38 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_yy = cbuffer.data(fd_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_yz = cbuffer.data(fd_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_x_0_yyy_zz = cbuffer.data(fd_geom_1010_off + 41 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xx = cbuffer.data(fd_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xy = cbuffer.data(fd_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_xz = cbuffer.data(fd_geom_1010_off + 44 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_yy = cbuffer.data(fd_geom_1010_off + 45 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_yz = cbuffer.data(fd_geom_1010_off + 46 * ccomps * dcomps);

            auto g_x_0_x_0_yyz_zz = cbuffer.data(fd_geom_1010_off + 47 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xx = cbuffer.data(fd_geom_1010_off + 48 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xy = cbuffer.data(fd_geom_1010_off + 49 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_xz = cbuffer.data(fd_geom_1010_off + 50 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_yy = cbuffer.data(fd_geom_1010_off + 51 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_yz = cbuffer.data(fd_geom_1010_off + 52 * ccomps * dcomps);

            auto g_x_0_x_0_yzz_zz = cbuffer.data(fd_geom_1010_off + 53 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xx = cbuffer.data(fd_geom_1010_off + 54 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xy = cbuffer.data(fd_geom_1010_off + 55 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_xz = cbuffer.data(fd_geom_1010_off + 56 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_yy = cbuffer.data(fd_geom_1010_off + 57 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_yz = cbuffer.data(fd_geom_1010_off + 58 * ccomps * dcomps);

            auto g_x_0_x_0_zzz_zz = cbuffer.data(fd_geom_1010_off + 59 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xx = cbuffer.data(fd_geom_1010_off + 60 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xy = cbuffer.data(fd_geom_1010_off + 61 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_xz = cbuffer.data(fd_geom_1010_off + 62 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_yy = cbuffer.data(fd_geom_1010_off + 63 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_yz = cbuffer.data(fd_geom_1010_off + 64 * ccomps * dcomps);

            auto g_x_0_y_0_xxx_zz = cbuffer.data(fd_geom_1010_off + 65 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xx = cbuffer.data(fd_geom_1010_off + 66 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xy = cbuffer.data(fd_geom_1010_off + 67 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_xz = cbuffer.data(fd_geom_1010_off + 68 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_yy = cbuffer.data(fd_geom_1010_off + 69 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_yz = cbuffer.data(fd_geom_1010_off + 70 * ccomps * dcomps);

            auto g_x_0_y_0_xxy_zz = cbuffer.data(fd_geom_1010_off + 71 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xx = cbuffer.data(fd_geom_1010_off + 72 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xy = cbuffer.data(fd_geom_1010_off + 73 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_xz = cbuffer.data(fd_geom_1010_off + 74 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_yy = cbuffer.data(fd_geom_1010_off + 75 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_yz = cbuffer.data(fd_geom_1010_off + 76 * ccomps * dcomps);

            auto g_x_0_y_0_xxz_zz = cbuffer.data(fd_geom_1010_off + 77 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xx = cbuffer.data(fd_geom_1010_off + 78 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xy = cbuffer.data(fd_geom_1010_off + 79 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_xz = cbuffer.data(fd_geom_1010_off + 80 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_yy = cbuffer.data(fd_geom_1010_off + 81 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_yz = cbuffer.data(fd_geom_1010_off + 82 * ccomps * dcomps);

            auto g_x_0_y_0_xyy_zz = cbuffer.data(fd_geom_1010_off + 83 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xx = cbuffer.data(fd_geom_1010_off + 84 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xy = cbuffer.data(fd_geom_1010_off + 85 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_xz = cbuffer.data(fd_geom_1010_off + 86 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_yy = cbuffer.data(fd_geom_1010_off + 87 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_yz = cbuffer.data(fd_geom_1010_off + 88 * ccomps * dcomps);

            auto g_x_0_y_0_xyz_zz = cbuffer.data(fd_geom_1010_off + 89 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xx = cbuffer.data(fd_geom_1010_off + 90 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xy = cbuffer.data(fd_geom_1010_off + 91 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_xz = cbuffer.data(fd_geom_1010_off + 92 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_yy = cbuffer.data(fd_geom_1010_off + 93 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_yz = cbuffer.data(fd_geom_1010_off + 94 * ccomps * dcomps);

            auto g_x_0_y_0_xzz_zz = cbuffer.data(fd_geom_1010_off + 95 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xx = cbuffer.data(fd_geom_1010_off + 96 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xy = cbuffer.data(fd_geom_1010_off + 97 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_xz = cbuffer.data(fd_geom_1010_off + 98 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_yy = cbuffer.data(fd_geom_1010_off + 99 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_yz = cbuffer.data(fd_geom_1010_off + 100 * ccomps * dcomps);

            auto g_x_0_y_0_yyy_zz = cbuffer.data(fd_geom_1010_off + 101 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xx = cbuffer.data(fd_geom_1010_off + 102 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xy = cbuffer.data(fd_geom_1010_off + 103 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_xz = cbuffer.data(fd_geom_1010_off + 104 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_yy = cbuffer.data(fd_geom_1010_off + 105 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_yz = cbuffer.data(fd_geom_1010_off + 106 * ccomps * dcomps);

            auto g_x_0_y_0_yyz_zz = cbuffer.data(fd_geom_1010_off + 107 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xx = cbuffer.data(fd_geom_1010_off + 108 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xy = cbuffer.data(fd_geom_1010_off + 109 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_xz = cbuffer.data(fd_geom_1010_off + 110 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_yy = cbuffer.data(fd_geom_1010_off + 111 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_yz = cbuffer.data(fd_geom_1010_off + 112 * ccomps * dcomps);

            auto g_x_0_y_0_yzz_zz = cbuffer.data(fd_geom_1010_off + 113 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xx = cbuffer.data(fd_geom_1010_off + 114 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xy = cbuffer.data(fd_geom_1010_off + 115 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_xz = cbuffer.data(fd_geom_1010_off + 116 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_yy = cbuffer.data(fd_geom_1010_off + 117 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_yz = cbuffer.data(fd_geom_1010_off + 118 * ccomps * dcomps);

            auto g_x_0_y_0_zzz_zz = cbuffer.data(fd_geom_1010_off + 119 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xx = cbuffer.data(fd_geom_1010_off + 120 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xy = cbuffer.data(fd_geom_1010_off + 121 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_xz = cbuffer.data(fd_geom_1010_off + 122 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_yy = cbuffer.data(fd_geom_1010_off + 123 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_yz = cbuffer.data(fd_geom_1010_off + 124 * ccomps * dcomps);

            auto g_x_0_z_0_xxx_zz = cbuffer.data(fd_geom_1010_off + 125 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xx = cbuffer.data(fd_geom_1010_off + 126 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xy = cbuffer.data(fd_geom_1010_off + 127 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_xz = cbuffer.data(fd_geom_1010_off + 128 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_yy = cbuffer.data(fd_geom_1010_off + 129 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_yz = cbuffer.data(fd_geom_1010_off + 130 * ccomps * dcomps);

            auto g_x_0_z_0_xxy_zz = cbuffer.data(fd_geom_1010_off + 131 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xx = cbuffer.data(fd_geom_1010_off + 132 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xy = cbuffer.data(fd_geom_1010_off + 133 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_xz = cbuffer.data(fd_geom_1010_off + 134 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_yy = cbuffer.data(fd_geom_1010_off + 135 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_yz = cbuffer.data(fd_geom_1010_off + 136 * ccomps * dcomps);

            auto g_x_0_z_0_xxz_zz = cbuffer.data(fd_geom_1010_off + 137 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xx = cbuffer.data(fd_geom_1010_off + 138 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xy = cbuffer.data(fd_geom_1010_off + 139 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_xz = cbuffer.data(fd_geom_1010_off + 140 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_yy = cbuffer.data(fd_geom_1010_off + 141 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_yz = cbuffer.data(fd_geom_1010_off + 142 * ccomps * dcomps);

            auto g_x_0_z_0_xyy_zz = cbuffer.data(fd_geom_1010_off + 143 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xx = cbuffer.data(fd_geom_1010_off + 144 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xy = cbuffer.data(fd_geom_1010_off + 145 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_xz = cbuffer.data(fd_geom_1010_off + 146 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_yy = cbuffer.data(fd_geom_1010_off + 147 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_yz = cbuffer.data(fd_geom_1010_off + 148 * ccomps * dcomps);

            auto g_x_0_z_0_xyz_zz = cbuffer.data(fd_geom_1010_off + 149 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xx = cbuffer.data(fd_geom_1010_off + 150 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xy = cbuffer.data(fd_geom_1010_off + 151 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_xz = cbuffer.data(fd_geom_1010_off + 152 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_yy = cbuffer.data(fd_geom_1010_off + 153 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_yz = cbuffer.data(fd_geom_1010_off + 154 * ccomps * dcomps);

            auto g_x_0_z_0_xzz_zz = cbuffer.data(fd_geom_1010_off + 155 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xx = cbuffer.data(fd_geom_1010_off + 156 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xy = cbuffer.data(fd_geom_1010_off + 157 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_xz = cbuffer.data(fd_geom_1010_off + 158 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_yy = cbuffer.data(fd_geom_1010_off + 159 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_yz = cbuffer.data(fd_geom_1010_off + 160 * ccomps * dcomps);

            auto g_x_0_z_0_yyy_zz = cbuffer.data(fd_geom_1010_off + 161 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xx = cbuffer.data(fd_geom_1010_off + 162 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xy = cbuffer.data(fd_geom_1010_off + 163 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_xz = cbuffer.data(fd_geom_1010_off + 164 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_yy = cbuffer.data(fd_geom_1010_off + 165 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_yz = cbuffer.data(fd_geom_1010_off + 166 * ccomps * dcomps);

            auto g_x_0_z_0_yyz_zz = cbuffer.data(fd_geom_1010_off + 167 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xx = cbuffer.data(fd_geom_1010_off + 168 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xy = cbuffer.data(fd_geom_1010_off + 169 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_xz = cbuffer.data(fd_geom_1010_off + 170 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_yy = cbuffer.data(fd_geom_1010_off + 171 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_yz = cbuffer.data(fd_geom_1010_off + 172 * ccomps * dcomps);

            auto g_x_0_z_0_yzz_zz = cbuffer.data(fd_geom_1010_off + 173 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xx = cbuffer.data(fd_geom_1010_off + 174 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xy = cbuffer.data(fd_geom_1010_off + 175 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_xz = cbuffer.data(fd_geom_1010_off + 176 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_yy = cbuffer.data(fd_geom_1010_off + 177 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_yz = cbuffer.data(fd_geom_1010_off + 178 * ccomps * dcomps);

            auto g_x_0_z_0_zzz_zz = cbuffer.data(fd_geom_1010_off + 179 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xx = cbuffer.data(fd_geom_1010_off + 180 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xy = cbuffer.data(fd_geom_1010_off + 181 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_xz = cbuffer.data(fd_geom_1010_off + 182 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_yy = cbuffer.data(fd_geom_1010_off + 183 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_yz = cbuffer.data(fd_geom_1010_off + 184 * ccomps * dcomps);

            auto g_y_0_x_0_xxx_zz = cbuffer.data(fd_geom_1010_off + 185 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xx = cbuffer.data(fd_geom_1010_off + 186 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xy = cbuffer.data(fd_geom_1010_off + 187 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_xz = cbuffer.data(fd_geom_1010_off + 188 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_yy = cbuffer.data(fd_geom_1010_off + 189 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_yz = cbuffer.data(fd_geom_1010_off + 190 * ccomps * dcomps);

            auto g_y_0_x_0_xxy_zz = cbuffer.data(fd_geom_1010_off + 191 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xx = cbuffer.data(fd_geom_1010_off + 192 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xy = cbuffer.data(fd_geom_1010_off + 193 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_xz = cbuffer.data(fd_geom_1010_off + 194 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_yy = cbuffer.data(fd_geom_1010_off + 195 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_yz = cbuffer.data(fd_geom_1010_off + 196 * ccomps * dcomps);

            auto g_y_0_x_0_xxz_zz = cbuffer.data(fd_geom_1010_off + 197 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xx = cbuffer.data(fd_geom_1010_off + 198 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xy = cbuffer.data(fd_geom_1010_off + 199 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_xz = cbuffer.data(fd_geom_1010_off + 200 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_yy = cbuffer.data(fd_geom_1010_off + 201 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_yz = cbuffer.data(fd_geom_1010_off + 202 * ccomps * dcomps);

            auto g_y_0_x_0_xyy_zz = cbuffer.data(fd_geom_1010_off + 203 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xx = cbuffer.data(fd_geom_1010_off + 204 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xy = cbuffer.data(fd_geom_1010_off + 205 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_xz = cbuffer.data(fd_geom_1010_off + 206 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_yy = cbuffer.data(fd_geom_1010_off + 207 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_yz = cbuffer.data(fd_geom_1010_off + 208 * ccomps * dcomps);

            auto g_y_0_x_0_xyz_zz = cbuffer.data(fd_geom_1010_off + 209 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xx = cbuffer.data(fd_geom_1010_off + 210 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xy = cbuffer.data(fd_geom_1010_off + 211 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_xz = cbuffer.data(fd_geom_1010_off + 212 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_yy = cbuffer.data(fd_geom_1010_off + 213 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_yz = cbuffer.data(fd_geom_1010_off + 214 * ccomps * dcomps);

            auto g_y_0_x_0_xzz_zz = cbuffer.data(fd_geom_1010_off + 215 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xx = cbuffer.data(fd_geom_1010_off + 216 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xy = cbuffer.data(fd_geom_1010_off + 217 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_xz = cbuffer.data(fd_geom_1010_off + 218 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_yy = cbuffer.data(fd_geom_1010_off + 219 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_yz = cbuffer.data(fd_geom_1010_off + 220 * ccomps * dcomps);

            auto g_y_0_x_0_yyy_zz = cbuffer.data(fd_geom_1010_off + 221 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xx = cbuffer.data(fd_geom_1010_off + 222 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xy = cbuffer.data(fd_geom_1010_off + 223 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_xz = cbuffer.data(fd_geom_1010_off + 224 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_yy = cbuffer.data(fd_geom_1010_off + 225 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_yz = cbuffer.data(fd_geom_1010_off + 226 * ccomps * dcomps);

            auto g_y_0_x_0_yyz_zz = cbuffer.data(fd_geom_1010_off + 227 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xx = cbuffer.data(fd_geom_1010_off + 228 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xy = cbuffer.data(fd_geom_1010_off + 229 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_xz = cbuffer.data(fd_geom_1010_off + 230 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_yy = cbuffer.data(fd_geom_1010_off + 231 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_yz = cbuffer.data(fd_geom_1010_off + 232 * ccomps * dcomps);

            auto g_y_0_x_0_yzz_zz = cbuffer.data(fd_geom_1010_off + 233 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xx = cbuffer.data(fd_geom_1010_off + 234 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xy = cbuffer.data(fd_geom_1010_off + 235 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_xz = cbuffer.data(fd_geom_1010_off + 236 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_yy = cbuffer.data(fd_geom_1010_off + 237 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_yz = cbuffer.data(fd_geom_1010_off + 238 * ccomps * dcomps);

            auto g_y_0_x_0_zzz_zz = cbuffer.data(fd_geom_1010_off + 239 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xx = cbuffer.data(fd_geom_1010_off + 240 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xy = cbuffer.data(fd_geom_1010_off + 241 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_xz = cbuffer.data(fd_geom_1010_off + 242 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_yy = cbuffer.data(fd_geom_1010_off + 243 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_yz = cbuffer.data(fd_geom_1010_off + 244 * ccomps * dcomps);

            auto g_y_0_y_0_xxx_zz = cbuffer.data(fd_geom_1010_off + 245 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xx = cbuffer.data(fd_geom_1010_off + 246 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xy = cbuffer.data(fd_geom_1010_off + 247 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_xz = cbuffer.data(fd_geom_1010_off + 248 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_yy = cbuffer.data(fd_geom_1010_off + 249 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_yz = cbuffer.data(fd_geom_1010_off + 250 * ccomps * dcomps);

            auto g_y_0_y_0_xxy_zz = cbuffer.data(fd_geom_1010_off + 251 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xx = cbuffer.data(fd_geom_1010_off + 252 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xy = cbuffer.data(fd_geom_1010_off + 253 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_xz = cbuffer.data(fd_geom_1010_off + 254 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_yy = cbuffer.data(fd_geom_1010_off + 255 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_yz = cbuffer.data(fd_geom_1010_off + 256 * ccomps * dcomps);

            auto g_y_0_y_0_xxz_zz = cbuffer.data(fd_geom_1010_off + 257 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xx = cbuffer.data(fd_geom_1010_off + 258 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xy = cbuffer.data(fd_geom_1010_off + 259 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_xz = cbuffer.data(fd_geom_1010_off + 260 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_yy = cbuffer.data(fd_geom_1010_off + 261 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_yz = cbuffer.data(fd_geom_1010_off + 262 * ccomps * dcomps);

            auto g_y_0_y_0_xyy_zz = cbuffer.data(fd_geom_1010_off + 263 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xx = cbuffer.data(fd_geom_1010_off + 264 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xy = cbuffer.data(fd_geom_1010_off + 265 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_xz = cbuffer.data(fd_geom_1010_off + 266 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_yy = cbuffer.data(fd_geom_1010_off + 267 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_yz = cbuffer.data(fd_geom_1010_off + 268 * ccomps * dcomps);

            auto g_y_0_y_0_xyz_zz = cbuffer.data(fd_geom_1010_off + 269 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xx = cbuffer.data(fd_geom_1010_off + 270 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xy = cbuffer.data(fd_geom_1010_off + 271 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_xz = cbuffer.data(fd_geom_1010_off + 272 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_yy = cbuffer.data(fd_geom_1010_off + 273 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_yz = cbuffer.data(fd_geom_1010_off + 274 * ccomps * dcomps);

            auto g_y_0_y_0_xzz_zz = cbuffer.data(fd_geom_1010_off + 275 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xx = cbuffer.data(fd_geom_1010_off + 276 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xy = cbuffer.data(fd_geom_1010_off + 277 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_xz = cbuffer.data(fd_geom_1010_off + 278 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_yy = cbuffer.data(fd_geom_1010_off + 279 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_yz = cbuffer.data(fd_geom_1010_off + 280 * ccomps * dcomps);

            auto g_y_0_y_0_yyy_zz = cbuffer.data(fd_geom_1010_off + 281 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xx = cbuffer.data(fd_geom_1010_off + 282 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xy = cbuffer.data(fd_geom_1010_off + 283 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_xz = cbuffer.data(fd_geom_1010_off + 284 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_yy = cbuffer.data(fd_geom_1010_off + 285 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_yz = cbuffer.data(fd_geom_1010_off + 286 * ccomps * dcomps);

            auto g_y_0_y_0_yyz_zz = cbuffer.data(fd_geom_1010_off + 287 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xx = cbuffer.data(fd_geom_1010_off + 288 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xy = cbuffer.data(fd_geom_1010_off + 289 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_xz = cbuffer.data(fd_geom_1010_off + 290 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_yy = cbuffer.data(fd_geom_1010_off + 291 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_yz = cbuffer.data(fd_geom_1010_off + 292 * ccomps * dcomps);

            auto g_y_0_y_0_yzz_zz = cbuffer.data(fd_geom_1010_off + 293 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xx = cbuffer.data(fd_geom_1010_off + 294 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xy = cbuffer.data(fd_geom_1010_off + 295 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_xz = cbuffer.data(fd_geom_1010_off + 296 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_yy = cbuffer.data(fd_geom_1010_off + 297 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_yz = cbuffer.data(fd_geom_1010_off + 298 * ccomps * dcomps);

            auto g_y_0_y_0_zzz_zz = cbuffer.data(fd_geom_1010_off + 299 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xx = cbuffer.data(fd_geom_1010_off + 300 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xy = cbuffer.data(fd_geom_1010_off + 301 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_xz = cbuffer.data(fd_geom_1010_off + 302 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_yy = cbuffer.data(fd_geom_1010_off + 303 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_yz = cbuffer.data(fd_geom_1010_off + 304 * ccomps * dcomps);

            auto g_y_0_z_0_xxx_zz = cbuffer.data(fd_geom_1010_off + 305 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xx = cbuffer.data(fd_geom_1010_off + 306 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xy = cbuffer.data(fd_geom_1010_off + 307 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_xz = cbuffer.data(fd_geom_1010_off + 308 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_yy = cbuffer.data(fd_geom_1010_off + 309 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_yz = cbuffer.data(fd_geom_1010_off + 310 * ccomps * dcomps);

            auto g_y_0_z_0_xxy_zz = cbuffer.data(fd_geom_1010_off + 311 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xx = cbuffer.data(fd_geom_1010_off + 312 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xy = cbuffer.data(fd_geom_1010_off + 313 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_xz = cbuffer.data(fd_geom_1010_off + 314 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_yy = cbuffer.data(fd_geom_1010_off + 315 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_yz = cbuffer.data(fd_geom_1010_off + 316 * ccomps * dcomps);

            auto g_y_0_z_0_xxz_zz = cbuffer.data(fd_geom_1010_off + 317 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xx = cbuffer.data(fd_geom_1010_off + 318 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xy = cbuffer.data(fd_geom_1010_off + 319 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_xz = cbuffer.data(fd_geom_1010_off + 320 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_yy = cbuffer.data(fd_geom_1010_off + 321 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_yz = cbuffer.data(fd_geom_1010_off + 322 * ccomps * dcomps);

            auto g_y_0_z_0_xyy_zz = cbuffer.data(fd_geom_1010_off + 323 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xx = cbuffer.data(fd_geom_1010_off + 324 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xy = cbuffer.data(fd_geom_1010_off + 325 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_xz = cbuffer.data(fd_geom_1010_off + 326 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_yy = cbuffer.data(fd_geom_1010_off + 327 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_yz = cbuffer.data(fd_geom_1010_off + 328 * ccomps * dcomps);

            auto g_y_0_z_0_xyz_zz = cbuffer.data(fd_geom_1010_off + 329 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xx = cbuffer.data(fd_geom_1010_off + 330 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xy = cbuffer.data(fd_geom_1010_off + 331 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_xz = cbuffer.data(fd_geom_1010_off + 332 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_yy = cbuffer.data(fd_geom_1010_off + 333 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_yz = cbuffer.data(fd_geom_1010_off + 334 * ccomps * dcomps);

            auto g_y_0_z_0_xzz_zz = cbuffer.data(fd_geom_1010_off + 335 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xx = cbuffer.data(fd_geom_1010_off + 336 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xy = cbuffer.data(fd_geom_1010_off + 337 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_xz = cbuffer.data(fd_geom_1010_off + 338 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_yy = cbuffer.data(fd_geom_1010_off + 339 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_yz = cbuffer.data(fd_geom_1010_off + 340 * ccomps * dcomps);

            auto g_y_0_z_0_yyy_zz = cbuffer.data(fd_geom_1010_off + 341 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xx = cbuffer.data(fd_geom_1010_off + 342 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xy = cbuffer.data(fd_geom_1010_off + 343 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_xz = cbuffer.data(fd_geom_1010_off + 344 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_yy = cbuffer.data(fd_geom_1010_off + 345 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_yz = cbuffer.data(fd_geom_1010_off + 346 * ccomps * dcomps);

            auto g_y_0_z_0_yyz_zz = cbuffer.data(fd_geom_1010_off + 347 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xx = cbuffer.data(fd_geom_1010_off + 348 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xy = cbuffer.data(fd_geom_1010_off + 349 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_xz = cbuffer.data(fd_geom_1010_off + 350 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_yy = cbuffer.data(fd_geom_1010_off + 351 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_yz = cbuffer.data(fd_geom_1010_off + 352 * ccomps * dcomps);

            auto g_y_0_z_0_yzz_zz = cbuffer.data(fd_geom_1010_off + 353 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xx = cbuffer.data(fd_geom_1010_off + 354 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xy = cbuffer.data(fd_geom_1010_off + 355 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_xz = cbuffer.data(fd_geom_1010_off + 356 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_yy = cbuffer.data(fd_geom_1010_off + 357 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_yz = cbuffer.data(fd_geom_1010_off + 358 * ccomps * dcomps);

            auto g_y_0_z_0_zzz_zz = cbuffer.data(fd_geom_1010_off + 359 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xx = cbuffer.data(fd_geom_1010_off + 360 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xy = cbuffer.data(fd_geom_1010_off + 361 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_xz = cbuffer.data(fd_geom_1010_off + 362 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_yy = cbuffer.data(fd_geom_1010_off + 363 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_yz = cbuffer.data(fd_geom_1010_off + 364 * ccomps * dcomps);

            auto g_z_0_x_0_xxx_zz = cbuffer.data(fd_geom_1010_off + 365 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xx = cbuffer.data(fd_geom_1010_off + 366 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xy = cbuffer.data(fd_geom_1010_off + 367 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_xz = cbuffer.data(fd_geom_1010_off + 368 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_yy = cbuffer.data(fd_geom_1010_off + 369 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_yz = cbuffer.data(fd_geom_1010_off + 370 * ccomps * dcomps);

            auto g_z_0_x_0_xxy_zz = cbuffer.data(fd_geom_1010_off + 371 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xx = cbuffer.data(fd_geom_1010_off + 372 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xy = cbuffer.data(fd_geom_1010_off + 373 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_xz = cbuffer.data(fd_geom_1010_off + 374 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_yy = cbuffer.data(fd_geom_1010_off + 375 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_yz = cbuffer.data(fd_geom_1010_off + 376 * ccomps * dcomps);

            auto g_z_0_x_0_xxz_zz = cbuffer.data(fd_geom_1010_off + 377 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xx = cbuffer.data(fd_geom_1010_off + 378 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xy = cbuffer.data(fd_geom_1010_off + 379 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_xz = cbuffer.data(fd_geom_1010_off + 380 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_yy = cbuffer.data(fd_geom_1010_off + 381 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_yz = cbuffer.data(fd_geom_1010_off + 382 * ccomps * dcomps);

            auto g_z_0_x_0_xyy_zz = cbuffer.data(fd_geom_1010_off + 383 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xx = cbuffer.data(fd_geom_1010_off + 384 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xy = cbuffer.data(fd_geom_1010_off + 385 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_xz = cbuffer.data(fd_geom_1010_off + 386 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_yy = cbuffer.data(fd_geom_1010_off + 387 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_yz = cbuffer.data(fd_geom_1010_off + 388 * ccomps * dcomps);

            auto g_z_0_x_0_xyz_zz = cbuffer.data(fd_geom_1010_off + 389 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xx = cbuffer.data(fd_geom_1010_off + 390 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xy = cbuffer.data(fd_geom_1010_off + 391 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_xz = cbuffer.data(fd_geom_1010_off + 392 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_yy = cbuffer.data(fd_geom_1010_off + 393 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_yz = cbuffer.data(fd_geom_1010_off + 394 * ccomps * dcomps);

            auto g_z_0_x_0_xzz_zz = cbuffer.data(fd_geom_1010_off + 395 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xx = cbuffer.data(fd_geom_1010_off + 396 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xy = cbuffer.data(fd_geom_1010_off + 397 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_xz = cbuffer.data(fd_geom_1010_off + 398 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_yy = cbuffer.data(fd_geom_1010_off + 399 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_yz = cbuffer.data(fd_geom_1010_off + 400 * ccomps * dcomps);

            auto g_z_0_x_0_yyy_zz = cbuffer.data(fd_geom_1010_off + 401 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xx = cbuffer.data(fd_geom_1010_off + 402 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xy = cbuffer.data(fd_geom_1010_off + 403 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_xz = cbuffer.data(fd_geom_1010_off + 404 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_yy = cbuffer.data(fd_geom_1010_off + 405 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_yz = cbuffer.data(fd_geom_1010_off + 406 * ccomps * dcomps);

            auto g_z_0_x_0_yyz_zz = cbuffer.data(fd_geom_1010_off + 407 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xx = cbuffer.data(fd_geom_1010_off + 408 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xy = cbuffer.data(fd_geom_1010_off + 409 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_xz = cbuffer.data(fd_geom_1010_off + 410 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_yy = cbuffer.data(fd_geom_1010_off + 411 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_yz = cbuffer.data(fd_geom_1010_off + 412 * ccomps * dcomps);

            auto g_z_0_x_0_yzz_zz = cbuffer.data(fd_geom_1010_off + 413 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xx = cbuffer.data(fd_geom_1010_off + 414 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xy = cbuffer.data(fd_geom_1010_off + 415 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_xz = cbuffer.data(fd_geom_1010_off + 416 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_yy = cbuffer.data(fd_geom_1010_off + 417 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_yz = cbuffer.data(fd_geom_1010_off + 418 * ccomps * dcomps);

            auto g_z_0_x_0_zzz_zz = cbuffer.data(fd_geom_1010_off + 419 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xx = cbuffer.data(fd_geom_1010_off + 420 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xy = cbuffer.data(fd_geom_1010_off + 421 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_xz = cbuffer.data(fd_geom_1010_off + 422 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_yy = cbuffer.data(fd_geom_1010_off + 423 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_yz = cbuffer.data(fd_geom_1010_off + 424 * ccomps * dcomps);

            auto g_z_0_y_0_xxx_zz = cbuffer.data(fd_geom_1010_off + 425 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xx = cbuffer.data(fd_geom_1010_off + 426 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xy = cbuffer.data(fd_geom_1010_off + 427 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_xz = cbuffer.data(fd_geom_1010_off + 428 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_yy = cbuffer.data(fd_geom_1010_off + 429 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_yz = cbuffer.data(fd_geom_1010_off + 430 * ccomps * dcomps);

            auto g_z_0_y_0_xxy_zz = cbuffer.data(fd_geom_1010_off + 431 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xx = cbuffer.data(fd_geom_1010_off + 432 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xy = cbuffer.data(fd_geom_1010_off + 433 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_xz = cbuffer.data(fd_geom_1010_off + 434 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_yy = cbuffer.data(fd_geom_1010_off + 435 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_yz = cbuffer.data(fd_geom_1010_off + 436 * ccomps * dcomps);

            auto g_z_0_y_0_xxz_zz = cbuffer.data(fd_geom_1010_off + 437 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xx = cbuffer.data(fd_geom_1010_off + 438 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xy = cbuffer.data(fd_geom_1010_off + 439 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_xz = cbuffer.data(fd_geom_1010_off + 440 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_yy = cbuffer.data(fd_geom_1010_off + 441 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_yz = cbuffer.data(fd_geom_1010_off + 442 * ccomps * dcomps);

            auto g_z_0_y_0_xyy_zz = cbuffer.data(fd_geom_1010_off + 443 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xx = cbuffer.data(fd_geom_1010_off + 444 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xy = cbuffer.data(fd_geom_1010_off + 445 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_xz = cbuffer.data(fd_geom_1010_off + 446 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_yy = cbuffer.data(fd_geom_1010_off + 447 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_yz = cbuffer.data(fd_geom_1010_off + 448 * ccomps * dcomps);

            auto g_z_0_y_0_xyz_zz = cbuffer.data(fd_geom_1010_off + 449 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xx = cbuffer.data(fd_geom_1010_off + 450 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xy = cbuffer.data(fd_geom_1010_off + 451 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_xz = cbuffer.data(fd_geom_1010_off + 452 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_yy = cbuffer.data(fd_geom_1010_off + 453 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_yz = cbuffer.data(fd_geom_1010_off + 454 * ccomps * dcomps);

            auto g_z_0_y_0_xzz_zz = cbuffer.data(fd_geom_1010_off + 455 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xx = cbuffer.data(fd_geom_1010_off + 456 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xy = cbuffer.data(fd_geom_1010_off + 457 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_xz = cbuffer.data(fd_geom_1010_off + 458 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_yy = cbuffer.data(fd_geom_1010_off + 459 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_yz = cbuffer.data(fd_geom_1010_off + 460 * ccomps * dcomps);

            auto g_z_0_y_0_yyy_zz = cbuffer.data(fd_geom_1010_off + 461 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xx = cbuffer.data(fd_geom_1010_off + 462 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xy = cbuffer.data(fd_geom_1010_off + 463 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_xz = cbuffer.data(fd_geom_1010_off + 464 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_yy = cbuffer.data(fd_geom_1010_off + 465 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_yz = cbuffer.data(fd_geom_1010_off + 466 * ccomps * dcomps);

            auto g_z_0_y_0_yyz_zz = cbuffer.data(fd_geom_1010_off + 467 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xx = cbuffer.data(fd_geom_1010_off + 468 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xy = cbuffer.data(fd_geom_1010_off + 469 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_xz = cbuffer.data(fd_geom_1010_off + 470 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_yy = cbuffer.data(fd_geom_1010_off + 471 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_yz = cbuffer.data(fd_geom_1010_off + 472 * ccomps * dcomps);

            auto g_z_0_y_0_yzz_zz = cbuffer.data(fd_geom_1010_off + 473 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xx = cbuffer.data(fd_geom_1010_off + 474 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xy = cbuffer.data(fd_geom_1010_off + 475 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_xz = cbuffer.data(fd_geom_1010_off + 476 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_yy = cbuffer.data(fd_geom_1010_off + 477 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_yz = cbuffer.data(fd_geom_1010_off + 478 * ccomps * dcomps);

            auto g_z_0_y_0_zzz_zz = cbuffer.data(fd_geom_1010_off + 479 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xx = cbuffer.data(fd_geom_1010_off + 480 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xy = cbuffer.data(fd_geom_1010_off + 481 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_xz = cbuffer.data(fd_geom_1010_off + 482 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_yy = cbuffer.data(fd_geom_1010_off + 483 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_yz = cbuffer.data(fd_geom_1010_off + 484 * ccomps * dcomps);

            auto g_z_0_z_0_xxx_zz = cbuffer.data(fd_geom_1010_off + 485 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xx = cbuffer.data(fd_geom_1010_off + 486 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xy = cbuffer.data(fd_geom_1010_off + 487 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_xz = cbuffer.data(fd_geom_1010_off + 488 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_yy = cbuffer.data(fd_geom_1010_off + 489 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_yz = cbuffer.data(fd_geom_1010_off + 490 * ccomps * dcomps);

            auto g_z_0_z_0_xxy_zz = cbuffer.data(fd_geom_1010_off + 491 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xx = cbuffer.data(fd_geom_1010_off + 492 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xy = cbuffer.data(fd_geom_1010_off + 493 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_xz = cbuffer.data(fd_geom_1010_off + 494 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_yy = cbuffer.data(fd_geom_1010_off + 495 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_yz = cbuffer.data(fd_geom_1010_off + 496 * ccomps * dcomps);

            auto g_z_0_z_0_xxz_zz = cbuffer.data(fd_geom_1010_off + 497 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xx = cbuffer.data(fd_geom_1010_off + 498 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xy = cbuffer.data(fd_geom_1010_off + 499 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_xz = cbuffer.data(fd_geom_1010_off + 500 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_yy = cbuffer.data(fd_geom_1010_off + 501 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_yz = cbuffer.data(fd_geom_1010_off + 502 * ccomps * dcomps);

            auto g_z_0_z_0_xyy_zz = cbuffer.data(fd_geom_1010_off + 503 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xx = cbuffer.data(fd_geom_1010_off + 504 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xy = cbuffer.data(fd_geom_1010_off + 505 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_xz = cbuffer.data(fd_geom_1010_off + 506 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_yy = cbuffer.data(fd_geom_1010_off + 507 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_yz = cbuffer.data(fd_geom_1010_off + 508 * ccomps * dcomps);

            auto g_z_0_z_0_xyz_zz = cbuffer.data(fd_geom_1010_off + 509 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xx = cbuffer.data(fd_geom_1010_off + 510 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xy = cbuffer.data(fd_geom_1010_off + 511 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_xz = cbuffer.data(fd_geom_1010_off + 512 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_yy = cbuffer.data(fd_geom_1010_off + 513 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_yz = cbuffer.data(fd_geom_1010_off + 514 * ccomps * dcomps);

            auto g_z_0_z_0_xzz_zz = cbuffer.data(fd_geom_1010_off + 515 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xx = cbuffer.data(fd_geom_1010_off + 516 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xy = cbuffer.data(fd_geom_1010_off + 517 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_xz = cbuffer.data(fd_geom_1010_off + 518 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_yy = cbuffer.data(fd_geom_1010_off + 519 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_yz = cbuffer.data(fd_geom_1010_off + 520 * ccomps * dcomps);

            auto g_z_0_z_0_yyy_zz = cbuffer.data(fd_geom_1010_off + 521 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xx = cbuffer.data(fd_geom_1010_off + 522 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xy = cbuffer.data(fd_geom_1010_off + 523 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_xz = cbuffer.data(fd_geom_1010_off + 524 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_yy = cbuffer.data(fd_geom_1010_off + 525 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_yz = cbuffer.data(fd_geom_1010_off + 526 * ccomps * dcomps);

            auto g_z_0_z_0_yyz_zz = cbuffer.data(fd_geom_1010_off + 527 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xx = cbuffer.data(fd_geom_1010_off + 528 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xy = cbuffer.data(fd_geom_1010_off + 529 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_xz = cbuffer.data(fd_geom_1010_off + 530 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_yy = cbuffer.data(fd_geom_1010_off + 531 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_yz = cbuffer.data(fd_geom_1010_off + 532 * ccomps * dcomps);

            auto g_z_0_z_0_yzz_zz = cbuffer.data(fd_geom_1010_off + 533 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xx = cbuffer.data(fd_geom_1010_off + 534 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xy = cbuffer.data(fd_geom_1010_off + 535 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_xz = cbuffer.data(fd_geom_1010_off + 536 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_yy = cbuffer.data(fd_geom_1010_off + 537 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_yz = cbuffer.data(fd_geom_1010_off + 538 * ccomps * dcomps);

            auto g_z_0_z_0_zzz_zz = cbuffer.data(fd_geom_1010_off + 539 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_gpxx

            const auto gp_geom_1010_off = idx_geom_1010_gpxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxx_x = cbuffer.data(gp_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_xxxx_y = cbuffer.data(gp_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_xxxx_z = cbuffer.data(gp_geom_1010_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_xxx_x, g_0_0_x_0_xxx_y, g_0_0_x_0_xxx_z, g_x_0_x_0_xxx_x, g_x_0_x_0_xxx_xx, g_x_0_x_0_xxx_xy, g_x_0_x_0_xxx_xz, g_x_0_x_0_xxx_y, g_x_0_x_0_xxx_z, g_x_0_x_0_xxxx_x, g_x_0_x_0_xxxx_y, g_x_0_x_0_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxx_x[k] = -g_0_0_x_0_xxx_x[k] - g_x_0_x_0_xxx_x[k] * ab_x + g_x_0_x_0_xxx_xx[k];

                g_x_0_x_0_xxxx_y[k] = -g_0_0_x_0_xxx_y[k] - g_x_0_x_0_xxx_y[k] * ab_x + g_x_0_x_0_xxx_xy[k];

                g_x_0_x_0_xxxx_z[k] = -g_0_0_x_0_xxx_z[k] - g_x_0_x_0_xxx_z[k] * ab_x + g_x_0_x_0_xxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxy_x = cbuffer.data(gp_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_xxxy_y = cbuffer.data(gp_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_xxxy_z = cbuffer.data(gp_geom_1010_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxx_x, g_x_0_x_0_xxx_xy, g_x_0_x_0_xxx_y, g_x_0_x_0_xxx_yy, g_x_0_x_0_xxx_yz, g_x_0_x_0_xxx_z, g_x_0_x_0_xxxy_x, g_x_0_x_0_xxxy_y, g_x_0_x_0_xxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxy_x[k] = -g_x_0_x_0_xxx_x[k] * ab_y + g_x_0_x_0_xxx_xy[k];

                g_x_0_x_0_xxxy_y[k] = -g_x_0_x_0_xxx_y[k] * ab_y + g_x_0_x_0_xxx_yy[k];

                g_x_0_x_0_xxxy_z[k] = -g_x_0_x_0_xxx_z[k] * ab_y + g_x_0_x_0_xxx_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxxz_x = cbuffer.data(gp_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_xxxz_y = cbuffer.data(gp_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_xxxz_z = cbuffer.data(gp_geom_1010_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxx_x, g_x_0_x_0_xxx_xz, g_x_0_x_0_xxx_y, g_x_0_x_0_xxx_yz, g_x_0_x_0_xxx_z, g_x_0_x_0_xxx_zz, g_x_0_x_0_xxxz_x, g_x_0_x_0_xxxz_y, g_x_0_x_0_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxxz_x[k] = -g_x_0_x_0_xxx_x[k] * ab_z + g_x_0_x_0_xxx_xz[k];

                g_x_0_x_0_xxxz_y[k] = -g_x_0_x_0_xxx_y[k] * ab_z + g_x_0_x_0_xxx_yz[k];

                g_x_0_x_0_xxxz_z[k] = -g_x_0_x_0_xxx_z[k] * ab_z + g_x_0_x_0_xxx_zz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxyy_x = cbuffer.data(gp_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_xxyy_y = cbuffer.data(gp_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_xxyy_z = cbuffer.data(gp_geom_1010_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxy_x, g_x_0_x_0_xxy_xy, g_x_0_x_0_xxy_y, g_x_0_x_0_xxy_yy, g_x_0_x_0_xxy_yz, g_x_0_x_0_xxy_z, g_x_0_x_0_xxyy_x, g_x_0_x_0_xxyy_y, g_x_0_x_0_xxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxyy_x[k] = -g_x_0_x_0_xxy_x[k] * ab_y + g_x_0_x_0_xxy_xy[k];

                g_x_0_x_0_xxyy_y[k] = -g_x_0_x_0_xxy_y[k] * ab_y + g_x_0_x_0_xxy_yy[k];

                g_x_0_x_0_xxyy_z[k] = -g_x_0_x_0_xxy_z[k] * ab_y + g_x_0_x_0_xxy_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxyz_x = cbuffer.data(gp_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_xxyz_y = cbuffer.data(gp_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_xxyz_z = cbuffer.data(gp_geom_1010_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxyz_x, g_x_0_x_0_xxyz_y, g_x_0_x_0_xxyz_z, g_x_0_x_0_xxz_x, g_x_0_x_0_xxz_xy, g_x_0_x_0_xxz_y, g_x_0_x_0_xxz_yy, g_x_0_x_0_xxz_yz, g_x_0_x_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxyz_x[k] = -g_x_0_x_0_xxz_x[k] * ab_y + g_x_0_x_0_xxz_xy[k];

                g_x_0_x_0_xxyz_y[k] = -g_x_0_x_0_xxz_y[k] * ab_y + g_x_0_x_0_xxz_yy[k];

                g_x_0_x_0_xxyz_z[k] = -g_x_0_x_0_xxz_z[k] * ab_y + g_x_0_x_0_xxz_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xxzz_x = cbuffer.data(gp_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_x_0_xxzz_y = cbuffer.data(gp_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_x_0_xxzz_z = cbuffer.data(gp_geom_1010_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xxz_x, g_x_0_x_0_xxz_xz, g_x_0_x_0_xxz_y, g_x_0_x_0_xxz_yz, g_x_0_x_0_xxz_z, g_x_0_x_0_xxz_zz, g_x_0_x_0_xxzz_x, g_x_0_x_0_xxzz_y, g_x_0_x_0_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xxzz_x[k] = -g_x_0_x_0_xxz_x[k] * ab_z + g_x_0_x_0_xxz_xz[k];

                g_x_0_x_0_xxzz_y[k] = -g_x_0_x_0_xxz_y[k] * ab_z + g_x_0_x_0_xxz_yz[k];

                g_x_0_x_0_xxzz_z[k] = -g_x_0_x_0_xxz_z[k] * ab_z + g_x_0_x_0_xxz_zz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyyy_x = cbuffer.data(gp_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_x_0_xyyy_y = cbuffer.data(gp_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_x_0_xyyy_z = cbuffer.data(gp_geom_1010_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyy_x, g_x_0_x_0_xyy_xy, g_x_0_x_0_xyy_y, g_x_0_x_0_xyy_yy, g_x_0_x_0_xyy_yz, g_x_0_x_0_xyy_z, g_x_0_x_0_xyyy_x, g_x_0_x_0_xyyy_y, g_x_0_x_0_xyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyyy_x[k] = -g_x_0_x_0_xyy_x[k] * ab_y + g_x_0_x_0_xyy_xy[k];

                g_x_0_x_0_xyyy_y[k] = -g_x_0_x_0_xyy_y[k] * ab_y + g_x_0_x_0_xyy_yy[k];

                g_x_0_x_0_xyyy_z[k] = -g_x_0_x_0_xyy_z[k] * ab_y + g_x_0_x_0_xyy_yz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyyz_x = cbuffer.data(gp_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_x_0_xyyz_y = cbuffer.data(gp_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_x_0_xyyz_z = cbuffer.data(gp_geom_1010_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyyz_x, g_x_0_x_0_xyyz_y, g_x_0_x_0_xyyz_z, g_x_0_x_0_xyz_x, g_x_0_x_0_xyz_xy, g_x_0_x_0_xyz_y, g_x_0_x_0_xyz_yy, g_x_0_x_0_xyz_yz, g_x_0_x_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyyz_x[k] = -g_x_0_x_0_xyz_x[k] * ab_y + g_x_0_x_0_xyz_xy[k];

                g_x_0_x_0_xyyz_y[k] = -g_x_0_x_0_xyz_y[k] * ab_y + g_x_0_x_0_xyz_yy[k];

                g_x_0_x_0_xyyz_z[k] = -g_x_0_x_0_xyz_z[k] * ab_y + g_x_0_x_0_xyz_yz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xyzz_x = cbuffer.data(gp_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_x_0_xyzz_y = cbuffer.data(gp_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_x_0_xyzz_z = cbuffer.data(gp_geom_1010_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xyzz_x, g_x_0_x_0_xyzz_y, g_x_0_x_0_xyzz_z, g_x_0_x_0_xzz_x, g_x_0_x_0_xzz_xy, g_x_0_x_0_xzz_y, g_x_0_x_0_xzz_yy, g_x_0_x_0_xzz_yz, g_x_0_x_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xyzz_x[k] = -g_x_0_x_0_xzz_x[k] * ab_y + g_x_0_x_0_xzz_xy[k];

                g_x_0_x_0_xyzz_y[k] = -g_x_0_x_0_xzz_y[k] * ab_y + g_x_0_x_0_xzz_yy[k];

                g_x_0_x_0_xyzz_z[k] = -g_x_0_x_0_xzz_z[k] * ab_y + g_x_0_x_0_xzz_yz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xzzz_x = cbuffer.data(gp_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_x_0_xzzz_y = cbuffer.data(gp_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_x_0_xzzz_z = cbuffer.data(gp_geom_1010_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_xzz_x, g_x_0_x_0_xzz_xz, g_x_0_x_0_xzz_y, g_x_0_x_0_xzz_yz, g_x_0_x_0_xzz_z, g_x_0_x_0_xzz_zz, g_x_0_x_0_xzzz_x, g_x_0_x_0_xzzz_y, g_x_0_x_0_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xzzz_x[k] = -g_x_0_x_0_xzz_x[k] * ab_z + g_x_0_x_0_xzz_xz[k];

                g_x_0_x_0_xzzz_y[k] = -g_x_0_x_0_xzz_y[k] * ab_z + g_x_0_x_0_xzz_yz[k];

                g_x_0_x_0_xzzz_z[k] = -g_x_0_x_0_xzz_z[k] * ab_z + g_x_0_x_0_xzz_zz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyyy_x = cbuffer.data(gp_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_x_0_yyyy_y = cbuffer.data(gp_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_x_0_yyyy_z = cbuffer.data(gp_geom_1010_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyy_x, g_x_0_x_0_yyy_xy, g_x_0_x_0_yyy_y, g_x_0_x_0_yyy_yy, g_x_0_x_0_yyy_yz, g_x_0_x_0_yyy_z, g_x_0_x_0_yyyy_x, g_x_0_x_0_yyyy_y, g_x_0_x_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyyy_x[k] = -g_x_0_x_0_yyy_x[k] * ab_y + g_x_0_x_0_yyy_xy[k];

                g_x_0_x_0_yyyy_y[k] = -g_x_0_x_0_yyy_y[k] * ab_y + g_x_0_x_0_yyy_yy[k];

                g_x_0_x_0_yyyy_z[k] = -g_x_0_x_0_yyy_z[k] * ab_y + g_x_0_x_0_yyy_yz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyyz_x = cbuffer.data(gp_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_x_0_yyyz_y = cbuffer.data(gp_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_x_0_yyyz_z = cbuffer.data(gp_geom_1010_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyyz_x, g_x_0_x_0_yyyz_y, g_x_0_x_0_yyyz_z, g_x_0_x_0_yyz_x, g_x_0_x_0_yyz_xy, g_x_0_x_0_yyz_y, g_x_0_x_0_yyz_yy, g_x_0_x_0_yyz_yz, g_x_0_x_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyyz_x[k] = -g_x_0_x_0_yyz_x[k] * ab_y + g_x_0_x_0_yyz_xy[k];

                g_x_0_x_0_yyyz_y[k] = -g_x_0_x_0_yyz_y[k] * ab_y + g_x_0_x_0_yyz_yy[k];

                g_x_0_x_0_yyyz_z[k] = -g_x_0_x_0_yyz_z[k] * ab_y + g_x_0_x_0_yyz_yz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yyzz_x = cbuffer.data(gp_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_x_0_yyzz_y = cbuffer.data(gp_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_x_0_yyzz_z = cbuffer.data(gp_geom_1010_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yyzz_x, g_x_0_x_0_yyzz_y, g_x_0_x_0_yyzz_z, g_x_0_x_0_yzz_x, g_x_0_x_0_yzz_xy, g_x_0_x_0_yzz_y, g_x_0_x_0_yzz_yy, g_x_0_x_0_yzz_yz, g_x_0_x_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yyzz_x[k] = -g_x_0_x_0_yzz_x[k] * ab_y + g_x_0_x_0_yzz_xy[k];

                g_x_0_x_0_yyzz_y[k] = -g_x_0_x_0_yzz_y[k] * ab_y + g_x_0_x_0_yzz_yy[k];

                g_x_0_x_0_yyzz_z[k] = -g_x_0_x_0_yzz_z[k] * ab_y + g_x_0_x_0_yzz_yz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yzzz_x = cbuffer.data(gp_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_x_0_yzzz_y = cbuffer.data(gp_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_x_0_yzzz_z = cbuffer.data(gp_geom_1010_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yzzz_x, g_x_0_x_0_yzzz_y, g_x_0_x_0_yzzz_z, g_x_0_x_0_zzz_x, g_x_0_x_0_zzz_xy, g_x_0_x_0_zzz_y, g_x_0_x_0_zzz_yy, g_x_0_x_0_zzz_yz, g_x_0_x_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yzzz_x[k] = -g_x_0_x_0_zzz_x[k] * ab_y + g_x_0_x_0_zzz_xy[k];

                g_x_0_x_0_yzzz_y[k] = -g_x_0_x_0_zzz_y[k] * ab_y + g_x_0_x_0_zzz_yy[k];

                g_x_0_x_0_yzzz_z[k] = -g_x_0_x_0_zzz_z[k] * ab_y + g_x_0_x_0_zzz_yz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_zzzz_x = cbuffer.data(gp_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_x_0_zzzz_y = cbuffer.data(gp_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_x_0_zzzz_z = cbuffer.data(gp_geom_1010_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_zzz_x, g_x_0_x_0_zzz_xz, g_x_0_x_0_zzz_y, g_x_0_x_0_zzz_yz, g_x_0_x_0_zzz_z, g_x_0_x_0_zzz_zz, g_x_0_x_0_zzzz_x, g_x_0_x_0_zzzz_y, g_x_0_x_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_zzzz_x[k] = -g_x_0_x_0_zzz_x[k] * ab_z + g_x_0_x_0_zzz_xz[k];

                g_x_0_x_0_zzzz_y[k] = -g_x_0_x_0_zzz_y[k] * ab_z + g_x_0_x_0_zzz_yz[k];

                g_x_0_x_0_zzzz_z[k] = -g_x_0_x_0_zzz_z[k] * ab_z + g_x_0_x_0_zzz_zz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxx_x = cbuffer.data(gp_geom_1010_off + 45 * ccomps * dcomps);

            auto g_x_0_y_0_xxxx_y = cbuffer.data(gp_geom_1010_off + 46 * ccomps * dcomps);

            auto g_x_0_y_0_xxxx_z = cbuffer.data(gp_geom_1010_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_xxx_x, g_0_0_y_0_xxx_y, g_0_0_y_0_xxx_z, g_x_0_y_0_xxx_x, g_x_0_y_0_xxx_xx, g_x_0_y_0_xxx_xy, g_x_0_y_0_xxx_xz, g_x_0_y_0_xxx_y, g_x_0_y_0_xxx_z, g_x_0_y_0_xxxx_x, g_x_0_y_0_xxxx_y, g_x_0_y_0_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxx_x[k] = -g_0_0_y_0_xxx_x[k] - g_x_0_y_0_xxx_x[k] * ab_x + g_x_0_y_0_xxx_xx[k];

                g_x_0_y_0_xxxx_y[k] = -g_0_0_y_0_xxx_y[k] - g_x_0_y_0_xxx_y[k] * ab_x + g_x_0_y_0_xxx_xy[k];

                g_x_0_y_0_xxxx_z[k] = -g_0_0_y_0_xxx_z[k] - g_x_0_y_0_xxx_z[k] * ab_x + g_x_0_y_0_xxx_xz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxy_x = cbuffer.data(gp_geom_1010_off + 48 * ccomps * dcomps);

            auto g_x_0_y_0_xxxy_y = cbuffer.data(gp_geom_1010_off + 49 * ccomps * dcomps);

            auto g_x_0_y_0_xxxy_z = cbuffer.data(gp_geom_1010_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxx_x, g_x_0_y_0_xxx_xy, g_x_0_y_0_xxx_y, g_x_0_y_0_xxx_yy, g_x_0_y_0_xxx_yz, g_x_0_y_0_xxx_z, g_x_0_y_0_xxxy_x, g_x_0_y_0_xxxy_y, g_x_0_y_0_xxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxy_x[k] = -g_x_0_y_0_xxx_x[k] * ab_y + g_x_0_y_0_xxx_xy[k];

                g_x_0_y_0_xxxy_y[k] = -g_x_0_y_0_xxx_y[k] * ab_y + g_x_0_y_0_xxx_yy[k];

                g_x_0_y_0_xxxy_z[k] = -g_x_0_y_0_xxx_z[k] * ab_y + g_x_0_y_0_xxx_yz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxxz_x = cbuffer.data(gp_geom_1010_off + 51 * ccomps * dcomps);

            auto g_x_0_y_0_xxxz_y = cbuffer.data(gp_geom_1010_off + 52 * ccomps * dcomps);

            auto g_x_0_y_0_xxxz_z = cbuffer.data(gp_geom_1010_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxx_x, g_x_0_y_0_xxx_xz, g_x_0_y_0_xxx_y, g_x_0_y_0_xxx_yz, g_x_0_y_0_xxx_z, g_x_0_y_0_xxx_zz, g_x_0_y_0_xxxz_x, g_x_0_y_0_xxxz_y, g_x_0_y_0_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxxz_x[k] = -g_x_0_y_0_xxx_x[k] * ab_z + g_x_0_y_0_xxx_xz[k];

                g_x_0_y_0_xxxz_y[k] = -g_x_0_y_0_xxx_y[k] * ab_z + g_x_0_y_0_xxx_yz[k];

                g_x_0_y_0_xxxz_z[k] = -g_x_0_y_0_xxx_z[k] * ab_z + g_x_0_y_0_xxx_zz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxyy_x = cbuffer.data(gp_geom_1010_off + 54 * ccomps * dcomps);

            auto g_x_0_y_0_xxyy_y = cbuffer.data(gp_geom_1010_off + 55 * ccomps * dcomps);

            auto g_x_0_y_0_xxyy_z = cbuffer.data(gp_geom_1010_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxy_x, g_x_0_y_0_xxy_xy, g_x_0_y_0_xxy_y, g_x_0_y_0_xxy_yy, g_x_0_y_0_xxy_yz, g_x_0_y_0_xxy_z, g_x_0_y_0_xxyy_x, g_x_0_y_0_xxyy_y, g_x_0_y_0_xxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxyy_x[k] = -g_x_0_y_0_xxy_x[k] * ab_y + g_x_0_y_0_xxy_xy[k];

                g_x_0_y_0_xxyy_y[k] = -g_x_0_y_0_xxy_y[k] * ab_y + g_x_0_y_0_xxy_yy[k];

                g_x_0_y_0_xxyy_z[k] = -g_x_0_y_0_xxy_z[k] * ab_y + g_x_0_y_0_xxy_yz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxyz_x = cbuffer.data(gp_geom_1010_off + 57 * ccomps * dcomps);

            auto g_x_0_y_0_xxyz_y = cbuffer.data(gp_geom_1010_off + 58 * ccomps * dcomps);

            auto g_x_0_y_0_xxyz_z = cbuffer.data(gp_geom_1010_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxyz_x, g_x_0_y_0_xxyz_y, g_x_0_y_0_xxyz_z, g_x_0_y_0_xxz_x, g_x_0_y_0_xxz_xy, g_x_0_y_0_xxz_y, g_x_0_y_0_xxz_yy, g_x_0_y_0_xxz_yz, g_x_0_y_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxyz_x[k] = -g_x_0_y_0_xxz_x[k] * ab_y + g_x_0_y_0_xxz_xy[k];

                g_x_0_y_0_xxyz_y[k] = -g_x_0_y_0_xxz_y[k] * ab_y + g_x_0_y_0_xxz_yy[k];

                g_x_0_y_0_xxyz_z[k] = -g_x_0_y_0_xxz_z[k] * ab_y + g_x_0_y_0_xxz_yz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xxzz_x = cbuffer.data(gp_geom_1010_off + 60 * ccomps * dcomps);

            auto g_x_0_y_0_xxzz_y = cbuffer.data(gp_geom_1010_off + 61 * ccomps * dcomps);

            auto g_x_0_y_0_xxzz_z = cbuffer.data(gp_geom_1010_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xxz_x, g_x_0_y_0_xxz_xz, g_x_0_y_0_xxz_y, g_x_0_y_0_xxz_yz, g_x_0_y_0_xxz_z, g_x_0_y_0_xxz_zz, g_x_0_y_0_xxzz_x, g_x_0_y_0_xxzz_y, g_x_0_y_0_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xxzz_x[k] = -g_x_0_y_0_xxz_x[k] * ab_z + g_x_0_y_0_xxz_xz[k];

                g_x_0_y_0_xxzz_y[k] = -g_x_0_y_0_xxz_y[k] * ab_z + g_x_0_y_0_xxz_yz[k];

                g_x_0_y_0_xxzz_z[k] = -g_x_0_y_0_xxz_z[k] * ab_z + g_x_0_y_0_xxz_zz[k];
            }

            /// Set up 63-66 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyyy_x = cbuffer.data(gp_geom_1010_off + 63 * ccomps * dcomps);

            auto g_x_0_y_0_xyyy_y = cbuffer.data(gp_geom_1010_off + 64 * ccomps * dcomps);

            auto g_x_0_y_0_xyyy_z = cbuffer.data(gp_geom_1010_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyy_x, g_x_0_y_0_xyy_xy, g_x_0_y_0_xyy_y, g_x_0_y_0_xyy_yy, g_x_0_y_0_xyy_yz, g_x_0_y_0_xyy_z, g_x_0_y_0_xyyy_x, g_x_0_y_0_xyyy_y, g_x_0_y_0_xyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyyy_x[k] = -g_x_0_y_0_xyy_x[k] * ab_y + g_x_0_y_0_xyy_xy[k];

                g_x_0_y_0_xyyy_y[k] = -g_x_0_y_0_xyy_y[k] * ab_y + g_x_0_y_0_xyy_yy[k];

                g_x_0_y_0_xyyy_z[k] = -g_x_0_y_0_xyy_z[k] * ab_y + g_x_0_y_0_xyy_yz[k];
            }

            /// Set up 66-69 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyyz_x = cbuffer.data(gp_geom_1010_off + 66 * ccomps * dcomps);

            auto g_x_0_y_0_xyyz_y = cbuffer.data(gp_geom_1010_off + 67 * ccomps * dcomps);

            auto g_x_0_y_0_xyyz_z = cbuffer.data(gp_geom_1010_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyyz_x, g_x_0_y_0_xyyz_y, g_x_0_y_0_xyyz_z, g_x_0_y_0_xyz_x, g_x_0_y_0_xyz_xy, g_x_0_y_0_xyz_y, g_x_0_y_0_xyz_yy, g_x_0_y_0_xyz_yz, g_x_0_y_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyyz_x[k] = -g_x_0_y_0_xyz_x[k] * ab_y + g_x_0_y_0_xyz_xy[k];

                g_x_0_y_0_xyyz_y[k] = -g_x_0_y_0_xyz_y[k] * ab_y + g_x_0_y_0_xyz_yy[k];

                g_x_0_y_0_xyyz_z[k] = -g_x_0_y_0_xyz_z[k] * ab_y + g_x_0_y_0_xyz_yz[k];
            }

            /// Set up 69-72 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xyzz_x = cbuffer.data(gp_geom_1010_off + 69 * ccomps * dcomps);

            auto g_x_0_y_0_xyzz_y = cbuffer.data(gp_geom_1010_off + 70 * ccomps * dcomps);

            auto g_x_0_y_0_xyzz_z = cbuffer.data(gp_geom_1010_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xyzz_x, g_x_0_y_0_xyzz_y, g_x_0_y_0_xyzz_z, g_x_0_y_0_xzz_x, g_x_0_y_0_xzz_xy, g_x_0_y_0_xzz_y, g_x_0_y_0_xzz_yy, g_x_0_y_0_xzz_yz, g_x_0_y_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xyzz_x[k] = -g_x_0_y_0_xzz_x[k] * ab_y + g_x_0_y_0_xzz_xy[k];

                g_x_0_y_0_xyzz_y[k] = -g_x_0_y_0_xzz_y[k] * ab_y + g_x_0_y_0_xzz_yy[k];

                g_x_0_y_0_xyzz_z[k] = -g_x_0_y_0_xzz_z[k] * ab_y + g_x_0_y_0_xzz_yz[k];
            }

            /// Set up 72-75 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xzzz_x = cbuffer.data(gp_geom_1010_off + 72 * ccomps * dcomps);

            auto g_x_0_y_0_xzzz_y = cbuffer.data(gp_geom_1010_off + 73 * ccomps * dcomps);

            auto g_x_0_y_0_xzzz_z = cbuffer.data(gp_geom_1010_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_xzz_x, g_x_0_y_0_xzz_xz, g_x_0_y_0_xzz_y, g_x_0_y_0_xzz_yz, g_x_0_y_0_xzz_z, g_x_0_y_0_xzz_zz, g_x_0_y_0_xzzz_x, g_x_0_y_0_xzzz_y, g_x_0_y_0_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xzzz_x[k] = -g_x_0_y_0_xzz_x[k] * ab_z + g_x_0_y_0_xzz_xz[k];

                g_x_0_y_0_xzzz_y[k] = -g_x_0_y_0_xzz_y[k] * ab_z + g_x_0_y_0_xzz_yz[k];

                g_x_0_y_0_xzzz_z[k] = -g_x_0_y_0_xzz_z[k] * ab_z + g_x_0_y_0_xzz_zz[k];
            }

            /// Set up 75-78 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyyy_x = cbuffer.data(gp_geom_1010_off + 75 * ccomps * dcomps);

            auto g_x_0_y_0_yyyy_y = cbuffer.data(gp_geom_1010_off + 76 * ccomps * dcomps);

            auto g_x_0_y_0_yyyy_z = cbuffer.data(gp_geom_1010_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyy_x, g_x_0_y_0_yyy_xy, g_x_0_y_0_yyy_y, g_x_0_y_0_yyy_yy, g_x_0_y_0_yyy_yz, g_x_0_y_0_yyy_z, g_x_0_y_0_yyyy_x, g_x_0_y_0_yyyy_y, g_x_0_y_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyyy_x[k] = -g_x_0_y_0_yyy_x[k] * ab_y + g_x_0_y_0_yyy_xy[k];

                g_x_0_y_0_yyyy_y[k] = -g_x_0_y_0_yyy_y[k] * ab_y + g_x_0_y_0_yyy_yy[k];

                g_x_0_y_0_yyyy_z[k] = -g_x_0_y_0_yyy_z[k] * ab_y + g_x_0_y_0_yyy_yz[k];
            }

            /// Set up 78-81 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyyz_x = cbuffer.data(gp_geom_1010_off + 78 * ccomps * dcomps);

            auto g_x_0_y_0_yyyz_y = cbuffer.data(gp_geom_1010_off + 79 * ccomps * dcomps);

            auto g_x_0_y_0_yyyz_z = cbuffer.data(gp_geom_1010_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyyz_x, g_x_0_y_0_yyyz_y, g_x_0_y_0_yyyz_z, g_x_0_y_0_yyz_x, g_x_0_y_0_yyz_xy, g_x_0_y_0_yyz_y, g_x_0_y_0_yyz_yy, g_x_0_y_0_yyz_yz, g_x_0_y_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyyz_x[k] = -g_x_0_y_0_yyz_x[k] * ab_y + g_x_0_y_0_yyz_xy[k];

                g_x_0_y_0_yyyz_y[k] = -g_x_0_y_0_yyz_y[k] * ab_y + g_x_0_y_0_yyz_yy[k];

                g_x_0_y_0_yyyz_z[k] = -g_x_0_y_0_yyz_z[k] * ab_y + g_x_0_y_0_yyz_yz[k];
            }

            /// Set up 81-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yyzz_x = cbuffer.data(gp_geom_1010_off + 81 * ccomps * dcomps);

            auto g_x_0_y_0_yyzz_y = cbuffer.data(gp_geom_1010_off + 82 * ccomps * dcomps);

            auto g_x_0_y_0_yyzz_z = cbuffer.data(gp_geom_1010_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yyzz_x, g_x_0_y_0_yyzz_y, g_x_0_y_0_yyzz_z, g_x_0_y_0_yzz_x, g_x_0_y_0_yzz_xy, g_x_0_y_0_yzz_y, g_x_0_y_0_yzz_yy, g_x_0_y_0_yzz_yz, g_x_0_y_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yyzz_x[k] = -g_x_0_y_0_yzz_x[k] * ab_y + g_x_0_y_0_yzz_xy[k];

                g_x_0_y_0_yyzz_y[k] = -g_x_0_y_0_yzz_y[k] * ab_y + g_x_0_y_0_yzz_yy[k];

                g_x_0_y_0_yyzz_z[k] = -g_x_0_y_0_yzz_z[k] * ab_y + g_x_0_y_0_yzz_yz[k];
            }

            /// Set up 84-87 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yzzz_x = cbuffer.data(gp_geom_1010_off + 84 * ccomps * dcomps);

            auto g_x_0_y_0_yzzz_y = cbuffer.data(gp_geom_1010_off + 85 * ccomps * dcomps);

            auto g_x_0_y_0_yzzz_z = cbuffer.data(gp_geom_1010_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yzzz_x, g_x_0_y_0_yzzz_y, g_x_0_y_0_yzzz_z, g_x_0_y_0_zzz_x, g_x_0_y_0_zzz_xy, g_x_0_y_0_zzz_y, g_x_0_y_0_zzz_yy, g_x_0_y_0_zzz_yz, g_x_0_y_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yzzz_x[k] = -g_x_0_y_0_zzz_x[k] * ab_y + g_x_0_y_0_zzz_xy[k];

                g_x_0_y_0_yzzz_y[k] = -g_x_0_y_0_zzz_y[k] * ab_y + g_x_0_y_0_zzz_yy[k];

                g_x_0_y_0_yzzz_z[k] = -g_x_0_y_0_zzz_z[k] * ab_y + g_x_0_y_0_zzz_yz[k];
            }

            /// Set up 87-90 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_zzzz_x = cbuffer.data(gp_geom_1010_off + 87 * ccomps * dcomps);

            auto g_x_0_y_0_zzzz_y = cbuffer.data(gp_geom_1010_off + 88 * ccomps * dcomps);

            auto g_x_0_y_0_zzzz_z = cbuffer.data(gp_geom_1010_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_zzz_x, g_x_0_y_0_zzz_xz, g_x_0_y_0_zzz_y, g_x_0_y_0_zzz_yz, g_x_0_y_0_zzz_z, g_x_0_y_0_zzz_zz, g_x_0_y_0_zzzz_x, g_x_0_y_0_zzzz_y, g_x_0_y_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_zzzz_x[k] = -g_x_0_y_0_zzz_x[k] * ab_z + g_x_0_y_0_zzz_xz[k];

                g_x_0_y_0_zzzz_y[k] = -g_x_0_y_0_zzz_y[k] * ab_z + g_x_0_y_0_zzz_yz[k];

                g_x_0_y_0_zzzz_z[k] = -g_x_0_y_0_zzz_z[k] * ab_z + g_x_0_y_0_zzz_zz[k];
            }

            /// Set up 90-93 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxx_x = cbuffer.data(gp_geom_1010_off + 90 * ccomps * dcomps);

            auto g_x_0_z_0_xxxx_y = cbuffer.data(gp_geom_1010_off + 91 * ccomps * dcomps);

            auto g_x_0_z_0_xxxx_z = cbuffer.data(gp_geom_1010_off + 92 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_xxx_x, g_0_0_z_0_xxx_y, g_0_0_z_0_xxx_z, g_x_0_z_0_xxx_x, g_x_0_z_0_xxx_xx, g_x_0_z_0_xxx_xy, g_x_0_z_0_xxx_xz, g_x_0_z_0_xxx_y, g_x_0_z_0_xxx_z, g_x_0_z_0_xxxx_x, g_x_0_z_0_xxxx_y, g_x_0_z_0_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxx_x[k] = -g_0_0_z_0_xxx_x[k] - g_x_0_z_0_xxx_x[k] * ab_x + g_x_0_z_0_xxx_xx[k];

                g_x_0_z_0_xxxx_y[k] = -g_0_0_z_0_xxx_y[k] - g_x_0_z_0_xxx_y[k] * ab_x + g_x_0_z_0_xxx_xy[k];

                g_x_0_z_0_xxxx_z[k] = -g_0_0_z_0_xxx_z[k] - g_x_0_z_0_xxx_z[k] * ab_x + g_x_0_z_0_xxx_xz[k];
            }

            /// Set up 93-96 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxy_x = cbuffer.data(gp_geom_1010_off + 93 * ccomps * dcomps);

            auto g_x_0_z_0_xxxy_y = cbuffer.data(gp_geom_1010_off + 94 * ccomps * dcomps);

            auto g_x_0_z_0_xxxy_z = cbuffer.data(gp_geom_1010_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxx_x, g_x_0_z_0_xxx_xy, g_x_0_z_0_xxx_y, g_x_0_z_0_xxx_yy, g_x_0_z_0_xxx_yz, g_x_0_z_0_xxx_z, g_x_0_z_0_xxxy_x, g_x_0_z_0_xxxy_y, g_x_0_z_0_xxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxy_x[k] = -g_x_0_z_0_xxx_x[k] * ab_y + g_x_0_z_0_xxx_xy[k];

                g_x_0_z_0_xxxy_y[k] = -g_x_0_z_0_xxx_y[k] * ab_y + g_x_0_z_0_xxx_yy[k];

                g_x_0_z_0_xxxy_z[k] = -g_x_0_z_0_xxx_z[k] * ab_y + g_x_0_z_0_xxx_yz[k];
            }

            /// Set up 96-99 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxxz_x = cbuffer.data(gp_geom_1010_off + 96 * ccomps * dcomps);

            auto g_x_0_z_0_xxxz_y = cbuffer.data(gp_geom_1010_off + 97 * ccomps * dcomps);

            auto g_x_0_z_0_xxxz_z = cbuffer.data(gp_geom_1010_off + 98 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxx_x, g_x_0_z_0_xxx_xz, g_x_0_z_0_xxx_y, g_x_0_z_0_xxx_yz, g_x_0_z_0_xxx_z, g_x_0_z_0_xxx_zz, g_x_0_z_0_xxxz_x, g_x_0_z_0_xxxz_y, g_x_0_z_0_xxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxxz_x[k] = -g_x_0_z_0_xxx_x[k] * ab_z + g_x_0_z_0_xxx_xz[k];

                g_x_0_z_0_xxxz_y[k] = -g_x_0_z_0_xxx_y[k] * ab_z + g_x_0_z_0_xxx_yz[k];

                g_x_0_z_0_xxxz_z[k] = -g_x_0_z_0_xxx_z[k] * ab_z + g_x_0_z_0_xxx_zz[k];
            }

            /// Set up 99-102 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxyy_x = cbuffer.data(gp_geom_1010_off + 99 * ccomps * dcomps);

            auto g_x_0_z_0_xxyy_y = cbuffer.data(gp_geom_1010_off + 100 * ccomps * dcomps);

            auto g_x_0_z_0_xxyy_z = cbuffer.data(gp_geom_1010_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxy_x, g_x_0_z_0_xxy_xy, g_x_0_z_0_xxy_y, g_x_0_z_0_xxy_yy, g_x_0_z_0_xxy_yz, g_x_0_z_0_xxy_z, g_x_0_z_0_xxyy_x, g_x_0_z_0_xxyy_y, g_x_0_z_0_xxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxyy_x[k] = -g_x_0_z_0_xxy_x[k] * ab_y + g_x_0_z_0_xxy_xy[k];

                g_x_0_z_0_xxyy_y[k] = -g_x_0_z_0_xxy_y[k] * ab_y + g_x_0_z_0_xxy_yy[k];

                g_x_0_z_0_xxyy_z[k] = -g_x_0_z_0_xxy_z[k] * ab_y + g_x_0_z_0_xxy_yz[k];
            }

            /// Set up 102-105 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxyz_x = cbuffer.data(gp_geom_1010_off + 102 * ccomps * dcomps);

            auto g_x_0_z_0_xxyz_y = cbuffer.data(gp_geom_1010_off + 103 * ccomps * dcomps);

            auto g_x_0_z_0_xxyz_z = cbuffer.data(gp_geom_1010_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxyz_x, g_x_0_z_0_xxyz_y, g_x_0_z_0_xxyz_z, g_x_0_z_0_xxz_x, g_x_0_z_0_xxz_xy, g_x_0_z_0_xxz_y, g_x_0_z_0_xxz_yy, g_x_0_z_0_xxz_yz, g_x_0_z_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxyz_x[k] = -g_x_0_z_0_xxz_x[k] * ab_y + g_x_0_z_0_xxz_xy[k];

                g_x_0_z_0_xxyz_y[k] = -g_x_0_z_0_xxz_y[k] * ab_y + g_x_0_z_0_xxz_yy[k];

                g_x_0_z_0_xxyz_z[k] = -g_x_0_z_0_xxz_z[k] * ab_y + g_x_0_z_0_xxz_yz[k];
            }

            /// Set up 105-108 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xxzz_x = cbuffer.data(gp_geom_1010_off + 105 * ccomps * dcomps);

            auto g_x_0_z_0_xxzz_y = cbuffer.data(gp_geom_1010_off + 106 * ccomps * dcomps);

            auto g_x_0_z_0_xxzz_z = cbuffer.data(gp_geom_1010_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xxz_x, g_x_0_z_0_xxz_xz, g_x_0_z_0_xxz_y, g_x_0_z_0_xxz_yz, g_x_0_z_0_xxz_z, g_x_0_z_0_xxz_zz, g_x_0_z_0_xxzz_x, g_x_0_z_0_xxzz_y, g_x_0_z_0_xxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xxzz_x[k] = -g_x_0_z_0_xxz_x[k] * ab_z + g_x_0_z_0_xxz_xz[k];

                g_x_0_z_0_xxzz_y[k] = -g_x_0_z_0_xxz_y[k] * ab_z + g_x_0_z_0_xxz_yz[k];

                g_x_0_z_0_xxzz_z[k] = -g_x_0_z_0_xxz_z[k] * ab_z + g_x_0_z_0_xxz_zz[k];
            }

            /// Set up 108-111 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyyy_x = cbuffer.data(gp_geom_1010_off + 108 * ccomps * dcomps);

            auto g_x_0_z_0_xyyy_y = cbuffer.data(gp_geom_1010_off + 109 * ccomps * dcomps);

            auto g_x_0_z_0_xyyy_z = cbuffer.data(gp_geom_1010_off + 110 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyy_x, g_x_0_z_0_xyy_xy, g_x_0_z_0_xyy_y, g_x_0_z_0_xyy_yy, g_x_0_z_0_xyy_yz, g_x_0_z_0_xyy_z, g_x_0_z_0_xyyy_x, g_x_0_z_0_xyyy_y, g_x_0_z_0_xyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyyy_x[k] = -g_x_0_z_0_xyy_x[k] * ab_y + g_x_0_z_0_xyy_xy[k];

                g_x_0_z_0_xyyy_y[k] = -g_x_0_z_0_xyy_y[k] * ab_y + g_x_0_z_0_xyy_yy[k];

                g_x_0_z_0_xyyy_z[k] = -g_x_0_z_0_xyy_z[k] * ab_y + g_x_0_z_0_xyy_yz[k];
            }

            /// Set up 111-114 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyyz_x = cbuffer.data(gp_geom_1010_off + 111 * ccomps * dcomps);

            auto g_x_0_z_0_xyyz_y = cbuffer.data(gp_geom_1010_off + 112 * ccomps * dcomps);

            auto g_x_0_z_0_xyyz_z = cbuffer.data(gp_geom_1010_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyyz_x, g_x_0_z_0_xyyz_y, g_x_0_z_0_xyyz_z, g_x_0_z_0_xyz_x, g_x_0_z_0_xyz_xy, g_x_0_z_0_xyz_y, g_x_0_z_0_xyz_yy, g_x_0_z_0_xyz_yz, g_x_0_z_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyyz_x[k] = -g_x_0_z_0_xyz_x[k] * ab_y + g_x_0_z_0_xyz_xy[k];

                g_x_0_z_0_xyyz_y[k] = -g_x_0_z_0_xyz_y[k] * ab_y + g_x_0_z_0_xyz_yy[k];

                g_x_0_z_0_xyyz_z[k] = -g_x_0_z_0_xyz_z[k] * ab_y + g_x_0_z_0_xyz_yz[k];
            }

            /// Set up 114-117 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xyzz_x = cbuffer.data(gp_geom_1010_off + 114 * ccomps * dcomps);

            auto g_x_0_z_0_xyzz_y = cbuffer.data(gp_geom_1010_off + 115 * ccomps * dcomps);

            auto g_x_0_z_0_xyzz_z = cbuffer.data(gp_geom_1010_off + 116 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xyzz_x, g_x_0_z_0_xyzz_y, g_x_0_z_0_xyzz_z, g_x_0_z_0_xzz_x, g_x_0_z_0_xzz_xy, g_x_0_z_0_xzz_y, g_x_0_z_0_xzz_yy, g_x_0_z_0_xzz_yz, g_x_0_z_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xyzz_x[k] = -g_x_0_z_0_xzz_x[k] * ab_y + g_x_0_z_0_xzz_xy[k];

                g_x_0_z_0_xyzz_y[k] = -g_x_0_z_0_xzz_y[k] * ab_y + g_x_0_z_0_xzz_yy[k];

                g_x_0_z_0_xyzz_z[k] = -g_x_0_z_0_xzz_z[k] * ab_y + g_x_0_z_0_xzz_yz[k];
            }

            /// Set up 117-120 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xzzz_x = cbuffer.data(gp_geom_1010_off + 117 * ccomps * dcomps);

            auto g_x_0_z_0_xzzz_y = cbuffer.data(gp_geom_1010_off + 118 * ccomps * dcomps);

            auto g_x_0_z_0_xzzz_z = cbuffer.data(gp_geom_1010_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_xzz_x, g_x_0_z_0_xzz_xz, g_x_0_z_0_xzz_y, g_x_0_z_0_xzz_yz, g_x_0_z_0_xzz_z, g_x_0_z_0_xzz_zz, g_x_0_z_0_xzzz_x, g_x_0_z_0_xzzz_y, g_x_0_z_0_xzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xzzz_x[k] = -g_x_0_z_0_xzz_x[k] * ab_z + g_x_0_z_0_xzz_xz[k];

                g_x_0_z_0_xzzz_y[k] = -g_x_0_z_0_xzz_y[k] * ab_z + g_x_0_z_0_xzz_yz[k];

                g_x_0_z_0_xzzz_z[k] = -g_x_0_z_0_xzz_z[k] * ab_z + g_x_0_z_0_xzz_zz[k];
            }

            /// Set up 120-123 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyyy_x = cbuffer.data(gp_geom_1010_off + 120 * ccomps * dcomps);

            auto g_x_0_z_0_yyyy_y = cbuffer.data(gp_geom_1010_off + 121 * ccomps * dcomps);

            auto g_x_0_z_0_yyyy_z = cbuffer.data(gp_geom_1010_off + 122 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyy_x, g_x_0_z_0_yyy_xy, g_x_0_z_0_yyy_y, g_x_0_z_0_yyy_yy, g_x_0_z_0_yyy_yz, g_x_0_z_0_yyy_z, g_x_0_z_0_yyyy_x, g_x_0_z_0_yyyy_y, g_x_0_z_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyyy_x[k] = -g_x_0_z_0_yyy_x[k] * ab_y + g_x_0_z_0_yyy_xy[k];

                g_x_0_z_0_yyyy_y[k] = -g_x_0_z_0_yyy_y[k] * ab_y + g_x_0_z_0_yyy_yy[k];

                g_x_0_z_0_yyyy_z[k] = -g_x_0_z_0_yyy_z[k] * ab_y + g_x_0_z_0_yyy_yz[k];
            }

            /// Set up 123-126 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyyz_x = cbuffer.data(gp_geom_1010_off + 123 * ccomps * dcomps);

            auto g_x_0_z_0_yyyz_y = cbuffer.data(gp_geom_1010_off + 124 * ccomps * dcomps);

            auto g_x_0_z_0_yyyz_z = cbuffer.data(gp_geom_1010_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyyz_x, g_x_0_z_0_yyyz_y, g_x_0_z_0_yyyz_z, g_x_0_z_0_yyz_x, g_x_0_z_0_yyz_xy, g_x_0_z_0_yyz_y, g_x_0_z_0_yyz_yy, g_x_0_z_0_yyz_yz, g_x_0_z_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyyz_x[k] = -g_x_0_z_0_yyz_x[k] * ab_y + g_x_0_z_0_yyz_xy[k];

                g_x_0_z_0_yyyz_y[k] = -g_x_0_z_0_yyz_y[k] * ab_y + g_x_0_z_0_yyz_yy[k];

                g_x_0_z_0_yyyz_z[k] = -g_x_0_z_0_yyz_z[k] * ab_y + g_x_0_z_0_yyz_yz[k];
            }

            /// Set up 126-129 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yyzz_x = cbuffer.data(gp_geom_1010_off + 126 * ccomps * dcomps);

            auto g_x_0_z_0_yyzz_y = cbuffer.data(gp_geom_1010_off + 127 * ccomps * dcomps);

            auto g_x_0_z_0_yyzz_z = cbuffer.data(gp_geom_1010_off + 128 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yyzz_x, g_x_0_z_0_yyzz_y, g_x_0_z_0_yyzz_z, g_x_0_z_0_yzz_x, g_x_0_z_0_yzz_xy, g_x_0_z_0_yzz_y, g_x_0_z_0_yzz_yy, g_x_0_z_0_yzz_yz, g_x_0_z_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yyzz_x[k] = -g_x_0_z_0_yzz_x[k] * ab_y + g_x_0_z_0_yzz_xy[k];

                g_x_0_z_0_yyzz_y[k] = -g_x_0_z_0_yzz_y[k] * ab_y + g_x_0_z_0_yzz_yy[k];

                g_x_0_z_0_yyzz_z[k] = -g_x_0_z_0_yzz_z[k] * ab_y + g_x_0_z_0_yzz_yz[k];
            }

            /// Set up 129-132 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yzzz_x = cbuffer.data(gp_geom_1010_off + 129 * ccomps * dcomps);

            auto g_x_0_z_0_yzzz_y = cbuffer.data(gp_geom_1010_off + 130 * ccomps * dcomps);

            auto g_x_0_z_0_yzzz_z = cbuffer.data(gp_geom_1010_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yzzz_x, g_x_0_z_0_yzzz_y, g_x_0_z_0_yzzz_z, g_x_0_z_0_zzz_x, g_x_0_z_0_zzz_xy, g_x_0_z_0_zzz_y, g_x_0_z_0_zzz_yy, g_x_0_z_0_zzz_yz, g_x_0_z_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yzzz_x[k] = -g_x_0_z_0_zzz_x[k] * ab_y + g_x_0_z_0_zzz_xy[k];

                g_x_0_z_0_yzzz_y[k] = -g_x_0_z_0_zzz_y[k] * ab_y + g_x_0_z_0_zzz_yy[k];

                g_x_0_z_0_yzzz_z[k] = -g_x_0_z_0_zzz_z[k] * ab_y + g_x_0_z_0_zzz_yz[k];
            }

            /// Set up 132-135 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_zzzz_x = cbuffer.data(gp_geom_1010_off + 132 * ccomps * dcomps);

            auto g_x_0_z_0_zzzz_y = cbuffer.data(gp_geom_1010_off + 133 * ccomps * dcomps);

            auto g_x_0_z_0_zzzz_z = cbuffer.data(gp_geom_1010_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_zzz_x, g_x_0_z_0_zzz_xz, g_x_0_z_0_zzz_y, g_x_0_z_0_zzz_yz, g_x_0_z_0_zzz_z, g_x_0_z_0_zzz_zz, g_x_0_z_0_zzzz_x, g_x_0_z_0_zzzz_y, g_x_0_z_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_zzzz_x[k] = -g_x_0_z_0_zzz_x[k] * ab_z + g_x_0_z_0_zzz_xz[k];

                g_x_0_z_0_zzzz_y[k] = -g_x_0_z_0_zzz_y[k] * ab_z + g_x_0_z_0_zzz_yz[k];

                g_x_0_z_0_zzzz_z[k] = -g_x_0_z_0_zzz_z[k] * ab_z + g_x_0_z_0_zzz_zz[k];
            }

            /// Set up 135-138 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxx_x = cbuffer.data(gp_geom_1010_off + 135 * ccomps * dcomps);

            auto g_y_0_x_0_xxxx_y = cbuffer.data(gp_geom_1010_off + 136 * ccomps * dcomps);

            auto g_y_0_x_0_xxxx_z = cbuffer.data(gp_geom_1010_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxx_x, g_y_0_x_0_xxx_xx, g_y_0_x_0_xxx_xy, g_y_0_x_0_xxx_xz, g_y_0_x_0_xxx_y, g_y_0_x_0_xxx_z, g_y_0_x_0_xxxx_x, g_y_0_x_0_xxxx_y, g_y_0_x_0_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxx_x[k] = -g_y_0_x_0_xxx_x[k] * ab_x + g_y_0_x_0_xxx_xx[k];

                g_y_0_x_0_xxxx_y[k] = -g_y_0_x_0_xxx_y[k] * ab_x + g_y_0_x_0_xxx_xy[k];

                g_y_0_x_0_xxxx_z[k] = -g_y_0_x_0_xxx_z[k] * ab_x + g_y_0_x_0_xxx_xz[k];
            }

            /// Set up 138-141 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxy_x = cbuffer.data(gp_geom_1010_off + 138 * ccomps * dcomps);

            auto g_y_0_x_0_xxxy_y = cbuffer.data(gp_geom_1010_off + 139 * ccomps * dcomps);

            auto g_y_0_x_0_xxxy_z = cbuffer.data(gp_geom_1010_off + 140 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxy_x, g_y_0_x_0_xxxy_y, g_y_0_x_0_xxxy_z, g_y_0_x_0_xxy_x, g_y_0_x_0_xxy_xx, g_y_0_x_0_xxy_xy, g_y_0_x_0_xxy_xz, g_y_0_x_0_xxy_y, g_y_0_x_0_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxy_x[k] = -g_y_0_x_0_xxy_x[k] * ab_x + g_y_0_x_0_xxy_xx[k];

                g_y_0_x_0_xxxy_y[k] = -g_y_0_x_0_xxy_y[k] * ab_x + g_y_0_x_0_xxy_xy[k];

                g_y_0_x_0_xxxy_z[k] = -g_y_0_x_0_xxy_z[k] * ab_x + g_y_0_x_0_xxy_xz[k];
            }

            /// Set up 141-144 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxxz_x = cbuffer.data(gp_geom_1010_off + 141 * ccomps * dcomps);

            auto g_y_0_x_0_xxxz_y = cbuffer.data(gp_geom_1010_off + 142 * ccomps * dcomps);

            auto g_y_0_x_0_xxxz_z = cbuffer.data(gp_geom_1010_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxxz_x, g_y_0_x_0_xxxz_y, g_y_0_x_0_xxxz_z, g_y_0_x_0_xxz_x, g_y_0_x_0_xxz_xx, g_y_0_x_0_xxz_xy, g_y_0_x_0_xxz_xz, g_y_0_x_0_xxz_y, g_y_0_x_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxxz_x[k] = -g_y_0_x_0_xxz_x[k] * ab_x + g_y_0_x_0_xxz_xx[k];

                g_y_0_x_0_xxxz_y[k] = -g_y_0_x_0_xxz_y[k] * ab_x + g_y_0_x_0_xxz_xy[k];

                g_y_0_x_0_xxxz_z[k] = -g_y_0_x_0_xxz_z[k] * ab_x + g_y_0_x_0_xxz_xz[k];
            }

            /// Set up 144-147 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxyy_x = cbuffer.data(gp_geom_1010_off + 144 * ccomps * dcomps);

            auto g_y_0_x_0_xxyy_y = cbuffer.data(gp_geom_1010_off + 145 * ccomps * dcomps);

            auto g_y_0_x_0_xxyy_z = cbuffer.data(gp_geom_1010_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxyy_x, g_y_0_x_0_xxyy_y, g_y_0_x_0_xxyy_z, g_y_0_x_0_xyy_x, g_y_0_x_0_xyy_xx, g_y_0_x_0_xyy_xy, g_y_0_x_0_xyy_xz, g_y_0_x_0_xyy_y, g_y_0_x_0_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxyy_x[k] = -g_y_0_x_0_xyy_x[k] * ab_x + g_y_0_x_0_xyy_xx[k];

                g_y_0_x_0_xxyy_y[k] = -g_y_0_x_0_xyy_y[k] * ab_x + g_y_0_x_0_xyy_xy[k];

                g_y_0_x_0_xxyy_z[k] = -g_y_0_x_0_xyy_z[k] * ab_x + g_y_0_x_0_xyy_xz[k];
            }

            /// Set up 147-150 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxyz_x = cbuffer.data(gp_geom_1010_off + 147 * ccomps * dcomps);

            auto g_y_0_x_0_xxyz_y = cbuffer.data(gp_geom_1010_off + 148 * ccomps * dcomps);

            auto g_y_0_x_0_xxyz_z = cbuffer.data(gp_geom_1010_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxyz_x, g_y_0_x_0_xxyz_y, g_y_0_x_0_xxyz_z, g_y_0_x_0_xyz_x, g_y_0_x_0_xyz_xx, g_y_0_x_0_xyz_xy, g_y_0_x_0_xyz_xz, g_y_0_x_0_xyz_y, g_y_0_x_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxyz_x[k] = -g_y_0_x_0_xyz_x[k] * ab_x + g_y_0_x_0_xyz_xx[k];

                g_y_0_x_0_xxyz_y[k] = -g_y_0_x_0_xyz_y[k] * ab_x + g_y_0_x_0_xyz_xy[k];

                g_y_0_x_0_xxyz_z[k] = -g_y_0_x_0_xyz_z[k] * ab_x + g_y_0_x_0_xyz_xz[k];
            }

            /// Set up 150-153 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xxzz_x = cbuffer.data(gp_geom_1010_off + 150 * ccomps * dcomps);

            auto g_y_0_x_0_xxzz_y = cbuffer.data(gp_geom_1010_off + 151 * ccomps * dcomps);

            auto g_y_0_x_0_xxzz_z = cbuffer.data(gp_geom_1010_off + 152 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xxzz_x, g_y_0_x_0_xxzz_y, g_y_0_x_0_xxzz_z, g_y_0_x_0_xzz_x, g_y_0_x_0_xzz_xx, g_y_0_x_0_xzz_xy, g_y_0_x_0_xzz_xz, g_y_0_x_0_xzz_y, g_y_0_x_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xxzz_x[k] = -g_y_0_x_0_xzz_x[k] * ab_x + g_y_0_x_0_xzz_xx[k];

                g_y_0_x_0_xxzz_y[k] = -g_y_0_x_0_xzz_y[k] * ab_x + g_y_0_x_0_xzz_xy[k];

                g_y_0_x_0_xxzz_z[k] = -g_y_0_x_0_xzz_z[k] * ab_x + g_y_0_x_0_xzz_xz[k];
            }

            /// Set up 153-156 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyyy_x = cbuffer.data(gp_geom_1010_off + 153 * ccomps * dcomps);

            auto g_y_0_x_0_xyyy_y = cbuffer.data(gp_geom_1010_off + 154 * ccomps * dcomps);

            auto g_y_0_x_0_xyyy_z = cbuffer.data(gp_geom_1010_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyyy_x, g_y_0_x_0_xyyy_y, g_y_0_x_0_xyyy_z, g_y_0_x_0_yyy_x, g_y_0_x_0_yyy_xx, g_y_0_x_0_yyy_xy, g_y_0_x_0_yyy_xz, g_y_0_x_0_yyy_y, g_y_0_x_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyyy_x[k] = -g_y_0_x_0_yyy_x[k] * ab_x + g_y_0_x_0_yyy_xx[k];

                g_y_0_x_0_xyyy_y[k] = -g_y_0_x_0_yyy_y[k] * ab_x + g_y_0_x_0_yyy_xy[k];

                g_y_0_x_0_xyyy_z[k] = -g_y_0_x_0_yyy_z[k] * ab_x + g_y_0_x_0_yyy_xz[k];
            }

            /// Set up 156-159 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyyz_x = cbuffer.data(gp_geom_1010_off + 156 * ccomps * dcomps);

            auto g_y_0_x_0_xyyz_y = cbuffer.data(gp_geom_1010_off + 157 * ccomps * dcomps);

            auto g_y_0_x_0_xyyz_z = cbuffer.data(gp_geom_1010_off + 158 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyyz_x, g_y_0_x_0_xyyz_y, g_y_0_x_0_xyyz_z, g_y_0_x_0_yyz_x, g_y_0_x_0_yyz_xx, g_y_0_x_0_yyz_xy, g_y_0_x_0_yyz_xz, g_y_0_x_0_yyz_y, g_y_0_x_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyyz_x[k] = -g_y_0_x_0_yyz_x[k] * ab_x + g_y_0_x_0_yyz_xx[k];

                g_y_0_x_0_xyyz_y[k] = -g_y_0_x_0_yyz_y[k] * ab_x + g_y_0_x_0_yyz_xy[k];

                g_y_0_x_0_xyyz_z[k] = -g_y_0_x_0_yyz_z[k] * ab_x + g_y_0_x_0_yyz_xz[k];
            }

            /// Set up 159-162 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xyzz_x = cbuffer.data(gp_geom_1010_off + 159 * ccomps * dcomps);

            auto g_y_0_x_0_xyzz_y = cbuffer.data(gp_geom_1010_off + 160 * ccomps * dcomps);

            auto g_y_0_x_0_xyzz_z = cbuffer.data(gp_geom_1010_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xyzz_x, g_y_0_x_0_xyzz_y, g_y_0_x_0_xyzz_z, g_y_0_x_0_yzz_x, g_y_0_x_0_yzz_xx, g_y_0_x_0_yzz_xy, g_y_0_x_0_yzz_xz, g_y_0_x_0_yzz_y, g_y_0_x_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xyzz_x[k] = -g_y_0_x_0_yzz_x[k] * ab_x + g_y_0_x_0_yzz_xx[k];

                g_y_0_x_0_xyzz_y[k] = -g_y_0_x_0_yzz_y[k] * ab_x + g_y_0_x_0_yzz_xy[k];

                g_y_0_x_0_xyzz_z[k] = -g_y_0_x_0_yzz_z[k] * ab_x + g_y_0_x_0_yzz_xz[k];
            }

            /// Set up 162-165 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xzzz_x = cbuffer.data(gp_geom_1010_off + 162 * ccomps * dcomps);

            auto g_y_0_x_0_xzzz_y = cbuffer.data(gp_geom_1010_off + 163 * ccomps * dcomps);

            auto g_y_0_x_0_xzzz_z = cbuffer.data(gp_geom_1010_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xzzz_x, g_y_0_x_0_xzzz_y, g_y_0_x_0_xzzz_z, g_y_0_x_0_zzz_x, g_y_0_x_0_zzz_xx, g_y_0_x_0_zzz_xy, g_y_0_x_0_zzz_xz, g_y_0_x_0_zzz_y, g_y_0_x_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xzzz_x[k] = -g_y_0_x_0_zzz_x[k] * ab_x + g_y_0_x_0_zzz_xx[k];

                g_y_0_x_0_xzzz_y[k] = -g_y_0_x_0_zzz_y[k] * ab_x + g_y_0_x_0_zzz_xy[k];

                g_y_0_x_0_xzzz_z[k] = -g_y_0_x_0_zzz_z[k] * ab_x + g_y_0_x_0_zzz_xz[k];
            }

            /// Set up 165-168 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyyy_x = cbuffer.data(gp_geom_1010_off + 165 * ccomps * dcomps);

            auto g_y_0_x_0_yyyy_y = cbuffer.data(gp_geom_1010_off + 166 * ccomps * dcomps);

            auto g_y_0_x_0_yyyy_z = cbuffer.data(gp_geom_1010_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_yyy_x, g_0_0_x_0_yyy_y, g_0_0_x_0_yyy_z, g_y_0_x_0_yyy_x, g_y_0_x_0_yyy_xy, g_y_0_x_0_yyy_y, g_y_0_x_0_yyy_yy, g_y_0_x_0_yyy_yz, g_y_0_x_0_yyy_z, g_y_0_x_0_yyyy_x, g_y_0_x_0_yyyy_y, g_y_0_x_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyyy_x[k] = -g_0_0_x_0_yyy_x[k] - g_y_0_x_0_yyy_x[k] * ab_y + g_y_0_x_0_yyy_xy[k];

                g_y_0_x_0_yyyy_y[k] = -g_0_0_x_0_yyy_y[k] - g_y_0_x_0_yyy_y[k] * ab_y + g_y_0_x_0_yyy_yy[k];

                g_y_0_x_0_yyyy_z[k] = -g_0_0_x_0_yyy_z[k] - g_y_0_x_0_yyy_z[k] * ab_y + g_y_0_x_0_yyy_yz[k];
            }

            /// Set up 168-171 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyyz_x = cbuffer.data(gp_geom_1010_off + 168 * ccomps * dcomps);

            auto g_y_0_x_0_yyyz_y = cbuffer.data(gp_geom_1010_off + 169 * ccomps * dcomps);

            auto g_y_0_x_0_yyyz_z = cbuffer.data(gp_geom_1010_off + 170 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yyy_x, g_y_0_x_0_yyy_xz, g_y_0_x_0_yyy_y, g_y_0_x_0_yyy_yz, g_y_0_x_0_yyy_z, g_y_0_x_0_yyy_zz, g_y_0_x_0_yyyz_x, g_y_0_x_0_yyyz_y, g_y_0_x_0_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyyz_x[k] = -g_y_0_x_0_yyy_x[k] * ab_z + g_y_0_x_0_yyy_xz[k];

                g_y_0_x_0_yyyz_y[k] = -g_y_0_x_0_yyy_y[k] * ab_z + g_y_0_x_0_yyy_yz[k];

                g_y_0_x_0_yyyz_z[k] = -g_y_0_x_0_yyy_z[k] * ab_z + g_y_0_x_0_yyy_zz[k];
            }

            /// Set up 171-174 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yyzz_x = cbuffer.data(gp_geom_1010_off + 171 * ccomps * dcomps);

            auto g_y_0_x_0_yyzz_y = cbuffer.data(gp_geom_1010_off + 172 * ccomps * dcomps);

            auto g_y_0_x_0_yyzz_z = cbuffer.data(gp_geom_1010_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yyz_x, g_y_0_x_0_yyz_xz, g_y_0_x_0_yyz_y, g_y_0_x_0_yyz_yz, g_y_0_x_0_yyz_z, g_y_0_x_0_yyz_zz, g_y_0_x_0_yyzz_x, g_y_0_x_0_yyzz_y, g_y_0_x_0_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yyzz_x[k] = -g_y_0_x_0_yyz_x[k] * ab_z + g_y_0_x_0_yyz_xz[k];

                g_y_0_x_0_yyzz_y[k] = -g_y_0_x_0_yyz_y[k] * ab_z + g_y_0_x_0_yyz_yz[k];

                g_y_0_x_0_yyzz_z[k] = -g_y_0_x_0_yyz_z[k] * ab_z + g_y_0_x_0_yyz_zz[k];
            }

            /// Set up 174-177 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yzzz_x = cbuffer.data(gp_geom_1010_off + 174 * ccomps * dcomps);

            auto g_y_0_x_0_yzzz_y = cbuffer.data(gp_geom_1010_off + 175 * ccomps * dcomps);

            auto g_y_0_x_0_yzzz_z = cbuffer.data(gp_geom_1010_off + 176 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_yzz_x, g_y_0_x_0_yzz_xz, g_y_0_x_0_yzz_y, g_y_0_x_0_yzz_yz, g_y_0_x_0_yzz_z, g_y_0_x_0_yzz_zz, g_y_0_x_0_yzzz_x, g_y_0_x_0_yzzz_y, g_y_0_x_0_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yzzz_x[k] = -g_y_0_x_0_yzz_x[k] * ab_z + g_y_0_x_0_yzz_xz[k];

                g_y_0_x_0_yzzz_y[k] = -g_y_0_x_0_yzz_y[k] * ab_z + g_y_0_x_0_yzz_yz[k];

                g_y_0_x_0_yzzz_z[k] = -g_y_0_x_0_yzz_z[k] * ab_z + g_y_0_x_0_yzz_zz[k];
            }

            /// Set up 177-180 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_zzzz_x = cbuffer.data(gp_geom_1010_off + 177 * ccomps * dcomps);

            auto g_y_0_x_0_zzzz_y = cbuffer.data(gp_geom_1010_off + 178 * ccomps * dcomps);

            auto g_y_0_x_0_zzzz_z = cbuffer.data(gp_geom_1010_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_zzz_x, g_y_0_x_0_zzz_xz, g_y_0_x_0_zzz_y, g_y_0_x_0_zzz_yz, g_y_0_x_0_zzz_z, g_y_0_x_0_zzz_zz, g_y_0_x_0_zzzz_x, g_y_0_x_0_zzzz_y, g_y_0_x_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_zzzz_x[k] = -g_y_0_x_0_zzz_x[k] * ab_z + g_y_0_x_0_zzz_xz[k];

                g_y_0_x_0_zzzz_y[k] = -g_y_0_x_0_zzz_y[k] * ab_z + g_y_0_x_0_zzz_yz[k];

                g_y_0_x_0_zzzz_z[k] = -g_y_0_x_0_zzz_z[k] * ab_z + g_y_0_x_0_zzz_zz[k];
            }

            /// Set up 180-183 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxx_x = cbuffer.data(gp_geom_1010_off + 180 * ccomps * dcomps);

            auto g_y_0_y_0_xxxx_y = cbuffer.data(gp_geom_1010_off + 181 * ccomps * dcomps);

            auto g_y_0_y_0_xxxx_z = cbuffer.data(gp_geom_1010_off + 182 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxx_x, g_y_0_y_0_xxx_xx, g_y_0_y_0_xxx_xy, g_y_0_y_0_xxx_xz, g_y_0_y_0_xxx_y, g_y_0_y_0_xxx_z, g_y_0_y_0_xxxx_x, g_y_0_y_0_xxxx_y, g_y_0_y_0_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxx_x[k] = -g_y_0_y_0_xxx_x[k] * ab_x + g_y_0_y_0_xxx_xx[k];

                g_y_0_y_0_xxxx_y[k] = -g_y_0_y_0_xxx_y[k] * ab_x + g_y_0_y_0_xxx_xy[k];

                g_y_0_y_0_xxxx_z[k] = -g_y_0_y_0_xxx_z[k] * ab_x + g_y_0_y_0_xxx_xz[k];
            }

            /// Set up 183-186 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxy_x = cbuffer.data(gp_geom_1010_off + 183 * ccomps * dcomps);

            auto g_y_0_y_0_xxxy_y = cbuffer.data(gp_geom_1010_off + 184 * ccomps * dcomps);

            auto g_y_0_y_0_xxxy_z = cbuffer.data(gp_geom_1010_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxy_x, g_y_0_y_0_xxxy_y, g_y_0_y_0_xxxy_z, g_y_0_y_0_xxy_x, g_y_0_y_0_xxy_xx, g_y_0_y_0_xxy_xy, g_y_0_y_0_xxy_xz, g_y_0_y_0_xxy_y, g_y_0_y_0_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxy_x[k] = -g_y_0_y_0_xxy_x[k] * ab_x + g_y_0_y_0_xxy_xx[k];

                g_y_0_y_0_xxxy_y[k] = -g_y_0_y_0_xxy_y[k] * ab_x + g_y_0_y_0_xxy_xy[k];

                g_y_0_y_0_xxxy_z[k] = -g_y_0_y_0_xxy_z[k] * ab_x + g_y_0_y_0_xxy_xz[k];
            }

            /// Set up 186-189 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxxz_x = cbuffer.data(gp_geom_1010_off + 186 * ccomps * dcomps);

            auto g_y_0_y_0_xxxz_y = cbuffer.data(gp_geom_1010_off + 187 * ccomps * dcomps);

            auto g_y_0_y_0_xxxz_z = cbuffer.data(gp_geom_1010_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxxz_x, g_y_0_y_0_xxxz_y, g_y_0_y_0_xxxz_z, g_y_0_y_0_xxz_x, g_y_0_y_0_xxz_xx, g_y_0_y_0_xxz_xy, g_y_0_y_0_xxz_xz, g_y_0_y_0_xxz_y, g_y_0_y_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxxz_x[k] = -g_y_0_y_0_xxz_x[k] * ab_x + g_y_0_y_0_xxz_xx[k];

                g_y_0_y_0_xxxz_y[k] = -g_y_0_y_0_xxz_y[k] * ab_x + g_y_0_y_0_xxz_xy[k];

                g_y_0_y_0_xxxz_z[k] = -g_y_0_y_0_xxz_z[k] * ab_x + g_y_0_y_0_xxz_xz[k];
            }

            /// Set up 189-192 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxyy_x = cbuffer.data(gp_geom_1010_off + 189 * ccomps * dcomps);

            auto g_y_0_y_0_xxyy_y = cbuffer.data(gp_geom_1010_off + 190 * ccomps * dcomps);

            auto g_y_0_y_0_xxyy_z = cbuffer.data(gp_geom_1010_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxyy_x, g_y_0_y_0_xxyy_y, g_y_0_y_0_xxyy_z, g_y_0_y_0_xyy_x, g_y_0_y_0_xyy_xx, g_y_0_y_0_xyy_xy, g_y_0_y_0_xyy_xz, g_y_0_y_0_xyy_y, g_y_0_y_0_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxyy_x[k] = -g_y_0_y_0_xyy_x[k] * ab_x + g_y_0_y_0_xyy_xx[k];

                g_y_0_y_0_xxyy_y[k] = -g_y_0_y_0_xyy_y[k] * ab_x + g_y_0_y_0_xyy_xy[k];

                g_y_0_y_0_xxyy_z[k] = -g_y_0_y_0_xyy_z[k] * ab_x + g_y_0_y_0_xyy_xz[k];
            }

            /// Set up 192-195 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxyz_x = cbuffer.data(gp_geom_1010_off + 192 * ccomps * dcomps);

            auto g_y_0_y_0_xxyz_y = cbuffer.data(gp_geom_1010_off + 193 * ccomps * dcomps);

            auto g_y_0_y_0_xxyz_z = cbuffer.data(gp_geom_1010_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxyz_x, g_y_0_y_0_xxyz_y, g_y_0_y_0_xxyz_z, g_y_0_y_0_xyz_x, g_y_0_y_0_xyz_xx, g_y_0_y_0_xyz_xy, g_y_0_y_0_xyz_xz, g_y_0_y_0_xyz_y, g_y_0_y_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxyz_x[k] = -g_y_0_y_0_xyz_x[k] * ab_x + g_y_0_y_0_xyz_xx[k];

                g_y_0_y_0_xxyz_y[k] = -g_y_0_y_0_xyz_y[k] * ab_x + g_y_0_y_0_xyz_xy[k];

                g_y_0_y_0_xxyz_z[k] = -g_y_0_y_0_xyz_z[k] * ab_x + g_y_0_y_0_xyz_xz[k];
            }

            /// Set up 195-198 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xxzz_x = cbuffer.data(gp_geom_1010_off + 195 * ccomps * dcomps);

            auto g_y_0_y_0_xxzz_y = cbuffer.data(gp_geom_1010_off + 196 * ccomps * dcomps);

            auto g_y_0_y_0_xxzz_z = cbuffer.data(gp_geom_1010_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xxzz_x, g_y_0_y_0_xxzz_y, g_y_0_y_0_xxzz_z, g_y_0_y_0_xzz_x, g_y_0_y_0_xzz_xx, g_y_0_y_0_xzz_xy, g_y_0_y_0_xzz_xz, g_y_0_y_0_xzz_y, g_y_0_y_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xxzz_x[k] = -g_y_0_y_0_xzz_x[k] * ab_x + g_y_0_y_0_xzz_xx[k];

                g_y_0_y_0_xxzz_y[k] = -g_y_0_y_0_xzz_y[k] * ab_x + g_y_0_y_0_xzz_xy[k];

                g_y_0_y_0_xxzz_z[k] = -g_y_0_y_0_xzz_z[k] * ab_x + g_y_0_y_0_xzz_xz[k];
            }

            /// Set up 198-201 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyyy_x = cbuffer.data(gp_geom_1010_off + 198 * ccomps * dcomps);

            auto g_y_0_y_0_xyyy_y = cbuffer.data(gp_geom_1010_off + 199 * ccomps * dcomps);

            auto g_y_0_y_0_xyyy_z = cbuffer.data(gp_geom_1010_off + 200 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyyy_x, g_y_0_y_0_xyyy_y, g_y_0_y_0_xyyy_z, g_y_0_y_0_yyy_x, g_y_0_y_0_yyy_xx, g_y_0_y_0_yyy_xy, g_y_0_y_0_yyy_xz, g_y_0_y_0_yyy_y, g_y_0_y_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyyy_x[k] = -g_y_0_y_0_yyy_x[k] * ab_x + g_y_0_y_0_yyy_xx[k];

                g_y_0_y_0_xyyy_y[k] = -g_y_0_y_0_yyy_y[k] * ab_x + g_y_0_y_0_yyy_xy[k];

                g_y_0_y_0_xyyy_z[k] = -g_y_0_y_0_yyy_z[k] * ab_x + g_y_0_y_0_yyy_xz[k];
            }

            /// Set up 201-204 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyyz_x = cbuffer.data(gp_geom_1010_off + 201 * ccomps * dcomps);

            auto g_y_0_y_0_xyyz_y = cbuffer.data(gp_geom_1010_off + 202 * ccomps * dcomps);

            auto g_y_0_y_0_xyyz_z = cbuffer.data(gp_geom_1010_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyyz_x, g_y_0_y_0_xyyz_y, g_y_0_y_0_xyyz_z, g_y_0_y_0_yyz_x, g_y_0_y_0_yyz_xx, g_y_0_y_0_yyz_xy, g_y_0_y_0_yyz_xz, g_y_0_y_0_yyz_y, g_y_0_y_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyyz_x[k] = -g_y_0_y_0_yyz_x[k] * ab_x + g_y_0_y_0_yyz_xx[k];

                g_y_0_y_0_xyyz_y[k] = -g_y_0_y_0_yyz_y[k] * ab_x + g_y_0_y_0_yyz_xy[k];

                g_y_0_y_0_xyyz_z[k] = -g_y_0_y_0_yyz_z[k] * ab_x + g_y_0_y_0_yyz_xz[k];
            }

            /// Set up 204-207 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xyzz_x = cbuffer.data(gp_geom_1010_off + 204 * ccomps * dcomps);

            auto g_y_0_y_0_xyzz_y = cbuffer.data(gp_geom_1010_off + 205 * ccomps * dcomps);

            auto g_y_0_y_0_xyzz_z = cbuffer.data(gp_geom_1010_off + 206 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xyzz_x, g_y_0_y_0_xyzz_y, g_y_0_y_0_xyzz_z, g_y_0_y_0_yzz_x, g_y_0_y_0_yzz_xx, g_y_0_y_0_yzz_xy, g_y_0_y_0_yzz_xz, g_y_0_y_0_yzz_y, g_y_0_y_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xyzz_x[k] = -g_y_0_y_0_yzz_x[k] * ab_x + g_y_0_y_0_yzz_xx[k];

                g_y_0_y_0_xyzz_y[k] = -g_y_0_y_0_yzz_y[k] * ab_x + g_y_0_y_0_yzz_xy[k];

                g_y_0_y_0_xyzz_z[k] = -g_y_0_y_0_yzz_z[k] * ab_x + g_y_0_y_0_yzz_xz[k];
            }

            /// Set up 207-210 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xzzz_x = cbuffer.data(gp_geom_1010_off + 207 * ccomps * dcomps);

            auto g_y_0_y_0_xzzz_y = cbuffer.data(gp_geom_1010_off + 208 * ccomps * dcomps);

            auto g_y_0_y_0_xzzz_z = cbuffer.data(gp_geom_1010_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xzzz_x, g_y_0_y_0_xzzz_y, g_y_0_y_0_xzzz_z, g_y_0_y_0_zzz_x, g_y_0_y_0_zzz_xx, g_y_0_y_0_zzz_xy, g_y_0_y_0_zzz_xz, g_y_0_y_0_zzz_y, g_y_0_y_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xzzz_x[k] = -g_y_0_y_0_zzz_x[k] * ab_x + g_y_0_y_0_zzz_xx[k];

                g_y_0_y_0_xzzz_y[k] = -g_y_0_y_0_zzz_y[k] * ab_x + g_y_0_y_0_zzz_xy[k];

                g_y_0_y_0_xzzz_z[k] = -g_y_0_y_0_zzz_z[k] * ab_x + g_y_0_y_0_zzz_xz[k];
            }

            /// Set up 210-213 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyyy_x = cbuffer.data(gp_geom_1010_off + 210 * ccomps * dcomps);

            auto g_y_0_y_0_yyyy_y = cbuffer.data(gp_geom_1010_off + 211 * ccomps * dcomps);

            auto g_y_0_y_0_yyyy_z = cbuffer.data(gp_geom_1010_off + 212 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_yyy_x, g_0_0_y_0_yyy_y, g_0_0_y_0_yyy_z, g_y_0_y_0_yyy_x, g_y_0_y_0_yyy_xy, g_y_0_y_0_yyy_y, g_y_0_y_0_yyy_yy, g_y_0_y_0_yyy_yz, g_y_0_y_0_yyy_z, g_y_0_y_0_yyyy_x, g_y_0_y_0_yyyy_y, g_y_0_y_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyyy_x[k] = -g_0_0_y_0_yyy_x[k] - g_y_0_y_0_yyy_x[k] * ab_y + g_y_0_y_0_yyy_xy[k];

                g_y_0_y_0_yyyy_y[k] = -g_0_0_y_0_yyy_y[k] - g_y_0_y_0_yyy_y[k] * ab_y + g_y_0_y_0_yyy_yy[k];

                g_y_0_y_0_yyyy_z[k] = -g_0_0_y_0_yyy_z[k] - g_y_0_y_0_yyy_z[k] * ab_y + g_y_0_y_0_yyy_yz[k];
            }

            /// Set up 213-216 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyyz_x = cbuffer.data(gp_geom_1010_off + 213 * ccomps * dcomps);

            auto g_y_0_y_0_yyyz_y = cbuffer.data(gp_geom_1010_off + 214 * ccomps * dcomps);

            auto g_y_0_y_0_yyyz_z = cbuffer.data(gp_geom_1010_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yyy_x, g_y_0_y_0_yyy_xz, g_y_0_y_0_yyy_y, g_y_0_y_0_yyy_yz, g_y_0_y_0_yyy_z, g_y_0_y_0_yyy_zz, g_y_0_y_0_yyyz_x, g_y_0_y_0_yyyz_y, g_y_0_y_0_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyyz_x[k] = -g_y_0_y_0_yyy_x[k] * ab_z + g_y_0_y_0_yyy_xz[k];

                g_y_0_y_0_yyyz_y[k] = -g_y_0_y_0_yyy_y[k] * ab_z + g_y_0_y_0_yyy_yz[k];

                g_y_0_y_0_yyyz_z[k] = -g_y_0_y_0_yyy_z[k] * ab_z + g_y_0_y_0_yyy_zz[k];
            }

            /// Set up 216-219 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yyzz_x = cbuffer.data(gp_geom_1010_off + 216 * ccomps * dcomps);

            auto g_y_0_y_0_yyzz_y = cbuffer.data(gp_geom_1010_off + 217 * ccomps * dcomps);

            auto g_y_0_y_0_yyzz_z = cbuffer.data(gp_geom_1010_off + 218 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yyz_x, g_y_0_y_0_yyz_xz, g_y_0_y_0_yyz_y, g_y_0_y_0_yyz_yz, g_y_0_y_0_yyz_z, g_y_0_y_0_yyz_zz, g_y_0_y_0_yyzz_x, g_y_0_y_0_yyzz_y, g_y_0_y_0_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yyzz_x[k] = -g_y_0_y_0_yyz_x[k] * ab_z + g_y_0_y_0_yyz_xz[k];

                g_y_0_y_0_yyzz_y[k] = -g_y_0_y_0_yyz_y[k] * ab_z + g_y_0_y_0_yyz_yz[k];

                g_y_0_y_0_yyzz_z[k] = -g_y_0_y_0_yyz_z[k] * ab_z + g_y_0_y_0_yyz_zz[k];
            }

            /// Set up 219-222 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yzzz_x = cbuffer.data(gp_geom_1010_off + 219 * ccomps * dcomps);

            auto g_y_0_y_0_yzzz_y = cbuffer.data(gp_geom_1010_off + 220 * ccomps * dcomps);

            auto g_y_0_y_0_yzzz_z = cbuffer.data(gp_geom_1010_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_yzz_x, g_y_0_y_0_yzz_xz, g_y_0_y_0_yzz_y, g_y_0_y_0_yzz_yz, g_y_0_y_0_yzz_z, g_y_0_y_0_yzz_zz, g_y_0_y_0_yzzz_x, g_y_0_y_0_yzzz_y, g_y_0_y_0_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yzzz_x[k] = -g_y_0_y_0_yzz_x[k] * ab_z + g_y_0_y_0_yzz_xz[k];

                g_y_0_y_0_yzzz_y[k] = -g_y_0_y_0_yzz_y[k] * ab_z + g_y_0_y_0_yzz_yz[k];

                g_y_0_y_0_yzzz_z[k] = -g_y_0_y_0_yzz_z[k] * ab_z + g_y_0_y_0_yzz_zz[k];
            }

            /// Set up 222-225 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_zzzz_x = cbuffer.data(gp_geom_1010_off + 222 * ccomps * dcomps);

            auto g_y_0_y_0_zzzz_y = cbuffer.data(gp_geom_1010_off + 223 * ccomps * dcomps);

            auto g_y_0_y_0_zzzz_z = cbuffer.data(gp_geom_1010_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_zzz_x, g_y_0_y_0_zzz_xz, g_y_0_y_0_zzz_y, g_y_0_y_0_zzz_yz, g_y_0_y_0_zzz_z, g_y_0_y_0_zzz_zz, g_y_0_y_0_zzzz_x, g_y_0_y_0_zzzz_y, g_y_0_y_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_zzzz_x[k] = -g_y_0_y_0_zzz_x[k] * ab_z + g_y_0_y_0_zzz_xz[k];

                g_y_0_y_0_zzzz_y[k] = -g_y_0_y_0_zzz_y[k] * ab_z + g_y_0_y_0_zzz_yz[k];

                g_y_0_y_0_zzzz_z[k] = -g_y_0_y_0_zzz_z[k] * ab_z + g_y_0_y_0_zzz_zz[k];
            }

            /// Set up 225-228 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxx_x = cbuffer.data(gp_geom_1010_off + 225 * ccomps * dcomps);

            auto g_y_0_z_0_xxxx_y = cbuffer.data(gp_geom_1010_off + 226 * ccomps * dcomps);

            auto g_y_0_z_0_xxxx_z = cbuffer.data(gp_geom_1010_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxx_x, g_y_0_z_0_xxx_xx, g_y_0_z_0_xxx_xy, g_y_0_z_0_xxx_xz, g_y_0_z_0_xxx_y, g_y_0_z_0_xxx_z, g_y_0_z_0_xxxx_x, g_y_0_z_0_xxxx_y, g_y_0_z_0_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxx_x[k] = -g_y_0_z_0_xxx_x[k] * ab_x + g_y_0_z_0_xxx_xx[k];

                g_y_0_z_0_xxxx_y[k] = -g_y_0_z_0_xxx_y[k] * ab_x + g_y_0_z_0_xxx_xy[k];

                g_y_0_z_0_xxxx_z[k] = -g_y_0_z_0_xxx_z[k] * ab_x + g_y_0_z_0_xxx_xz[k];
            }

            /// Set up 228-231 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxy_x = cbuffer.data(gp_geom_1010_off + 228 * ccomps * dcomps);

            auto g_y_0_z_0_xxxy_y = cbuffer.data(gp_geom_1010_off + 229 * ccomps * dcomps);

            auto g_y_0_z_0_xxxy_z = cbuffer.data(gp_geom_1010_off + 230 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxy_x, g_y_0_z_0_xxxy_y, g_y_0_z_0_xxxy_z, g_y_0_z_0_xxy_x, g_y_0_z_0_xxy_xx, g_y_0_z_0_xxy_xy, g_y_0_z_0_xxy_xz, g_y_0_z_0_xxy_y, g_y_0_z_0_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxy_x[k] = -g_y_0_z_0_xxy_x[k] * ab_x + g_y_0_z_0_xxy_xx[k];

                g_y_0_z_0_xxxy_y[k] = -g_y_0_z_0_xxy_y[k] * ab_x + g_y_0_z_0_xxy_xy[k];

                g_y_0_z_0_xxxy_z[k] = -g_y_0_z_0_xxy_z[k] * ab_x + g_y_0_z_0_xxy_xz[k];
            }

            /// Set up 231-234 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxxz_x = cbuffer.data(gp_geom_1010_off + 231 * ccomps * dcomps);

            auto g_y_0_z_0_xxxz_y = cbuffer.data(gp_geom_1010_off + 232 * ccomps * dcomps);

            auto g_y_0_z_0_xxxz_z = cbuffer.data(gp_geom_1010_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxxz_x, g_y_0_z_0_xxxz_y, g_y_0_z_0_xxxz_z, g_y_0_z_0_xxz_x, g_y_0_z_0_xxz_xx, g_y_0_z_0_xxz_xy, g_y_0_z_0_xxz_xz, g_y_0_z_0_xxz_y, g_y_0_z_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxxz_x[k] = -g_y_0_z_0_xxz_x[k] * ab_x + g_y_0_z_0_xxz_xx[k];

                g_y_0_z_0_xxxz_y[k] = -g_y_0_z_0_xxz_y[k] * ab_x + g_y_0_z_0_xxz_xy[k];

                g_y_0_z_0_xxxz_z[k] = -g_y_0_z_0_xxz_z[k] * ab_x + g_y_0_z_0_xxz_xz[k];
            }

            /// Set up 234-237 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxyy_x = cbuffer.data(gp_geom_1010_off + 234 * ccomps * dcomps);

            auto g_y_0_z_0_xxyy_y = cbuffer.data(gp_geom_1010_off + 235 * ccomps * dcomps);

            auto g_y_0_z_0_xxyy_z = cbuffer.data(gp_geom_1010_off + 236 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxyy_x, g_y_0_z_0_xxyy_y, g_y_0_z_0_xxyy_z, g_y_0_z_0_xyy_x, g_y_0_z_0_xyy_xx, g_y_0_z_0_xyy_xy, g_y_0_z_0_xyy_xz, g_y_0_z_0_xyy_y, g_y_0_z_0_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxyy_x[k] = -g_y_0_z_0_xyy_x[k] * ab_x + g_y_0_z_0_xyy_xx[k];

                g_y_0_z_0_xxyy_y[k] = -g_y_0_z_0_xyy_y[k] * ab_x + g_y_0_z_0_xyy_xy[k];

                g_y_0_z_0_xxyy_z[k] = -g_y_0_z_0_xyy_z[k] * ab_x + g_y_0_z_0_xyy_xz[k];
            }

            /// Set up 237-240 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxyz_x = cbuffer.data(gp_geom_1010_off + 237 * ccomps * dcomps);

            auto g_y_0_z_0_xxyz_y = cbuffer.data(gp_geom_1010_off + 238 * ccomps * dcomps);

            auto g_y_0_z_0_xxyz_z = cbuffer.data(gp_geom_1010_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxyz_x, g_y_0_z_0_xxyz_y, g_y_0_z_0_xxyz_z, g_y_0_z_0_xyz_x, g_y_0_z_0_xyz_xx, g_y_0_z_0_xyz_xy, g_y_0_z_0_xyz_xz, g_y_0_z_0_xyz_y, g_y_0_z_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxyz_x[k] = -g_y_0_z_0_xyz_x[k] * ab_x + g_y_0_z_0_xyz_xx[k];

                g_y_0_z_0_xxyz_y[k] = -g_y_0_z_0_xyz_y[k] * ab_x + g_y_0_z_0_xyz_xy[k];

                g_y_0_z_0_xxyz_z[k] = -g_y_0_z_0_xyz_z[k] * ab_x + g_y_0_z_0_xyz_xz[k];
            }

            /// Set up 240-243 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xxzz_x = cbuffer.data(gp_geom_1010_off + 240 * ccomps * dcomps);

            auto g_y_0_z_0_xxzz_y = cbuffer.data(gp_geom_1010_off + 241 * ccomps * dcomps);

            auto g_y_0_z_0_xxzz_z = cbuffer.data(gp_geom_1010_off + 242 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xxzz_x, g_y_0_z_0_xxzz_y, g_y_0_z_0_xxzz_z, g_y_0_z_0_xzz_x, g_y_0_z_0_xzz_xx, g_y_0_z_0_xzz_xy, g_y_0_z_0_xzz_xz, g_y_0_z_0_xzz_y, g_y_0_z_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xxzz_x[k] = -g_y_0_z_0_xzz_x[k] * ab_x + g_y_0_z_0_xzz_xx[k];

                g_y_0_z_0_xxzz_y[k] = -g_y_0_z_0_xzz_y[k] * ab_x + g_y_0_z_0_xzz_xy[k];

                g_y_0_z_0_xxzz_z[k] = -g_y_0_z_0_xzz_z[k] * ab_x + g_y_0_z_0_xzz_xz[k];
            }

            /// Set up 243-246 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyyy_x = cbuffer.data(gp_geom_1010_off + 243 * ccomps * dcomps);

            auto g_y_0_z_0_xyyy_y = cbuffer.data(gp_geom_1010_off + 244 * ccomps * dcomps);

            auto g_y_0_z_0_xyyy_z = cbuffer.data(gp_geom_1010_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyyy_x, g_y_0_z_0_xyyy_y, g_y_0_z_0_xyyy_z, g_y_0_z_0_yyy_x, g_y_0_z_0_yyy_xx, g_y_0_z_0_yyy_xy, g_y_0_z_0_yyy_xz, g_y_0_z_0_yyy_y, g_y_0_z_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyyy_x[k] = -g_y_0_z_0_yyy_x[k] * ab_x + g_y_0_z_0_yyy_xx[k];

                g_y_0_z_0_xyyy_y[k] = -g_y_0_z_0_yyy_y[k] * ab_x + g_y_0_z_0_yyy_xy[k];

                g_y_0_z_0_xyyy_z[k] = -g_y_0_z_0_yyy_z[k] * ab_x + g_y_0_z_0_yyy_xz[k];
            }

            /// Set up 246-249 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyyz_x = cbuffer.data(gp_geom_1010_off + 246 * ccomps * dcomps);

            auto g_y_0_z_0_xyyz_y = cbuffer.data(gp_geom_1010_off + 247 * ccomps * dcomps);

            auto g_y_0_z_0_xyyz_z = cbuffer.data(gp_geom_1010_off + 248 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyyz_x, g_y_0_z_0_xyyz_y, g_y_0_z_0_xyyz_z, g_y_0_z_0_yyz_x, g_y_0_z_0_yyz_xx, g_y_0_z_0_yyz_xy, g_y_0_z_0_yyz_xz, g_y_0_z_0_yyz_y, g_y_0_z_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyyz_x[k] = -g_y_0_z_0_yyz_x[k] * ab_x + g_y_0_z_0_yyz_xx[k];

                g_y_0_z_0_xyyz_y[k] = -g_y_0_z_0_yyz_y[k] * ab_x + g_y_0_z_0_yyz_xy[k];

                g_y_0_z_0_xyyz_z[k] = -g_y_0_z_0_yyz_z[k] * ab_x + g_y_0_z_0_yyz_xz[k];
            }

            /// Set up 249-252 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xyzz_x = cbuffer.data(gp_geom_1010_off + 249 * ccomps * dcomps);

            auto g_y_0_z_0_xyzz_y = cbuffer.data(gp_geom_1010_off + 250 * ccomps * dcomps);

            auto g_y_0_z_0_xyzz_z = cbuffer.data(gp_geom_1010_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xyzz_x, g_y_0_z_0_xyzz_y, g_y_0_z_0_xyzz_z, g_y_0_z_0_yzz_x, g_y_0_z_0_yzz_xx, g_y_0_z_0_yzz_xy, g_y_0_z_0_yzz_xz, g_y_0_z_0_yzz_y, g_y_0_z_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xyzz_x[k] = -g_y_0_z_0_yzz_x[k] * ab_x + g_y_0_z_0_yzz_xx[k];

                g_y_0_z_0_xyzz_y[k] = -g_y_0_z_0_yzz_y[k] * ab_x + g_y_0_z_0_yzz_xy[k];

                g_y_0_z_0_xyzz_z[k] = -g_y_0_z_0_yzz_z[k] * ab_x + g_y_0_z_0_yzz_xz[k];
            }

            /// Set up 252-255 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xzzz_x = cbuffer.data(gp_geom_1010_off + 252 * ccomps * dcomps);

            auto g_y_0_z_0_xzzz_y = cbuffer.data(gp_geom_1010_off + 253 * ccomps * dcomps);

            auto g_y_0_z_0_xzzz_z = cbuffer.data(gp_geom_1010_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xzzz_x, g_y_0_z_0_xzzz_y, g_y_0_z_0_xzzz_z, g_y_0_z_0_zzz_x, g_y_0_z_0_zzz_xx, g_y_0_z_0_zzz_xy, g_y_0_z_0_zzz_xz, g_y_0_z_0_zzz_y, g_y_0_z_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xzzz_x[k] = -g_y_0_z_0_zzz_x[k] * ab_x + g_y_0_z_0_zzz_xx[k];

                g_y_0_z_0_xzzz_y[k] = -g_y_0_z_0_zzz_y[k] * ab_x + g_y_0_z_0_zzz_xy[k];

                g_y_0_z_0_xzzz_z[k] = -g_y_0_z_0_zzz_z[k] * ab_x + g_y_0_z_0_zzz_xz[k];
            }

            /// Set up 255-258 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyyy_x = cbuffer.data(gp_geom_1010_off + 255 * ccomps * dcomps);

            auto g_y_0_z_0_yyyy_y = cbuffer.data(gp_geom_1010_off + 256 * ccomps * dcomps);

            auto g_y_0_z_0_yyyy_z = cbuffer.data(gp_geom_1010_off + 257 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_yyy_x, g_0_0_z_0_yyy_y, g_0_0_z_0_yyy_z, g_y_0_z_0_yyy_x, g_y_0_z_0_yyy_xy, g_y_0_z_0_yyy_y, g_y_0_z_0_yyy_yy, g_y_0_z_0_yyy_yz, g_y_0_z_0_yyy_z, g_y_0_z_0_yyyy_x, g_y_0_z_0_yyyy_y, g_y_0_z_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyyy_x[k] = -g_0_0_z_0_yyy_x[k] - g_y_0_z_0_yyy_x[k] * ab_y + g_y_0_z_0_yyy_xy[k];

                g_y_0_z_0_yyyy_y[k] = -g_0_0_z_0_yyy_y[k] - g_y_0_z_0_yyy_y[k] * ab_y + g_y_0_z_0_yyy_yy[k];

                g_y_0_z_0_yyyy_z[k] = -g_0_0_z_0_yyy_z[k] - g_y_0_z_0_yyy_z[k] * ab_y + g_y_0_z_0_yyy_yz[k];
            }

            /// Set up 258-261 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyyz_x = cbuffer.data(gp_geom_1010_off + 258 * ccomps * dcomps);

            auto g_y_0_z_0_yyyz_y = cbuffer.data(gp_geom_1010_off + 259 * ccomps * dcomps);

            auto g_y_0_z_0_yyyz_z = cbuffer.data(gp_geom_1010_off + 260 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yyy_x, g_y_0_z_0_yyy_xz, g_y_0_z_0_yyy_y, g_y_0_z_0_yyy_yz, g_y_0_z_0_yyy_z, g_y_0_z_0_yyy_zz, g_y_0_z_0_yyyz_x, g_y_0_z_0_yyyz_y, g_y_0_z_0_yyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyyz_x[k] = -g_y_0_z_0_yyy_x[k] * ab_z + g_y_0_z_0_yyy_xz[k];

                g_y_0_z_0_yyyz_y[k] = -g_y_0_z_0_yyy_y[k] * ab_z + g_y_0_z_0_yyy_yz[k];

                g_y_0_z_0_yyyz_z[k] = -g_y_0_z_0_yyy_z[k] * ab_z + g_y_0_z_0_yyy_zz[k];
            }

            /// Set up 261-264 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yyzz_x = cbuffer.data(gp_geom_1010_off + 261 * ccomps * dcomps);

            auto g_y_0_z_0_yyzz_y = cbuffer.data(gp_geom_1010_off + 262 * ccomps * dcomps);

            auto g_y_0_z_0_yyzz_z = cbuffer.data(gp_geom_1010_off + 263 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yyz_x, g_y_0_z_0_yyz_xz, g_y_0_z_0_yyz_y, g_y_0_z_0_yyz_yz, g_y_0_z_0_yyz_z, g_y_0_z_0_yyz_zz, g_y_0_z_0_yyzz_x, g_y_0_z_0_yyzz_y, g_y_0_z_0_yyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yyzz_x[k] = -g_y_0_z_0_yyz_x[k] * ab_z + g_y_0_z_0_yyz_xz[k];

                g_y_0_z_0_yyzz_y[k] = -g_y_0_z_0_yyz_y[k] * ab_z + g_y_0_z_0_yyz_yz[k];

                g_y_0_z_0_yyzz_z[k] = -g_y_0_z_0_yyz_z[k] * ab_z + g_y_0_z_0_yyz_zz[k];
            }

            /// Set up 264-267 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yzzz_x = cbuffer.data(gp_geom_1010_off + 264 * ccomps * dcomps);

            auto g_y_0_z_0_yzzz_y = cbuffer.data(gp_geom_1010_off + 265 * ccomps * dcomps);

            auto g_y_0_z_0_yzzz_z = cbuffer.data(gp_geom_1010_off + 266 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_yzz_x, g_y_0_z_0_yzz_xz, g_y_0_z_0_yzz_y, g_y_0_z_0_yzz_yz, g_y_0_z_0_yzz_z, g_y_0_z_0_yzz_zz, g_y_0_z_0_yzzz_x, g_y_0_z_0_yzzz_y, g_y_0_z_0_yzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yzzz_x[k] = -g_y_0_z_0_yzz_x[k] * ab_z + g_y_0_z_0_yzz_xz[k];

                g_y_0_z_0_yzzz_y[k] = -g_y_0_z_0_yzz_y[k] * ab_z + g_y_0_z_0_yzz_yz[k];

                g_y_0_z_0_yzzz_z[k] = -g_y_0_z_0_yzz_z[k] * ab_z + g_y_0_z_0_yzz_zz[k];
            }

            /// Set up 267-270 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_zzzz_x = cbuffer.data(gp_geom_1010_off + 267 * ccomps * dcomps);

            auto g_y_0_z_0_zzzz_y = cbuffer.data(gp_geom_1010_off + 268 * ccomps * dcomps);

            auto g_y_0_z_0_zzzz_z = cbuffer.data(gp_geom_1010_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_zzz_x, g_y_0_z_0_zzz_xz, g_y_0_z_0_zzz_y, g_y_0_z_0_zzz_yz, g_y_0_z_0_zzz_z, g_y_0_z_0_zzz_zz, g_y_0_z_0_zzzz_x, g_y_0_z_0_zzzz_y, g_y_0_z_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_zzzz_x[k] = -g_y_0_z_0_zzz_x[k] * ab_z + g_y_0_z_0_zzz_xz[k];

                g_y_0_z_0_zzzz_y[k] = -g_y_0_z_0_zzz_y[k] * ab_z + g_y_0_z_0_zzz_yz[k];

                g_y_0_z_0_zzzz_z[k] = -g_y_0_z_0_zzz_z[k] * ab_z + g_y_0_z_0_zzz_zz[k];
            }

            /// Set up 270-273 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxx_x = cbuffer.data(gp_geom_1010_off + 270 * ccomps * dcomps);

            auto g_z_0_x_0_xxxx_y = cbuffer.data(gp_geom_1010_off + 271 * ccomps * dcomps);

            auto g_z_0_x_0_xxxx_z = cbuffer.data(gp_geom_1010_off + 272 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxx_x, g_z_0_x_0_xxx_xx, g_z_0_x_0_xxx_xy, g_z_0_x_0_xxx_xz, g_z_0_x_0_xxx_y, g_z_0_x_0_xxx_z, g_z_0_x_0_xxxx_x, g_z_0_x_0_xxxx_y, g_z_0_x_0_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxx_x[k] = -g_z_0_x_0_xxx_x[k] * ab_x + g_z_0_x_0_xxx_xx[k];

                g_z_0_x_0_xxxx_y[k] = -g_z_0_x_0_xxx_y[k] * ab_x + g_z_0_x_0_xxx_xy[k];

                g_z_0_x_0_xxxx_z[k] = -g_z_0_x_0_xxx_z[k] * ab_x + g_z_0_x_0_xxx_xz[k];
            }

            /// Set up 273-276 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxy_x = cbuffer.data(gp_geom_1010_off + 273 * ccomps * dcomps);

            auto g_z_0_x_0_xxxy_y = cbuffer.data(gp_geom_1010_off + 274 * ccomps * dcomps);

            auto g_z_0_x_0_xxxy_z = cbuffer.data(gp_geom_1010_off + 275 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxy_x, g_z_0_x_0_xxxy_y, g_z_0_x_0_xxxy_z, g_z_0_x_0_xxy_x, g_z_0_x_0_xxy_xx, g_z_0_x_0_xxy_xy, g_z_0_x_0_xxy_xz, g_z_0_x_0_xxy_y, g_z_0_x_0_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxy_x[k] = -g_z_0_x_0_xxy_x[k] * ab_x + g_z_0_x_0_xxy_xx[k];

                g_z_0_x_0_xxxy_y[k] = -g_z_0_x_0_xxy_y[k] * ab_x + g_z_0_x_0_xxy_xy[k];

                g_z_0_x_0_xxxy_z[k] = -g_z_0_x_0_xxy_z[k] * ab_x + g_z_0_x_0_xxy_xz[k];
            }

            /// Set up 276-279 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxxz_x = cbuffer.data(gp_geom_1010_off + 276 * ccomps * dcomps);

            auto g_z_0_x_0_xxxz_y = cbuffer.data(gp_geom_1010_off + 277 * ccomps * dcomps);

            auto g_z_0_x_0_xxxz_z = cbuffer.data(gp_geom_1010_off + 278 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxxz_x, g_z_0_x_0_xxxz_y, g_z_0_x_0_xxxz_z, g_z_0_x_0_xxz_x, g_z_0_x_0_xxz_xx, g_z_0_x_0_xxz_xy, g_z_0_x_0_xxz_xz, g_z_0_x_0_xxz_y, g_z_0_x_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxxz_x[k] = -g_z_0_x_0_xxz_x[k] * ab_x + g_z_0_x_0_xxz_xx[k];

                g_z_0_x_0_xxxz_y[k] = -g_z_0_x_0_xxz_y[k] * ab_x + g_z_0_x_0_xxz_xy[k];

                g_z_0_x_0_xxxz_z[k] = -g_z_0_x_0_xxz_z[k] * ab_x + g_z_0_x_0_xxz_xz[k];
            }

            /// Set up 279-282 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxyy_x = cbuffer.data(gp_geom_1010_off + 279 * ccomps * dcomps);

            auto g_z_0_x_0_xxyy_y = cbuffer.data(gp_geom_1010_off + 280 * ccomps * dcomps);

            auto g_z_0_x_0_xxyy_z = cbuffer.data(gp_geom_1010_off + 281 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxyy_x, g_z_0_x_0_xxyy_y, g_z_0_x_0_xxyy_z, g_z_0_x_0_xyy_x, g_z_0_x_0_xyy_xx, g_z_0_x_0_xyy_xy, g_z_0_x_0_xyy_xz, g_z_0_x_0_xyy_y, g_z_0_x_0_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxyy_x[k] = -g_z_0_x_0_xyy_x[k] * ab_x + g_z_0_x_0_xyy_xx[k];

                g_z_0_x_0_xxyy_y[k] = -g_z_0_x_0_xyy_y[k] * ab_x + g_z_0_x_0_xyy_xy[k];

                g_z_0_x_0_xxyy_z[k] = -g_z_0_x_0_xyy_z[k] * ab_x + g_z_0_x_0_xyy_xz[k];
            }

            /// Set up 282-285 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxyz_x = cbuffer.data(gp_geom_1010_off + 282 * ccomps * dcomps);

            auto g_z_0_x_0_xxyz_y = cbuffer.data(gp_geom_1010_off + 283 * ccomps * dcomps);

            auto g_z_0_x_0_xxyz_z = cbuffer.data(gp_geom_1010_off + 284 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxyz_x, g_z_0_x_0_xxyz_y, g_z_0_x_0_xxyz_z, g_z_0_x_0_xyz_x, g_z_0_x_0_xyz_xx, g_z_0_x_0_xyz_xy, g_z_0_x_0_xyz_xz, g_z_0_x_0_xyz_y, g_z_0_x_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxyz_x[k] = -g_z_0_x_0_xyz_x[k] * ab_x + g_z_0_x_0_xyz_xx[k];

                g_z_0_x_0_xxyz_y[k] = -g_z_0_x_0_xyz_y[k] * ab_x + g_z_0_x_0_xyz_xy[k];

                g_z_0_x_0_xxyz_z[k] = -g_z_0_x_0_xyz_z[k] * ab_x + g_z_0_x_0_xyz_xz[k];
            }

            /// Set up 285-288 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xxzz_x = cbuffer.data(gp_geom_1010_off + 285 * ccomps * dcomps);

            auto g_z_0_x_0_xxzz_y = cbuffer.data(gp_geom_1010_off + 286 * ccomps * dcomps);

            auto g_z_0_x_0_xxzz_z = cbuffer.data(gp_geom_1010_off + 287 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xxzz_x, g_z_0_x_0_xxzz_y, g_z_0_x_0_xxzz_z, g_z_0_x_0_xzz_x, g_z_0_x_0_xzz_xx, g_z_0_x_0_xzz_xy, g_z_0_x_0_xzz_xz, g_z_0_x_0_xzz_y, g_z_0_x_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xxzz_x[k] = -g_z_0_x_0_xzz_x[k] * ab_x + g_z_0_x_0_xzz_xx[k];

                g_z_0_x_0_xxzz_y[k] = -g_z_0_x_0_xzz_y[k] * ab_x + g_z_0_x_0_xzz_xy[k];

                g_z_0_x_0_xxzz_z[k] = -g_z_0_x_0_xzz_z[k] * ab_x + g_z_0_x_0_xzz_xz[k];
            }

            /// Set up 288-291 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyyy_x = cbuffer.data(gp_geom_1010_off + 288 * ccomps * dcomps);

            auto g_z_0_x_0_xyyy_y = cbuffer.data(gp_geom_1010_off + 289 * ccomps * dcomps);

            auto g_z_0_x_0_xyyy_z = cbuffer.data(gp_geom_1010_off + 290 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyyy_x, g_z_0_x_0_xyyy_y, g_z_0_x_0_xyyy_z, g_z_0_x_0_yyy_x, g_z_0_x_0_yyy_xx, g_z_0_x_0_yyy_xy, g_z_0_x_0_yyy_xz, g_z_0_x_0_yyy_y, g_z_0_x_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyyy_x[k] = -g_z_0_x_0_yyy_x[k] * ab_x + g_z_0_x_0_yyy_xx[k];

                g_z_0_x_0_xyyy_y[k] = -g_z_0_x_0_yyy_y[k] * ab_x + g_z_0_x_0_yyy_xy[k];

                g_z_0_x_0_xyyy_z[k] = -g_z_0_x_0_yyy_z[k] * ab_x + g_z_0_x_0_yyy_xz[k];
            }

            /// Set up 291-294 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyyz_x = cbuffer.data(gp_geom_1010_off + 291 * ccomps * dcomps);

            auto g_z_0_x_0_xyyz_y = cbuffer.data(gp_geom_1010_off + 292 * ccomps * dcomps);

            auto g_z_0_x_0_xyyz_z = cbuffer.data(gp_geom_1010_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyyz_x, g_z_0_x_0_xyyz_y, g_z_0_x_0_xyyz_z, g_z_0_x_0_yyz_x, g_z_0_x_0_yyz_xx, g_z_0_x_0_yyz_xy, g_z_0_x_0_yyz_xz, g_z_0_x_0_yyz_y, g_z_0_x_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyyz_x[k] = -g_z_0_x_0_yyz_x[k] * ab_x + g_z_0_x_0_yyz_xx[k];

                g_z_0_x_0_xyyz_y[k] = -g_z_0_x_0_yyz_y[k] * ab_x + g_z_0_x_0_yyz_xy[k];

                g_z_0_x_0_xyyz_z[k] = -g_z_0_x_0_yyz_z[k] * ab_x + g_z_0_x_0_yyz_xz[k];
            }

            /// Set up 294-297 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xyzz_x = cbuffer.data(gp_geom_1010_off + 294 * ccomps * dcomps);

            auto g_z_0_x_0_xyzz_y = cbuffer.data(gp_geom_1010_off + 295 * ccomps * dcomps);

            auto g_z_0_x_0_xyzz_z = cbuffer.data(gp_geom_1010_off + 296 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xyzz_x, g_z_0_x_0_xyzz_y, g_z_0_x_0_xyzz_z, g_z_0_x_0_yzz_x, g_z_0_x_0_yzz_xx, g_z_0_x_0_yzz_xy, g_z_0_x_0_yzz_xz, g_z_0_x_0_yzz_y, g_z_0_x_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xyzz_x[k] = -g_z_0_x_0_yzz_x[k] * ab_x + g_z_0_x_0_yzz_xx[k];

                g_z_0_x_0_xyzz_y[k] = -g_z_0_x_0_yzz_y[k] * ab_x + g_z_0_x_0_yzz_xy[k];

                g_z_0_x_0_xyzz_z[k] = -g_z_0_x_0_yzz_z[k] * ab_x + g_z_0_x_0_yzz_xz[k];
            }

            /// Set up 297-300 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xzzz_x = cbuffer.data(gp_geom_1010_off + 297 * ccomps * dcomps);

            auto g_z_0_x_0_xzzz_y = cbuffer.data(gp_geom_1010_off + 298 * ccomps * dcomps);

            auto g_z_0_x_0_xzzz_z = cbuffer.data(gp_geom_1010_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xzzz_x, g_z_0_x_0_xzzz_y, g_z_0_x_0_xzzz_z, g_z_0_x_0_zzz_x, g_z_0_x_0_zzz_xx, g_z_0_x_0_zzz_xy, g_z_0_x_0_zzz_xz, g_z_0_x_0_zzz_y, g_z_0_x_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xzzz_x[k] = -g_z_0_x_0_zzz_x[k] * ab_x + g_z_0_x_0_zzz_xx[k];

                g_z_0_x_0_xzzz_y[k] = -g_z_0_x_0_zzz_y[k] * ab_x + g_z_0_x_0_zzz_xy[k];

                g_z_0_x_0_xzzz_z[k] = -g_z_0_x_0_zzz_z[k] * ab_x + g_z_0_x_0_zzz_xz[k];
            }

            /// Set up 300-303 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyyy_x = cbuffer.data(gp_geom_1010_off + 300 * ccomps * dcomps);

            auto g_z_0_x_0_yyyy_y = cbuffer.data(gp_geom_1010_off + 301 * ccomps * dcomps);

            auto g_z_0_x_0_yyyy_z = cbuffer.data(gp_geom_1010_off + 302 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyy_x, g_z_0_x_0_yyy_xy, g_z_0_x_0_yyy_y, g_z_0_x_0_yyy_yy, g_z_0_x_0_yyy_yz, g_z_0_x_0_yyy_z, g_z_0_x_0_yyyy_x, g_z_0_x_0_yyyy_y, g_z_0_x_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyyy_x[k] = -g_z_0_x_0_yyy_x[k] * ab_y + g_z_0_x_0_yyy_xy[k];

                g_z_0_x_0_yyyy_y[k] = -g_z_0_x_0_yyy_y[k] * ab_y + g_z_0_x_0_yyy_yy[k];

                g_z_0_x_0_yyyy_z[k] = -g_z_0_x_0_yyy_z[k] * ab_y + g_z_0_x_0_yyy_yz[k];
            }

            /// Set up 303-306 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyyz_x = cbuffer.data(gp_geom_1010_off + 303 * ccomps * dcomps);

            auto g_z_0_x_0_yyyz_y = cbuffer.data(gp_geom_1010_off + 304 * ccomps * dcomps);

            auto g_z_0_x_0_yyyz_z = cbuffer.data(gp_geom_1010_off + 305 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyyz_x, g_z_0_x_0_yyyz_y, g_z_0_x_0_yyyz_z, g_z_0_x_0_yyz_x, g_z_0_x_0_yyz_xy, g_z_0_x_0_yyz_y, g_z_0_x_0_yyz_yy, g_z_0_x_0_yyz_yz, g_z_0_x_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyyz_x[k] = -g_z_0_x_0_yyz_x[k] * ab_y + g_z_0_x_0_yyz_xy[k];

                g_z_0_x_0_yyyz_y[k] = -g_z_0_x_0_yyz_y[k] * ab_y + g_z_0_x_0_yyz_yy[k];

                g_z_0_x_0_yyyz_z[k] = -g_z_0_x_0_yyz_z[k] * ab_y + g_z_0_x_0_yyz_yz[k];
            }

            /// Set up 306-309 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yyzz_x = cbuffer.data(gp_geom_1010_off + 306 * ccomps * dcomps);

            auto g_z_0_x_0_yyzz_y = cbuffer.data(gp_geom_1010_off + 307 * ccomps * dcomps);

            auto g_z_0_x_0_yyzz_z = cbuffer.data(gp_geom_1010_off + 308 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yyzz_x, g_z_0_x_0_yyzz_y, g_z_0_x_0_yyzz_z, g_z_0_x_0_yzz_x, g_z_0_x_0_yzz_xy, g_z_0_x_0_yzz_y, g_z_0_x_0_yzz_yy, g_z_0_x_0_yzz_yz, g_z_0_x_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yyzz_x[k] = -g_z_0_x_0_yzz_x[k] * ab_y + g_z_0_x_0_yzz_xy[k];

                g_z_0_x_0_yyzz_y[k] = -g_z_0_x_0_yzz_y[k] * ab_y + g_z_0_x_0_yzz_yy[k];

                g_z_0_x_0_yyzz_z[k] = -g_z_0_x_0_yzz_z[k] * ab_y + g_z_0_x_0_yzz_yz[k];
            }

            /// Set up 309-312 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yzzz_x = cbuffer.data(gp_geom_1010_off + 309 * ccomps * dcomps);

            auto g_z_0_x_0_yzzz_y = cbuffer.data(gp_geom_1010_off + 310 * ccomps * dcomps);

            auto g_z_0_x_0_yzzz_z = cbuffer.data(gp_geom_1010_off + 311 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yzzz_x, g_z_0_x_0_yzzz_y, g_z_0_x_0_yzzz_z, g_z_0_x_0_zzz_x, g_z_0_x_0_zzz_xy, g_z_0_x_0_zzz_y, g_z_0_x_0_zzz_yy, g_z_0_x_0_zzz_yz, g_z_0_x_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yzzz_x[k] = -g_z_0_x_0_zzz_x[k] * ab_y + g_z_0_x_0_zzz_xy[k];

                g_z_0_x_0_yzzz_y[k] = -g_z_0_x_0_zzz_y[k] * ab_y + g_z_0_x_0_zzz_yy[k];

                g_z_0_x_0_yzzz_z[k] = -g_z_0_x_0_zzz_z[k] * ab_y + g_z_0_x_0_zzz_yz[k];
            }

            /// Set up 312-315 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_zzzz_x = cbuffer.data(gp_geom_1010_off + 312 * ccomps * dcomps);

            auto g_z_0_x_0_zzzz_y = cbuffer.data(gp_geom_1010_off + 313 * ccomps * dcomps);

            auto g_z_0_x_0_zzzz_z = cbuffer.data(gp_geom_1010_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_zzz_x, g_0_0_x_0_zzz_y, g_0_0_x_0_zzz_z, g_z_0_x_0_zzz_x, g_z_0_x_0_zzz_xz, g_z_0_x_0_zzz_y, g_z_0_x_0_zzz_yz, g_z_0_x_0_zzz_z, g_z_0_x_0_zzz_zz, g_z_0_x_0_zzzz_x, g_z_0_x_0_zzzz_y, g_z_0_x_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_zzzz_x[k] = -g_0_0_x_0_zzz_x[k] - g_z_0_x_0_zzz_x[k] * ab_z + g_z_0_x_0_zzz_xz[k];

                g_z_0_x_0_zzzz_y[k] = -g_0_0_x_0_zzz_y[k] - g_z_0_x_0_zzz_y[k] * ab_z + g_z_0_x_0_zzz_yz[k];

                g_z_0_x_0_zzzz_z[k] = -g_0_0_x_0_zzz_z[k] - g_z_0_x_0_zzz_z[k] * ab_z + g_z_0_x_0_zzz_zz[k];
            }

            /// Set up 315-318 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxx_x = cbuffer.data(gp_geom_1010_off + 315 * ccomps * dcomps);

            auto g_z_0_y_0_xxxx_y = cbuffer.data(gp_geom_1010_off + 316 * ccomps * dcomps);

            auto g_z_0_y_0_xxxx_z = cbuffer.data(gp_geom_1010_off + 317 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxx_x, g_z_0_y_0_xxx_xx, g_z_0_y_0_xxx_xy, g_z_0_y_0_xxx_xz, g_z_0_y_0_xxx_y, g_z_0_y_0_xxx_z, g_z_0_y_0_xxxx_x, g_z_0_y_0_xxxx_y, g_z_0_y_0_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxx_x[k] = -g_z_0_y_0_xxx_x[k] * ab_x + g_z_0_y_0_xxx_xx[k];

                g_z_0_y_0_xxxx_y[k] = -g_z_0_y_0_xxx_y[k] * ab_x + g_z_0_y_0_xxx_xy[k];

                g_z_0_y_0_xxxx_z[k] = -g_z_0_y_0_xxx_z[k] * ab_x + g_z_0_y_0_xxx_xz[k];
            }

            /// Set up 318-321 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxy_x = cbuffer.data(gp_geom_1010_off + 318 * ccomps * dcomps);

            auto g_z_0_y_0_xxxy_y = cbuffer.data(gp_geom_1010_off + 319 * ccomps * dcomps);

            auto g_z_0_y_0_xxxy_z = cbuffer.data(gp_geom_1010_off + 320 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxy_x, g_z_0_y_0_xxxy_y, g_z_0_y_0_xxxy_z, g_z_0_y_0_xxy_x, g_z_0_y_0_xxy_xx, g_z_0_y_0_xxy_xy, g_z_0_y_0_xxy_xz, g_z_0_y_0_xxy_y, g_z_0_y_0_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxy_x[k] = -g_z_0_y_0_xxy_x[k] * ab_x + g_z_0_y_0_xxy_xx[k];

                g_z_0_y_0_xxxy_y[k] = -g_z_0_y_0_xxy_y[k] * ab_x + g_z_0_y_0_xxy_xy[k];

                g_z_0_y_0_xxxy_z[k] = -g_z_0_y_0_xxy_z[k] * ab_x + g_z_0_y_0_xxy_xz[k];
            }

            /// Set up 321-324 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxxz_x = cbuffer.data(gp_geom_1010_off + 321 * ccomps * dcomps);

            auto g_z_0_y_0_xxxz_y = cbuffer.data(gp_geom_1010_off + 322 * ccomps * dcomps);

            auto g_z_0_y_0_xxxz_z = cbuffer.data(gp_geom_1010_off + 323 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxxz_x, g_z_0_y_0_xxxz_y, g_z_0_y_0_xxxz_z, g_z_0_y_0_xxz_x, g_z_0_y_0_xxz_xx, g_z_0_y_0_xxz_xy, g_z_0_y_0_xxz_xz, g_z_0_y_0_xxz_y, g_z_0_y_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxxz_x[k] = -g_z_0_y_0_xxz_x[k] * ab_x + g_z_0_y_0_xxz_xx[k];

                g_z_0_y_0_xxxz_y[k] = -g_z_0_y_0_xxz_y[k] * ab_x + g_z_0_y_0_xxz_xy[k];

                g_z_0_y_0_xxxz_z[k] = -g_z_0_y_0_xxz_z[k] * ab_x + g_z_0_y_0_xxz_xz[k];
            }

            /// Set up 324-327 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxyy_x = cbuffer.data(gp_geom_1010_off + 324 * ccomps * dcomps);

            auto g_z_0_y_0_xxyy_y = cbuffer.data(gp_geom_1010_off + 325 * ccomps * dcomps);

            auto g_z_0_y_0_xxyy_z = cbuffer.data(gp_geom_1010_off + 326 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxyy_x, g_z_0_y_0_xxyy_y, g_z_0_y_0_xxyy_z, g_z_0_y_0_xyy_x, g_z_0_y_0_xyy_xx, g_z_0_y_0_xyy_xy, g_z_0_y_0_xyy_xz, g_z_0_y_0_xyy_y, g_z_0_y_0_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxyy_x[k] = -g_z_0_y_0_xyy_x[k] * ab_x + g_z_0_y_0_xyy_xx[k];

                g_z_0_y_0_xxyy_y[k] = -g_z_0_y_0_xyy_y[k] * ab_x + g_z_0_y_0_xyy_xy[k];

                g_z_0_y_0_xxyy_z[k] = -g_z_0_y_0_xyy_z[k] * ab_x + g_z_0_y_0_xyy_xz[k];
            }

            /// Set up 327-330 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxyz_x = cbuffer.data(gp_geom_1010_off + 327 * ccomps * dcomps);

            auto g_z_0_y_0_xxyz_y = cbuffer.data(gp_geom_1010_off + 328 * ccomps * dcomps);

            auto g_z_0_y_0_xxyz_z = cbuffer.data(gp_geom_1010_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxyz_x, g_z_0_y_0_xxyz_y, g_z_0_y_0_xxyz_z, g_z_0_y_0_xyz_x, g_z_0_y_0_xyz_xx, g_z_0_y_0_xyz_xy, g_z_0_y_0_xyz_xz, g_z_0_y_0_xyz_y, g_z_0_y_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxyz_x[k] = -g_z_0_y_0_xyz_x[k] * ab_x + g_z_0_y_0_xyz_xx[k];

                g_z_0_y_0_xxyz_y[k] = -g_z_0_y_0_xyz_y[k] * ab_x + g_z_0_y_0_xyz_xy[k];

                g_z_0_y_0_xxyz_z[k] = -g_z_0_y_0_xyz_z[k] * ab_x + g_z_0_y_0_xyz_xz[k];
            }

            /// Set up 330-333 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xxzz_x = cbuffer.data(gp_geom_1010_off + 330 * ccomps * dcomps);

            auto g_z_0_y_0_xxzz_y = cbuffer.data(gp_geom_1010_off + 331 * ccomps * dcomps);

            auto g_z_0_y_0_xxzz_z = cbuffer.data(gp_geom_1010_off + 332 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xxzz_x, g_z_0_y_0_xxzz_y, g_z_0_y_0_xxzz_z, g_z_0_y_0_xzz_x, g_z_0_y_0_xzz_xx, g_z_0_y_0_xzz_xy, g_z_0_y_0_xzz_xz, g_z_0_y_0_xzz_y, g_z_0_y_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xxzz_x[k] = -g_z_0_y_0_xzz_x[k] * ab_x + g_z_0_y_0_xzz_xx[k];

                g_z_0_y_0_xxzz_y[k] = -g_z_0_y_0_xzz_y[k] * ab_x + g_z_0_y_0_xzz_xy[k];

                g_z_0_y_0_xxzz_z[k] = -g_z_0_y_0_xzz_z[k] * ab_x + g_z_0_y_0_xzz_xz[k];
            }

            /// Set up 333-336 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyyy_x = cbuffer.data(gp_geom_1010_off + 333 * ccomps * dcomps);

            auto g_z_0_y_0_xyyy_y = cbuffer.data(gp_geom_1010_off + 334 * ccomps * dcomps);

            auto g_z_0_y_0_xyyy_z = cbuffer.data(gp_geom_1010_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyyy_x, g_z_0_y_0_xyyy_y, g_z_0_y_0_xyyy_z, g_z_0_y_0_yyy_x, g_z_0_y_0_yyy_xx, g_z_0_y_0_yyy_xy, g_z_0_y_0_yyy_xz, g_z_0_y_0_yyy_y, g_z_0_y_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyyy_x[k] = -g_z_0_y_0_yyy_x[k] * ab_x + g_z_0_y_0_yyy_xx[k];

                g_z_0_y_0_xyyy_y[k] = -g_z_0_y_0_yyy_y[k] * ab_x + g_z_0_y_0_yyy_xy[k];

                g_z_0_y_0_xyyy_z[k] = -g_z_0_y_0_yyy_z[k] * ab_x + g_z_0_y_0_yyy_xz[k];
            }

            /// Set up 336-339 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyyz_x = cbuffer.data(gp_geom_1010_off + 336 * ccomps * dcomps);

            auto g_z_0_y_0_xyyz_y = cbuffer.data(gp_geom_1010_off + 337 * ccomps * dcomps);

            auto g_z_0_y_0_xyyz_z = cbuffer.data(gp_geom_1010_off + 338 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyyz_x, g_z_0_y_0_xyyz_y, g_z_0_y_0_xyyz_z, g_z_0_y_0_yyz_x, g_z_0_y_0_yyz_xx, g_z_0_y_0_yyz_xy, g_z_0_y_0_yyz_xz, g_z_0_y_0_yyz_y, g_z_0_y_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyyz_x[k] = -g_z_0_y_0_yyz_x[k] * ab_x + g_z_0_y_0_yyz_xx[k];

                g_z_0_y_0_xyyz_y[k] = -g_z_0_y_0_yyz_y[k] * ab_x + g_z_0_y_0_yyz_xy[k];

                g_z_0_y_0_xyyz_z[k] = -g_z_0_y_0_yyz_z[k] * ab_x + g_z_0_y_0_yyz_xz[k];
            }

            /// Set up 339-342 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xyzz_x = cbuffer.data(gp_geom_1010_off + 339 * ccomps * dcomps);

            auto g_z_0_y_0_xyzz_y = cbuffer.data(gp_geom_1010_off + 340 * ccomps * dcomps);

            auto g_z_0_y_0_xyzz_z = cbuffer.data(gp_geom_1010_off + 341 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xyzz_x, g_z_0_y_0_xyzz_y, g_z_0_y_0_xyzz_z, g_z_0_y_0_yzz_x, g_z_0_y_0_yzz_xx, g_z_0_y_0_yzz_xy, g_z_0_y_0_yzz_xz, g_z_0_y_0_yzz_y, g_z_0_y_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xyzz_x[k] = -g_z_0_y_0_yzz_x[k] * ab_x + g_z_0_y_0_yzz_xx[k];

                g_z_0_y_0_xyzz_y[k] = -g_z_0_y_0_yzz_y[k] * ab_x + g_z_0_y_0_yzz_xy[k];

                g_z_0_y_0_xyzz_z[k] = -g_z_0_y_0_yzz_z[k] * ab_x + g_z_0_y_0_yzz_xz[k];
            }

            /// Set up 342-345 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xzzz_x = cbuffer.data(gp_geom_1010_off + 342 * ccomps * dcomps);

            auto g_z_0_y_0_xzzz_y = cbuffer.data(gp_geom_1010_off + 343 * ccomps * dcomps);

            auto g_z_0_y_0_xzzz_z = cbuffer.data(gp_geom_1010_off + 344 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xzzz_x, g_z_0_y_0_xzzz_y, g_z_0_y_0_xzzz_z, g_z_0_y_0_zzz_x, g_z_0_y_0_zzz_xx, g_z_0_y_0_zzz_xy, g_z_0_y_0_zzz_xz, g_z_0_y_0_zzz_y, g_z_0_y_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xzzz_x[k] = -g_z_0_y_0_zzz_x[k] * ab_x + g_z_0_y_0_zzz_xx[k];

                g_z_0_y_0_xzzz_y[k] = -g_z_0_y_0_zzz_y[k] * ab_x + g_z_0_y_0_zzz_xy[k];

                g_z_0_y_0_xzzz_z[k] = -g_z_0_y_0_zzz_z[k] * ab_x + g_z_0_y_0_zzz_xz[k];
            }

            /// Set up 345-348 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyyy_x = cbuffer.data(gp_geom_1010_off + 345 * ccomps * dcomps);

            auto g_z_0_y_0_yyyy_y = cbuffer.data(gp_geom_1010_off + 346 * ccomps * dcomps);

            auto g_z_0_y_0_yyyy_z = cbuffer.data(gp_geom_1010_off + 347 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyy_x, g_z_0_y_0_yyy_xy, g_z_0_y_0_yyy_y, g_z_0_y_0_yyy_yy, g_z_0_y_0_yyy_yz, g_z_0_y_0_yyy_z, g_z_0_y_0_yyyy_x, g_z_0_y_0_yyyy_y, g_z_0_y_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyyy_x[k] = -g_z_0_y_0_yyy_x[k] * ab_y + g_z_0_y_0_yyy_xy[k];

                g_z_0_y_0_yyyy_y[k] = -g_z_0_y_0_yyy_y[k] * ab_y + g_z_0_y_0_yyy_yy[k];

                g_z_0_y_0_yyyy_z[k] = -g_z_0_y_0_yyy_z[k] * ab_y + g_z_0_y_0_yyy_yz[k];
            }

            /// Set up 348-351 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyyz_x = cbuffer.data(gp_geom_1010_off + 348 * ccomps * dcomps);

            auto g_z_0_y_0_yyyz_y = cbuffer.data(gp_geom_1010_off + 349 * ccomps * dcomps);

            auto g_z_0_y_0_yyyz_z = cbuffer.data(gp_geom_1010_off + 350 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyyz_x, g_z_0_y_0_yyyz_y, g_z_0_y_0_yyyz_z, g_z_0_y_0_yyz_x, g_z_0_y_0_yyz_xy, g_z_0_y_0_yyz_y, g_z_0_y_0_yyz_yy, g_z_0_y_0_yyz_yz, g_z_0_y_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyyz_x[k] = -g_z_0_y_0_yyz_x[k] * ab_y + g_z_0_y_0_yyz_xy[k];

                g_z_0_y_0_yyyz_y[k] = -g_z_0_y_0_yyz_y[k] * ab_y + g_z_0_y_0_yyz_yy[k];

                g_z_0_y_0_yyyz_z[k] = -g_z_0_y_0_yyz_z[k] * ab_y + g_z_0_y_0_yyz_yz[k];
            }

            /// Set up 351-354 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yyzz_x = cbuffer.data(gp_geom_1010_off + 351 * ccomps * dcomps);

            auto g_z_0_y_0_yyzz_y = cbuffer.data(gp_geom_1010_off + 352 * ccomps * dcomps);

            auto g_z_0_y_0_yyzz_z = cbuffer.data(gp_geom_1010_off + 353 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yyzz_x, g_z_0_y_0_yyzz_y, g_z_0_y_0_yyzz_z, g_z_0_y_0_yzz_x, g_z_0_y_0_yzz_xy, g_z_0_y_0_yzz_y, g_z_0_y_0_yzz_yy, g_z_0_y_0_yzz_yz, g_z_0_y_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yyzz_x[k] = -g_z_0_y_0_yzz_x[k] * ab_y + g_z_0_y_0_yzz_xy[k];

                g_z_0_y_0_yyzz_y[k] = -g_z_0_y_0_yzz_y[k] * ab_y + g_z_0_y_0_yzz_yy[k];

                g_z_0_y_0_yyzz_z[k] = -g_z_0_y_0_yzz_z[k] * ab_y + g_z_0_y_0_yzz_yz[k];
            }

            /// Set up 354-357 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yzzz_x = cbuffer.data(gp_geom_1010_off + 354 * ccomps * dcomps);

            auto g_z_0_y_0_yzzz_y = cbuffer.data(gp_geom_1010_off + 355 * ccomps * dcomps);

            auto g_z_0_y_0_yzzz_z = cbuffer.data(gp_geom_1010_off + 356 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yzzz_x, g_z_0_y_0_yzzz_y, g_z_0_y_0_yzzz_z, g_z_0_y_0_zzz_x, g_z_0_y_0_zzz_xy, g_z_0_y_0_zzz_y, g_z_0_y_0_zzz_yy, g_z_0_y_0_zzz_yz, g_z_0_y_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yzzz_x[k] = -g_z_0_y_0_zzz_x[k] * ab_y + g_z_0_y_0_zzz_xy[k];

                g_z_0_y_0_yzzz_y[k] = -g_z_0_y_0_zzz_y[k] * ab_y + g_z_0_y_0_zzz_yy[k];

                g_z_0_y_0_yzzz_z[k] = -g_z_0_y_0_zzz_z[k] * ab_y + g_z_0_y_0_zzz_yz[k];
            }

            /// Set up 357-360 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_zzzz_x = cbuffer.data(gp_geom_1010_off + 357 * ccomps * dcomps);

            auto g_z_0_y_0_zzzz_y = cbuffer.data(gp_geom_1010_off + 358 * ccomps * dcomps);

            auto g_z_0_y_0_zzzz_z = cbuffer.data(gp_geom_1010_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_zzz_x, g_0_0_y_0_zzz_y, g_0_0_y_0_zzz_z, g_z_0_y_0_zzz_x, g_z_0_y_0_zzz_xz, g_z_0_y_0_zzz_y, g_z_0_y_0_zzz_yz, g_z_0_y_0_zzz_z, g_z_0_y_0_zzz_zz, g_z_0_y_0_zzzz_x, g_z_0_y_0_zzzz_y, g_z_0_y_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_zzzz_x[k] = -g_0_0_y_0_zzz_x[k] - g_z_0_y_0_zzz_x[k] * ab_z + g_z_0_y_0_zzz_xz[k];

                g_z_0_y_0_zzzz_y[k] = -g_0_0_y_0_zzz_y[k] - g_z_0_y_0_zzz_y[k] * ab_z + g_z_0_y_0_zzz_yz[k];

                g_z_0_y_0_zzzz_z[k] = -g_0_0_y_0_zzz_z[k] - g_z_0_y_0_zzz_z[k] * ab_z + g_z_0_y_0_zzz_zz[k];
            }

            /// Set up 360-363 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxx_x = cbuffer.data(gp_geom_1010_off + 360 * ccomps * dcomps);

            auto g_z_0_z_0_xxxx_y = cbuffer.data(gp_geom_1010_off + 361 * ccomps * dcomps);

            auto g_z_0_z_0_xxxx_z = cbuffer.data(gp_geom_1010_off + 362 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxx_x, g_z_0_z_0_xxx_xx, g_z_0_z_0_xxx_xy, g_z_0_z_0_xxx_xz, g_z_0_z_0_xxx_y, g_z_0_z_0_xxx_z, g_z_0_z_0_xxxx_x, g_z_0_z_0_xxxx_y, g_z_0_z_0_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxx_x[k] = -g_z_0_z_0_xxx_x[k] * ab_x + g_z_0_z_0_xxx_xx[k];

                g_z_0_z_0_xxxx_y[k] = -g_z_0_z_0_xxx_y[k] * ab_x + g_z_0_z_0_xxx_xy[k];

                g_z_0_z_0_xxxx_z[k] = -g_z_0_z_0_xxx_z[k] * ab_x + g_z_0_z_0_xxx_xz[k];
            }

            /// Set up 363-366 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxy_x = cbuffer.data(gp_geom_1010_off + 363 * ccomps * dcomps);

            auto g_z_0_z_0_xxxy_y = cbuffer.data(gp_geom_1010_off + 364 * ccomps * dcomps);

            auto g_z_0_z_0_xxxy_z = cbuffer.data(gp_geom_1010_off + 365 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxy_x, g_z_0_z_0_xxxy_y, g_z_0_z_0_xxxy_z, g_z_0_z_0_xxy_x, g_z_0_z_0_xxy_xx, g_z_0_z_0_xxy_xy, g_z_0_z_0_xxy_xz, g_z_0_z_0_xxy_y, g_z_0_z_0_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxy_x[k] = -g_z_0_z_0_xxy_x[k] * ab_x + g_z_0_z_0_xxy_xx[k];

                g_z_0_z_0_xxxy_y[k] = -g_z_0_z_0_xxy_y[k] * ab_x + g_z_0_z_0_xxy_xy[k];

                g_z_0_z_0_xxxy_z[k] = -g_z_0_z_0_xxy_z[k] * ab_x + g_z_0_z_0_xxy_xz[k];
            }

            /// Set up 366-369 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxxz_x = cbuffer.data(gp_geom_1010_off + 366 * ccomps * dcomps);

            auto g_z_0_z_0_xxxz_y = cbuffer.data(gp_geom_1010_off + 367 * ccomps * dcomps);

            auto g_z_0_z_0_xxxz_z = cbuffer.data(gp_geom_1010_off + 368 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxxz_x, g_z_0_z_0_xxxz_y, g_z_0_z_0_xxxz_z, g_z_0_z_0_xxz_x, g_z_0_z_0_xxz_xx, g_z_0_z_0_xxz_xy, g_z_0_z_0_xxz_xz, g_z_0_z_0_xxz_y, g_z_0_z_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxxz_x[k] = -g_z_0_z_0_xxz_x[k] * ab_x + g_z_0_z_0_xxz_xx[k];

                g_z_0_z_0_xxxz_y[k] = -g_z_0_z_0_xxz_y[k] * ab_x + g_z_0_z_0_xxz_xy[k];

                g_z_0_z_0_xxxz_z[k] = -g_z_0_z_0_xxz_z[k] * ab_x + g_z_0_z_0_xxz_xz[k];
            }

            /// Set up 369-372 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxyy_x = cbuffer.data(gp_geom_1010_off + 369 * ccomps * dcomps);

            auto g_z_0_z_0_xxyy_y = cbuffer.data(gp_geom_1010_off + 370 * ccomps * dcomps);

            auto g_z_0_z_0_xxyy_z = cbuffer.data(gp_geom_1010_off + 371 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxyy_x, g_z_0_z_0_xxyy_y, g_z_0_z_0_xxyy_z, g_z_0_z_0_xyy_x, g_z_0_z_0_xyy_xx, g_z_0_z_0_xyy_xy, g_z_0_z_0_xyy_xz, g_z_0_z_0_xyy_y, g_z_0_z_0_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxyy_x[k] = -g_z_0_z_0_xyy_x[k] * ab_x + g_z_0_z_0_xyy_xx[k];

                g_z_0_z_0_xxyy_y[k] = -g_z_0_z_0_xyy_y[k] * ab_x + g_z_0_z_0_xyy_xy[k];

                g_z_0_z_0_xxyy_z[k] = -g_z_0_z_0_xyy_z[k] * ab_x + g_z_0_z_0_xyy_xz[k];
            }

            /// Set up 372-375 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxyz_x = cbuffer.data(gp_geom_1010_off + 372 * ccomps * dcomps);

            auto g_z_0_z_0_xxyz_y = cbuffer.data(gp_geom_1010_off + 373 * ccomps * dcomps);

            auto g_z_0_z_0_xxyz_z = cbuffer.data(gp_geom_1010_off + 374 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxyz_x, g_z_0_z_0_xxyz_y, g_z_0_z_0_xxyz_z, g_z_0_z_0_xyz_x, g_z_0_z_0_xyz_xx, g_z_0_z_0_xyz_xy, g_z_0_z_0_xyz_xz, g_z_0_z_0_xyz_y, g_z_0_z_0_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxyz_x[k] = -g_z_0_z_0_xyz_x[k] * ab_x + g_z_0_z_0_xyz_xx[k];

                g_z_0_z_0_xxyz_y[k] = -g_z_0_z_0_xyz_y[k] * ab_x + g_z_0_z_0_xyz_xy[k];

                g_z_0_z_0_xxyz_z[k] = -g_z_0_z_0_xyz_z[k] * ab_x + g_z_0_z_0_xyz_xz[k];
            }

            /// Set up 375-378 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xxzz_x = cbuffer.data(gp_geom_1010_off + 375 * ccomps * dcomps);

            auto g_z_0_z_0_xxzz_y = cbuffer.data(gp_geom_1010_off + 376 * ccomps * dcomps);

            auto g_z_0_z_0_xxzz_z = cbuffer.data(gp_geom_1010_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xxzz_x, g_z_0_z_0_xxzz_y, g_z_0_z_0_xxzz_z, g_z_0_z_0_xzz_x, g_z_0_z_0_xzz_xx, g_z_0_z_0_xzz_xy, g_z_0_z_0_xzz_xz, g_z_0_z_0_xzz_y, g_z_0_z_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xxzz_x[k] = -g_z_0_z_0_xzz_x[k] * ab_x + g_z_0_z_0_xzz_xx[k];

                g_z_0_z_0_xxzz_y[k] = -g_z_0_z_0_xzz_y[k] * ab_x + g_z_0_z_0_xzz_xy[k];

                g_z_0_z_0_xxzz_z[k] = -g_z_0_z_0_xzz_z[k] * ab_x + g_z_0_z_0_xzz_xz[k];
            }

            /// Set up 378-381 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyyy_x = cbuffer.data(gp_geom_1010_off + 378 * ccomps * dcomps);

            auto g_z_0_z_0_xyyy_y = cbuffer.data(gp_geom_1010_off + 379 * ccomps * dcomps);

            auto g_z_0_z_0_xyyy_z = cbuffer.data(gp_geom_1010_off + 380 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyyy_x, g_z_0_z_0_xyyy_y, g_z_0_z_0_xyyy_z, g_z_0_z_0_yyy_x, g_z_0_z_0_yyy_xx, g_z_0_z_0_yyy_xy, g_z_0_z_0_yyy_xz, g_z_0_z_0_yyy_y, g_z_0_z_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyyy_x[k] = -g_z_0_z_0_yyy_x[k] * ab_x + g_z_0_z_0_yyy_xx[k];

                g_z_0_z_0_xyyy_y[k] = -g_z_0_z_0_yyy_y[k] * ab_x + g_z_0_z_0_yyy_xy[k];

                g_z_0_z_0_xyyy_z[k] = -g_z_0_z_0_yyy_z[k] * ab_x + g_z_0_z_0_yyy_xz[k];
            }

            /// Set up 381-384 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyyz_x = cbuffer.data(gp_geom_1010_off + 381 * ccomps * dcomps);

            auto g_z_0_z_0_xyyz_y = cbuffer.data(gp_geom_1010_off + 382 * ccomps * dcomps);

            auto g_z_0_z_0_xyyz_z = cbuffer.data(gp_geom_1010_off + 383 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyyz_x, g_z_0_z_0_xyyz_y, g_z_0_z_0_xyyz_z, g_z_0_z_0_yyz_x, g_z_0_z_0_yyz_xx, g_z_0_z_0_yyz_xy, g_z_0_z_0_yyz_xz, g_z_0_z_0_yyz_y, g_z_0_z_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyyz_x[k] = -g_z_0_z_0_yyz_x[k] * ab_x + g_z_0_z_0_yyz_xx[k];

                g_z_0_z_0_xyyz_y[k] = -g_z_0_z_0_yyz_y[k] * ab_x + g_z_0_z_0_yyz_xy[k];

                g_z_0_z_0_xyyz_z[k] = -g_z_0_z_0_yyz_z[k] * ab_x + g_z_0_z_0_yyz_xz[k];
            }

            /// Set up 384-387 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xyzz_x = cbuffer.data(gp_geom_1010_off + 384 * ccomps * dcomps);

            auto g_z_0_z_0_xyzz_y = cbuffer.data(gp_geom_1010_off + 385 * ccomps * dcomps);

            auto g_z_0_z_0_xyzz_z = cbuffer.data(gp_geom_1010_off + 386 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xyzz_x, g_z_0_z_0_xyzz_y, g_z_0_z_0_xyzz_z, g_z_0_z_0_yzz_x, g_z_0_z_0_yzz_xx, g_z_0_z_0_yzz_xy, g_z_0_z_0_yzz_xz, g_z_0_z_0_yzz_y, g_z_0_z_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xyzz_x[k] = -g_z_0_z_0_yzz_x[k] * ab_x + g_z_0_z_0_yzz_xx[k];

                g_z_0_z_0_xyzz_y[k] = -g_z_0_z_0_yzz_y[k] * ab_x + g_z_0_z_0_yzz_xy[k];

                g_z_0_z_0_xyzz_z[k] = -g_z_0_z_0_yzz_z[k] * ab_x + g_z_0_z_0_yzz_xz[k];
            }

            /// Set up 387-390 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xzzz_x = cbuffer.data(gp_geom_1010_off + 387 * ccomps * dcomps);

            auto g_z_0_z_0_xzzz_y = cbuffer.data(gp_geom_1010_off + 388 * ccomps * dcomps);

            auto g_z_0_z_0_xzzz_z = cbuffer.data(gp_geom_1010_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xzzz_x, g_z_0_z_0_xzzz_y, g_z_0_z_0_xzzz_z, g_z_0_z_0_zzz_x, g_z_0_z_0_zzz_xx, g_z_0_z_0_zzz_xy, g_z_0_z_0_zzz_xz, g_z_0_z_0_zzz_y, g_z_0_z_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xzzz_x[k] = -g_z_0_z_0_zzz_x[k] * ab_x + g_z_0_z_0_zzz_xx[k];

                g_z_0_z_0_xzzz_y[k] = -g_z_0_z_0_zzz_y[k] * ab_x + g_z_0_z_0_zzz_xy[k];

                g_z_0_z_0_xzzz_z[k] = -g_z_0_z_0_zzz_z[k] * ab_x + g_z_0_z_0_zzz_xz[k];
            }

            /// Set up 390-393 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyyy_x = cbuffer.data(gp_geom_1010_off + 390 * ccomps * dcomps);

            auto g_z_0_z_0_yyyy_y = cbuffer.data(gp_geom_1010_off + 391 * ccomps * dcomps);

            auto g_z_0_z_0_yyyy_z = cbuffer.data(gp_geom_1010_off + 392 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyy_x, g_z_0_z_0_yyy_xy, g_z_0_z_0_yyy_y, g_z_0_z_0_yyy_yy, g_z_0_z_0_yyy_yz, g_z_0_z_0_yyy_z, g_z_0_z_0_yyyy_x, g_z_0_z_0_yyyy_y, g_z_0_z_0_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyyy_x[k] = -g_z_0_z_0_yyy_x[k] * ab_y + g_z_0_z_0_yyy_xy[k];

                g_z_0_z_0_yyyy_y[k] = -g_z_0_z_0_yyy_y[k] * ab_y + g_z_0_z_0_yyy_yy[k];

                g_z_0_z_0_yyyy_z[k] = -g_z_0_z_0_yyy_z[k] * ab_y + g_z_0_z_0_yyy_yz[k];
            }

            /// Set up 393-396 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyyz_x = cbuffer.data(gp_geom_1010_off + 393 * ccomps * dcomps);

            auto g_z_0_z_0_yyyz_y = cbuffer.data(gp_geom_1010_off + 394 * ccomps * dcomps);

            auto g_z_0_z_0_yyyz_z = cbuffer.data(gp_geom_1010_off + 395 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyyz_x, g_z_0_z_0_yyyz_y, g_z_0_z_0_yyyz_z, g_z_0_z_0_yyz_x, g_z_0_z_0_yyz_xy, g_z_0_z_0_yyz_y, g_z_0_z_0_yyz_yy, g_z_0_z_0_yyz_yz, g_z_0_z_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyyz_x[k] = -g_z_0_z_0_yyz_x[k] * ab_y + g_z_0_z_0_yyz_xy[k];

                g_z_0_z_0_yyyz_y[k] = -g_z_0_z_0_yyz_y[k] * ab_y + g_z_0_z_0_yyz_yy[k];

                g_z_0_z_0_yyyz_z[k] = -g_z_0_z_0_yyz_z[k] * ab_y + g_z_0_z_0_yyz_yz[k];
            }

            /// Set up 396-399 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yyzz_x = cbuffer.data(gp_geom_1010_off + 396 * ccomps * dcomps);

            auto g_z_0_z_0_yyzz_y = cbuffer.data(gp_geom_1010_off + 397 * ccomps * dcomps);

            auto g_z_0_z_0_yyzz_z = cbuffer.data(gp_geom_1010_off + 398 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yyzz_x, g_z_0_z_0_yyzz_y, g_z_0_z_0_yyzz_z, g_z_0_z_0_yzz_x, g_z_0_z_0_yzz_xy, g_z_0_z_0_yzz_y, g_z_0_z_0_yzz_yy, g_z_0_z_0_yzz_yz, g_z_0_z_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yyzz_x[k] = -g_z_0_z_0_yzz_x[k] * ab_y + g_z_0_z_0_yzz_xy[k];

                g_z_0_z_0_yyzz_y[k] = -g_z_0_z_0_yzz_y[k] * ab_y + g_z_0_z_0_yzz_yy[k];

                g_z_0_z_0_yyzz_z[k] = -g_z_0_z_0_yzz_z[k] * ab_y + g_z_0_z_0_yzz_yz[k];
            }

            /// Set up 399-402 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yzzz_x = cbuffer.data(gp_geom_1010_off + 399 * ccomps * dcomps);

            auto g_z_0_z_0_yzzz_y = cbuffer.data(gp_geom_1010_off + 400 * ccomps * dcomps);

            auto g_z_0_z_0_yzzz_z = cbuffer.data(gp_geom_1010_off + 401 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yzzz_x, g_z_0_z_0_yzzz_y, g_z_0_z_0_yzzz_z, g_z_0_z_0_zzz_x, g_z_0_z_0_zzz_xy, g_z_0_z_0_zzz_y, g_z_0_z_0_zzz_yy, g_z_0_z_0_zzz_yz, g_z_0_z_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yzzz_x[k] = -g_z_0_z_0_zzz_x[k] * ab_y + g_z_0_z_0_zzz_xy[k];

                g_z_0_z_0_yzzz_y[k] = -g_z_0_z_0_zzz_y[k] * ab_y + g_z_0_z_0_zzz_yy[k];

                g_z_0_z_0_yzzz_z[k] = -g_z_0_z_0_zzz_z[k] * ab_y + g_z_0_z_0_zzz_yz[k];
            }

            /// Set up 402-405 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_zzzz_x = cbuffer.data(gp_geom_1010_off + 402 * ccomps * dcomps);

            auto g_z_0_z_0_zzzz_y = cbuffer.data(gp_geom_1010_off + 403 * ccomps * dcomps);

            auto g_z_0_z_0_zzzz_z = cbuffer.data(gp_geom_1010_off + 404 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_zzz_x, g_0_0_z_0_zzz_y, g_0_0_z_0_zzz_z, g_z_0_z_0_zzz_x, g_z_0_z_0_zzz_xz, g_z_0_z_0_zzz_y, g_z_0_z_0_zzz_yz, g_z_0_z_0_zzz_z, g_z_0_z_0_zzz_zz, g_z_0_z_0_zzzz_x, g_z_0_z_0_zzzz_y, g_z_0_z_0_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_zzzz_x[k] = -g_0_0_z_0_zzz_x[k] - g_z_0_z_0_zzz_x[k] * ab_z + g_z_0_z_0_zzz_xz[k];

                g_z_0_z_0_zzzz_y[k] = -g_0_0_z_0_zzz_y[k] - g_z_0_z_0_zzz_y[k] * ab_z + g_z_0_z_0_zzz_yz[k];

                g_z_0_z_0_zzzz_z[k] = -g_0_0_z_0_zzz_z[k] - g_z_0_z_0_zzz_z[k] * ab_z + g_z_0_z_0_zzz_zz[k];
            }
        }
    }
}

} // erirec namespace

