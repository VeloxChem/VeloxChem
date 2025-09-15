#include "ElectronRepulsionGeom2000ContrRecFSXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_fsxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_fsxx,
                                            const size_t idx_geom_10_dsxx,
                                            const size_t idx_geom_20_dsxx,
                                            const size_t idx_geom_20_dpxx,
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
            /// Set up components of auxilary buffer : DSSS

            const auto ds_geom_10_off = idx_geom_10_dsxx + i * dcomps + j;

            auto g_x_0_xx_0 = cbuffer.data(ds_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xy_0 = cbuffer.data(ds_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xz_0 = cbuffer.data(ds_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_yy_0 = cbuffer.data(ds_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_yz_0 = cbuffer.data(ds_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_zz_0 = cbuffer.data(ds_geom_10_off + 5 * ccomps * dcomps);

            auto g_y_0_xx_0 = cbuffer.data(ds_geom_10_off + 6 * ccomps * dcomps);

            auto g_y_0_xy_0 = cbuffer.data(ds_geom_10_off + 7 * ccomps * dcomps);

            auto g_y_0_xz_0 = cbuffer.data(ds_geom_10_off + 8 * ccomps * dcomps);

            auto g_y_0_yy_0 = cbuffer.data(ds_geom_10_off + 9 * ccomps * dcomps);

            auto g_y_0_yz_0 = cbuffer.data(ds_geom_10_off + 10 * ccomps * dcomps);

            auto g_y_0_zz_0 = cbuffer.data(ds_geom_10_off + 11 * ccomps * dcomps);

            auto g_z_0_xx_0 = cbuffer.data(ds_geom_10_off + 12 * ccomps * dcomps);

            auto g_z_0_xy_0 = cbuffer.data(ds_geom_10_off + 13 * ccomps * dcomps);

            auto g_z_0_xz_0 = cbuffer.data(ds_geom_10_off + 14 * ccomps * dcomps);

            auto g_z_0_yy_0 = cbuffer.data(ds_geom_10_off + 15 * ccomps * dcomps);

            auto g_z_0_yz_0 = cbuffer.data(ds_geom_10_off + 16 * ccomps * dcomps);

            auto g_z_0_zz_0 = cbuffer.data(ds_geom_10_off + 17 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DSSS

            const auto ds_geom_20_off = idx_geom_20_dsxx + i * dcomps + j;

            auto g_xx_0_xx_0 = cbuffer.data(ds_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xy_0 = cbuffer.data(ds_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xz_0 = cbuffer.data(ds_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_yy_0 = cbuffer.data(ds_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_yz_0 = cbuffer.data(ds_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_zz_0 = cbuffer.data(ds_geom_20_off + 5 * ccomps * dcomps);

            auto g_xy_0_xx_0 = cbuffer.data(ds_geom_20_off + 6 * ccomps * dcomps);

            auto g_xy_0_xy_0 = cbuffer.data(ds_geom_20_off + 7 * ccomps * dcomps);

            auto g_xy_0_xz_0 = cbuffer.data(ds_geom_20_off + 8 * ccomps * dcomps);

            auto g_xy_0_yy_0 = cbuffer.data(ds_geom_20_off + 9 * ccomps * dcomps);

            auto g_xy_0_yz_0 = cbuffer.data(ds_geom_20_off + 10 * ccomps * dcomps);

            auto g_xy_0_zz_0 = cbuffer.data(ds_geom_20_off + 11 * ccomps * dcomps);

            auto g_xz_0_xx_0 = cbuffer.data(ds_geom_20_off + 12 * ccomps * dcomps);

            auto g_xz_0_xy_0 = cbuffer.data(ds_geom_20_off + 13 * ccomps * dcomps);

            auto g_xz_0_xz_0 = cbuffer.data(ds_geom_20_off + 14 * ccomps * dcomps);

            auto g_xz_0_yy_0 = cbuffer.data(ds_geom_20_off + 15 * ccomps * dcomps);

            auto g_xz_0_yz_0 = cbuffer.data(ds_geom_20_off + 16 * ccomps * dcomps);

            auto g_xz_0_zz_0 = cbuffer.data(ds_geom_20_off + 17 * ccomps * dcomps);

            auto g_yy_0_xx_0 = cbuffer.data(ds_geom_20_off + 18 * ccomps * dcomps);

            auto g_yy_0_xy_0 = cbuffer.data(ds_geom_20_off + 19 * ccomps * dcomps);

            auto g_yy_0_xz_0 = cbuffer.data(ds_geom_20_off + 20 * ccomps * dcomps);

            auto g_yy_0_yy_0 = cbuffer.data(ds_geom_20_off + 21 * ccomps * dcomps);

            auto g_yy_0_yz_0 = cbuffer.data(ds_geom_20_off + 22 * ccomps * dcomps);

            auto g_yy_0_zz_0 = cbuffer.data(ds_geom_20_off + 23 * ccomps * dcomps);

            auto g_yz_0_xx_0 = cbuffer.data(ds_geom_20_off + 24 * ccomps * dcomps);

            auto g_yz_0_xy_0 = cbuffer.data(ds_geom_20_off + 25 * ccomps * dcomps);

            auto g_yz_0_xz_0 = cbuffer.data(ds_geom_20_off + 26 * ccomps * dcomps);

            auto g_yz_0_yy_0 = cbuffer.data(ds_geom_20_off + 27 * ccomps * dcomps);

            auto g_yz_0_yz_0 = cbuffer.data(ds_geom_20_off + 28 * ccomps * dcomps);

            auto g_yz_0_zz_0 = cbuffer.data(ds_geom_20_off + 29 * ccomps * dcomps);

            auto g_zz_0_xx_0 = cbuffer.data(ds_geom_20_off + 30 * ccomps * dcomps);

            auto g_zz_0_xy_0 = cbuffer.data(ds_geom_20_off + 31 * ccomps * dcomps);

            auto g_zz_0_xz_0 = cbuffer.data(ds_geom_20_off + 32 * ccomps * dcomps);

            auto g_zz_0_yy_0 = cbuffer.data(ds_geom_20_off + 33 * ccomps * dcomps);

            auto g_zz_0_yz_0 = cbuffer.data(ds_geom_20_off + 34 * ccomps * dcomps);

            auto g_zz_0_zz_0 = cbuffer.data(ds_geom_20_off + 35 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DPSS

            const auto dp_geom_20_off = idx_geom_20_dpxx + i * dcomps + j;

            auto g_xx_0_xx_x = cbuffer.data(dp_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xx_y = cbuffer.data(dp_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xx_z = cbuffer.data(dp_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xy_x = cbuffer.data(dp_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xy_y = cbuffer.data(dp_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xy_z = cbuffer.data(dp_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xz_x = cbuffer.data(dp_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xz_y = cbuffer.data(dp_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xz_z = cbuffer.data(dp_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_yy_x = cbuffer.data(dp_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_yy_y = cbuffer.data(dp_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_yy_z = cbuffer.data(dp_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_yz_x = cbuffer.data(dp_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_yz_y = cbuffer.data(dp_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_yz_z = cbuffer.data(dp_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_zz_x = cbuffer.data(dp_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_zz_y = cbuffer.data(dp_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_zz_z = cbuffer.data(dp_geom_20_off + 17 * ccomps * dcomps);

            auto g_xy_0_xx_x = cbuffer.data(dp_geom_20_off + 18 * ccomps * dcomps);

            auto g_xy_0_xx_y = cbuffer.data(dp_geom_20_off + 19 * ccomps * dcomps);

            auto g_xy_0_xx_z = cbuffer.data(dp_geom_20_off + 20 * ccomps * dcomps);

            auto g_xy_0_xy_x = cbuffer.data(dp_geom_20_off + 21 * ccomps * dcomps);

            auto g_xy_0_xy_y = cbuffer.data(dp_geom_20_off + 22 * ccomps * dcomps);

            auto g_xy_0_xy_z = cbuffer.data(dp_geom_20_off + 23 * ccomps * dcomps);

            auto g_xy_0_xz_x = cbuffer.data(dp_geom_20_off + 24 * ccomps * dcomps);

            auto g_xy_0_xz_y = cbuffer.data(dp_geom_20_off + 25 * ccomps * dcomps);

            auto g_xy_0_xz_z = cbuffer.data(dp_geom_20_off + 26 * ccomps * dcomps);

            auto g_xy_0_yy_x = cbuffer.data(dp_geom_20_off + 27 * ccomps * dcomps);

            auto g_xy_0_yy_y = cbuffer.data(dp_geom_20_off + 28 * ccomps * dcomps);

            auto g_xy_0_yy_z = cbuffer.data(dp_geom_20_off + 29 * ccomps * dcomps);

            auto g_xy_0_yz_x = cbuffer.data(dp_geom_20_off + 30 * ccomps * dcomps);

            auto g_xy_0_yz_y = cbuffer.data(dp_geom_20_off + 31 * ccomps * dcomps);

            auto g_xy_0_yz_z = cbuffer.data(dp_geom_20_off + 32 * ccomps * dcomps);

            auto g_xy_0_zz_x = cbuffer.data(dp_geom_20_off + 33 * ccomps * dcomps);

            auto g_xy_0_zz_y = cbuffer.data(dp_geom_20_off + 34 * ccomps * dcomps);

            auto g_xy_0_zz_z = cbuffer.data(dp_geom_20_off + 35 * ccomps * dcomps);

            auto g_xz_0_xx_x = cbuffer.data(dp_geom_20_off + 36 * ccomps * dcomps);

            auto g_xz_0_xx_y = cbuffer.data(dp_geom_20_off + 37 * ccomps * dcomps);

            auto g_xz_0_xx_z = cbuffer.data(dp_geom_20_off + 38 * ccomps * dcomps);

            auto g_xz_0_xy_x = cbuffer.data(dp_geom_20_off + 39 * ccomps * dcomps);

            auto g_xz_0_xy_y = cbuffer.data(dp_geom_20_off + 40 * ccomps * dcomps);

            auto g_xz_0_xy_z = cbuffer.data(dp_geom_20_off + 41 * ccomps * dcomps);

            auto g_xz_0_xz_x = cbuffer.data(dp_geom_20_off + 42 * ccomps * dcomps);

            auto g_xz_0_xz_y = cbuffer.data(dp_geom_20_off + 43 * ccomps * dcomps);

            auto g_xz_0_xz_z = cbuffer.data(dp_geom_20_off + 44 * ccomps * dcomps);

            auto g_xz_0_yy_x = cbuffer.data(dp_geom_20_off + 45 * ccomps * dcomps);

            auto g_xz_0_yy_y = cbuffer.data(dp_geom_20_off + 46 * ccomps * dcomps);

            auto g_xz_0_yy_z = cbuffer.data(dp_geom_20_off + 47 * ccomps * dcomps);

            auto g_xz_0_yz_x = cbuffer.data(dp_geom_20_off + 48 * ccomps * dcomps);

            auto g_xz_0_yz_y = cbuffer.data(dp_geom_20_off + 49 * ccomps * dcomps);

            auto g_xz_0_yz_z = cbuffer.data(dp_geom_20_off + 50 * ccomps * dcomps);

            auto g_xz_0_zz_x = cbuffer.data(dp_geom_20_off + 51 * ccomps * dcomps);

            auto g_xz_0_zz_y = cbuffer.data(dp_geom_20_off + 52 * ccomps * dcomps);

            auto g_xz_0_zz_z = cbuffer.data(dp_geom_20_off + 53 * ccomps * dcomps);

            auto g_yy_0_xx_x = cbuffer.data(dp_geom_20_off + 54 * ccomps * dcomps);

            auto g_yy_0_xx_y = cbuffer.data(dp_geom_20_off + 55 * ccomps * dcomps);

            auto g_yy_0_xx_z = cbuffer.data(dp_geom_20_off + 56 * ccomps * dcomps);

            auto g_yy_0_xy_x = cbuffer.data(dp_geom_20_off + 57 * ccomps * dcomps);

            auto g_yy_0_xy_y = cbuffer.data(dp_geom_20_off + 58 * ccomps * dcomps);

            auto g_yy_0_xy_z = cbuffer.data(dp_geom_20_off + 59 * ccomps * dcomps);

            auto g_yy_0_xz_x = cbuffer.data(dp_geom_20_off + 60 * ccomps * dcomps);

            auto g_yy_0_xz_y = cbuffer.data(dp_geom_20_off + 61 * ccomps * dcomps);

            auto g_yy_0_xz_z = cbuffer.data(dp_geom_20_off + 62 * ccomps * dcomps);

            auto g_yy_0_yy_x = cbuffer.data(dp_geom_20_off + 63 * ccomps * dcomps);

            auto g_yy_0_yy_y = cbuffer.data(dp_geom_20_off + 64 * ccomps * dcomps);

            auto g_yy_0_yy_z = cbuffer.data(dp_geom_20_off + 65 * ccomps * dcomps);

            auto g_yy_0_yz_x = cbuffer.data(dp_geom_20_off + 66 * ccomps * dcomps);

            auto g_yy_0_yz_y = cbuffer.data(dp_geom_20_off + 67 * ccomps * dcomps);

            auto g_yy_0_yz_z = cbuffer.data(dp_geom_20_off + 68 * ccomps * dcomps);

            auto g_yy_0_zz_x = cbuffer.data(dp_geom_20_off + 69 * ccomps * dcomps);

            auto g_yy_0_zz_y = cbuffer.data(dp_geom_20_off + 70 * ccomps * dcomps);

            auto g_yy_0_zz_z = cbuffer.data(dp_geom_20_off + 71 * ccomps * dcomps);

            auto g_yz_0_xx_x = cbuffer.data(dp_geom_20_off + 72 * ccomps * dcomps);

            auto g_yz_0_xx_y = cbuffer.data(dp_geom_20_off + 73 * ccomps * dcomps);

            auto g_yz_0_xx_z = cbuffer.data(dp_geom_20_off + 74 * ccomps * dcomps);

            auto g_yz_0_xy_x = cbuffer.data(dp_geom_20_off + 75 * ccomps * dcomps);

            auto g_yz_0_xy_y = cbuffer.data(dp_geom_20_off + 76 * ccomps * dcomps);

            auto g_yz_0_xy_z = cbuffer.data(dp_geom_20_off + 77 * ccomps * dcomps);

            auto g_yz_0_xz_x = cbuffer.data(dp_geom_20_off + 78 * ccomps * dcomps);

            auto g_yz_0_xz_y = cbuffer.data(dp_geom_20_off + 79 * ccomps * dcomps);

            auto g_yz_0_xz_z = cbuffer.data(dp_geom_20_off + 80 * ccomps * dcomps);

            auto g_yz_0_yy_x = cbuffer.data(dp_geom_20_off + 81 * ccomps * dcomps);

            auto g_yz_0_yy_y = cbuffer.data(dp_geom_20_off + 82 * ccomps * dcomps);

            auto g_yz_0_yy_z = cbuffer.data(dp_geom_20_off + 83 * ccomps * dcomps);

            auto g_yz_0_yz_x = cbuffer.data(dp_geom_20_off + 84 * ccomps * dcomps);

            auto g_yz_0_yz_y = cbuffer.data(dp_geom_20_off + 85 * ccomps * dcomps);

            auto g_yz_0_yz_z = cbuffer.data(dp_geom_20_off + 86 * ccomps * dcomps);

            auto g_yz_0_zz_x = cbuffer.data(dp_geom_20_off + 87 * ccomps * dcomps);

            auto g_yz_0_zz_y = cbuffer.data(dp_geom_20_off + 88 * ccomps * dcomps);

            auto g_yz_0_zz_z = cbuffer.data(dp_geom_20_off + 89 * ccomps * dcomps);

            auto g_zz_0_xx_x = cbuffer.data(dp_geom_20_off + 90 * ccomps * dcomps);

            auto g_zz_0_xx_y = cbuffer.data(dp_geom_20_off + 91 * ccomps * dcomps);

            auto g_zz_0_xx_z = cbuffer.data(dp_geom_20_off + 92 * ccomps * dcomps);

            auto g_zz_0_xy_x = cbuffer.data(dp_geom_20_off + 93 * ccomps * dcomps);

            auto g_zz_0_xy_y = cbuffer.data(dp_geom_20_off + 94 * ccomps * dcomps);

            auto g_zz_0_xy_z = cbuffer.data(dp_geom_20_off + 95 * ccomps * dcomps);

            auto g_zz_0_xz_x = cbuffer.data(dp_geom_20_off + 96 * ccomps * dcomps);

            auto g_zz_0_xz_y = cbuffer.data(dp_geom_20_off + 97 * ccomps * dcomps);

            auto g_zz_0_xz_z = cbuffer.data(dp_geom_20_off + 98 * ccomps * dcomps);

            auto g_zz_0_yy_x = cbuffer.data(dp_geom_20_off + 99 * ccomps * dcomps);

            auto g_zz_0_yy_y = cbuffer.data(dp_geom_20_off + 100 * ccomps * dcomps);

            auto g_zz_0_yy_z = cbuffer.data(dp_geom_20_off + 101 * ccomps * dcomps);

            auto g_zz_0_yz_x = cbuffer.data(dp_geom_20_off + 102 * ccomps * dcomps);

            auto g_zz_0_yz_y = cbuffer.data(dp_geom_20_off + 103 * ccomps * dcomps);

            auto g_zz_0_yz_z = cbuffer.data(dp_geom_20_off + 104 * ccomps * dcomps);

            auto g_zz_0_zz_x = cbuffer.data(dp_geom_20_off + 105 * ccomps * dcomps);

            auto g_zz_0_zz_y = cbuffer.data(dp_geom_20_off + 106 * ccomps * dcomps);

            auto g_zz_0_zz_z = cbuffer.data(dp_geom_20_off + 107 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_fsxx

            const auto fs_geom_20_off = idx_geom_20_fsxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxx_0 = cbuffer.data(fs_geom_20_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xx_0, g_xx_0_xx_0, g_xx_0_xx_x, g_xx_0_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxx_0[k] = -2.0 * g_x_0_xx_0[k] - g_xx_0_xx_0[k] * ab_x + g_xx_0_xx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxy_0 = cbuffer.data(fs_geom_20_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xx_0, g_xx_0_xx_y, g_xx_0_xxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxy_0[k] = -g_xx_0_xx_0[k] * ab_y + g_xx_0_xx_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxz_0 = cbuffer.data(fs_geom_20_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xx_0, g_xx_0_xx_z, g_xx_0_xxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxz_0[k] = -g_xx_0_xx_0[k] * ab_z + g_xx_0_xx_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyy_0 = cbuffer.data(fs_geom_20_off + 3 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xy_0, g_xx_0_xy_y, g_xx_0_xyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyy_0[k] = -g_xx_0_xy_0[k] * ab_y + g_xx_0_xy_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyz_0 = cbuffer.data(fs_geom_20_off + 4 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyz_0, g_xx_0_xz_0, g_xx_0_xz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyz_0[k] = -g_xx_0_xz_0[k] * ab_y + g_xx_0_xz_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xzz_0 = cbuffer.data(fs_geom_20_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xz_0, g_xx_0_xz_z, g_xx_0_xzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xzz_0[k] = -g_xx_0_xz_0[k] * ab_z + g_xx_0_xz_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyy_0 = cbuffer.data(fs_geom_20_off + 6 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yy_0, g_xx_0_yy_y, g_xx_0_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyy_0[k] = -g_xx_0_yy_0[k] * ab_y + g_xx_0_yy_y[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyz_0 = cbuffer.data(fs_geom_20_off + 7 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyz_0, g_xx_0_yz_0, g_xx_0_yz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyz_0[k] = -g_xx_0_yz_0[k] * ab_y + g_xx_0_yz_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yzz_0 = cbuffer.data(fs_geom_20_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yzz_0, g_xx_0_zz_0, g_xx_0_zz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yzz_0[k] = -g_xx_0_zz_0[k] * ab_y + g_xx_0_zz_y[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_xx_0_zzz_0 = cbuffer.data(fs_geom_20_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_zz_0, g_xx_0_zz_z, g_xx_0_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_zzz_0[k] = -g_xx_0_zz_0[k] * ab_z + g_xx_0_zz_z[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxx_0 = cbuffer.data(fs_geom_20_off + 10 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xx_0, g_xy_0_xx_x, g_xy_0_xxx_0, g_y_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxx_0[k] = -g_y_0_xx_0[k] - g_xy_0_xx_0[k] * ab_x + g_xy_0_xx_x[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxy_0 = cbuffer.data(fs_geom_20_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxy_0, g_xy_0_xy_0, g_xy_0_xy_x, g_y_0_xy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxy_0[k] = -g_y_0_xy_0[k] - g_xy_0_xy_0[k] * ab_x + g_xy_0_xy_x[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxz_0 = cbuffer.data(fs_geom_20_off + 12 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xx_0, g_xy_0_xx_z, g_xy_0_xxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxz_0[k] = -g_xy_0_xx_0[k] * ab_z + g_xy_0_xx_z[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyy_0 = cbuffer.data(fs_geom_20_off + 13 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyy_0, g_xy_0_yy_0, g_xy_0_yy_x, g_y_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyy_0[k] = -g_y_0_yy_0[k] - g_xy_0_yy_0[k] * ab_x + g_xy_0_yy_x[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyz_0 = cbuffer.data(fs_geom_20_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xy_0, g_xy_0_xy_z, g_xy_0_xyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyz_0[k] = -g_xy_0_xy_0[k] * ab_z + g_xy_0_xy_z[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xzz_0 = cbuffer.data(fs_geom_20_off + 15 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xz_0, g_xy_0_xz_z, g_xy_0_xzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xzz_0[k] = -g_xy_0_xz_0[k] * ab_z + g_xy_0_xz_z[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyy_0 = cbuffer.data(fs_geom_20_off + 16 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yy_0, g_xy_0_yy_0, g_xy_0_yy_y, g_xy_0_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyy_0[k] = -g_x_0_yy_0[k] - g_xy_0_yy_0[k] * ab_y + g_xy_0_yy_y[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyz_0 = cbuffer.data(fs_geom_20_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yy_0, g_xy_0_yy_z, g_xy_0_yyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyz_0[k] = -g_xy_0_yy_0[k] * ab_z + g_xy_0_yy_z[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yzz_0 = cbuffer.data(fs_geom_20_off + 18 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yz_0, g_xy_0_yz_z, g_xy_0_yzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yzz_0[k] = -g_xy_0_yz_0[k] * ab_z + g_xy_0_yz_z[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_xy_0_zzz_0 = cbuffer.data(fs_geom_20_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_zz_0, g_xy_0_zz_z, g_xy_0_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_zzz_0[k] = -g_xy_0_zz_0[k] * ab_z + g_xy_0_zz_z[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxx_0 = cbuffer.data(fs_geom_20_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xx_0, g_xz_0_xx_x, g_xz_0_xxx_0, g_z_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxx_0[k] = -g_z_0_xx_0[k] - g_xz_0_xx_0[k] * ab_x + g_xz_0_xx_x[k];
            }

            /// Set up 21-22 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxy_0 = cbuffer.data(fs_geom_20_off + 21 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xx_0, g_xz_0_xx_y, g_xz_0_xxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxy_0[k] = -g_xz_0_xx_0[k] * ab_y + g_xz_0_xx_y[k];
            }

            /// Set up 22-23 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxz_0 = cbuffer.data(fs_geom_20_off + 22 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxz_0, g_xz_0_xz_0, g_xz_0_xz_x, g_z_0_xz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxz_0[k] = -g_z_0_xz_0[k] - g_xz_0_xz_0[k] * ab_x + g_xz_0_xz_x[k];
            }

            /// Set up 23-24 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyy_0 = cbuffer.data(fs_geom_20_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xy_0, g_xz_0_xy_y, g_xz_0_xyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyy_0[k] = -g_xz_0_xy_0[k] * ab_y + g_xz_0_xy_y[k];
            }

            /// Set up 24-25 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyz_0 = cbuffer.data(fs_geom_20_off + 24 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyz_0, g_xz_0_xz_0, g_xz_0_xz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyz_0[k] = -g_xz_0_xz_0[k] * ab_y + g_xz_0_xz_y[k];
            }

            /// Set up 25-26 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xzz_0 = cbuffer.data(fs_geom_20_off + 25 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xzz_0, g_xz_0_zz_0, g_xz_0_zz_x, g_z_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xzz_0[k] = -g_z_0_zz_0[k] - g_xz_0_zz_0[k] * ab_x + g_xz_0_zz_x[k];
            }

            /// Set up 26-27 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyy_0 = cbuffer.data(fs_geom_20_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yy_0, g_xz_0_yy_y, g_xz_0_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyy_0[k] = -g_xz_0_yy_0[k] * ab_y + g_xz_0_yy_y[k];
            }

            /// Set up 27-28 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyz_0 = cbuffer.data(fs_geom_20_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyz_0, g_xz_0_yz_0, g_xz_0_yz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyz_0[k] = -g_xz_0_yz_0[k] * ab_y + g_xz_0_yz_y[k];
            }

            /// Set up 28-29 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yzz_0 = cbuffer.data(fs_geom_20_off + 28 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yzz_0, g_xz_0_zz_0, g_xz_0_zz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yzz_0[k] = -g_xz_0_zz_0[k] * ab_y + g_xz_0_zz_y[k];
            }

            /// Set up 29-30 components of targeted buffer : cbuffer.data(

            auto g_xz_0_zzz_0 = cbuffer.data(fs_geom_20_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zz_0, g_xz_0_zz_0, g_xz_0_zz_z, g_xz_0_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_zzz_0[k] = -g_x_0_zz_0[k] - g_xz_0_zz_0[k] * ab_z + g_xz_0_zz_z[k];
            }

            /// Set up 30-31 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxx_0 = cbuffer.data(fs_geom_20_off + 30 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xx_0, g_yy_0_xx_x, g_yy_0_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxx_0[k] = -g_yy_0_xx_0[k] * ab_x + g_yy_0_xx_x[k];
            }

            /// Set up 31-32 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxy_0 = cbuffer.data(fs_geom_20_off + 31 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxy_0, g_yy_0_xy_0, g_yy_0_xy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxy_0[k] = -g_yy_0_xy_0[k] * ab_x + g_yy_0_xy_x[k];
            }

            /// Set up 32-33 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxz_0 = cbuffer.data(fs_geom_20_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxz_0, g_yy_0_xz_0, g_yy_0_xz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxz_0[k] = -g_yy_0_xz_0[k] * ab_x + g_yy_0_xz_x[k];
            }

            /// Set up 33-34 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyy_0 = cbuffer.data(fs_geom_20_off + 33 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyy_0, g_yy_0_yy_0, g_yy_0_yy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyy_0[k] = -g_yy_0_yy_0[k] * ab_x + g_yy_0_yy_x[k];
            }

            /// Set up 34-35 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyz_0 = cbuffer.data(fs_geom_20_off + 34 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyz_0, g_yy_0_yz_0, g_yy_0_yz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyz_0[k] = -g_yy_0_yz_0[k] * ab_x + g_yy_0_yz_x[k];
            }

            /// Set up 35-36 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xzz_0 = cbuffer.data(fs_geom_20_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xzz_0, g_yy_0_zz_0, g_yy_0_zz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xzz_0[k] = -g_yy_0_zz_0[k] * ab_x + g_yy_0_zz_x[k];
            }

            /// Set up 36-37 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyy_0 = cbuffer.data(fs_geom_20_off + 36 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yy_0, g_yy_0_yy_0, g_yy_0_yy_y, g_yy_0_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyy_0[k] = -2.0 * g_y_0_yy_0[k] - g_yy_0_yy_0[k] * ab_y + g_yy_0_yy_y[k];
            }

            /// Set up 37-38 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyz_0 = cbuffer.data(fs_geom_20_off + 37 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yy_0, g_yy_0_yy_z, g_yy_0_yyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyz_0[k] = -g_yy_0_yy_0[k] * ab_z + g_yy_0_yy_z[k];
            }

            /// Set up 38-39 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yzz_0 = cbuffer.data(fs_geom_20_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yz_0, g_yy_0_yz_z, g_yy_0_yzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yzz_0[k] = -g_yy_0_yz_0[k] * ab_z + g_yy_0_yz_z[k];
            }

            /// Set up 39-40 components of targeted buffer : cbuffer.data(

            auto g_yy_0_zzz_0 = cbuffer.data(fs_geom_20_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_zz_0, g_yy_0_zz_z, g_yy_0_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_zzz_0[k] = -g_yy_0_zz_0[k] * ab_z + g_yy_0_zz_z[k];
            }

            /// Set up 40-41 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxx_0 = cbuffer.data(fs_geom_20_off + 40 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xx_0, g_yz_0_xx_x, g_yz_0_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxx_0[k] = -g_yz_0_xx_0[k] * ab_x + g_yz_0_xx_x[k];
            }

            /// Set up 41-42 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxy_0 = cbuffer.data(fs_geom_20_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxy_0, g_yz_0_xy_0, g_yz_0_xy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxy_0[k] = -g_yz_0_xy_0[k] * ab_x + g_yz_0_xy_x[k];
            }

            /// Set up 42-43 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxz_0 = cbuffer.data(fs_geom_20_off + 42 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxz_0, g_yz_0_xz_0, g_yz_0_xz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxz_0[k] = -g_yz_0_xz_0[k] * ab_x + g_yz_0_xz_x[k];
            }

            /// Set up 43-44 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyy_0 = cbuffer.data(fs_geom_20_off + 43 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyy_0, g_yz_0_yy_0, g_yz_0_yy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyy_0[k] = -g_yz_0_yy_0[k] * ab_x + g_yz_0_yy_x[k];
            }

            /// Set up 44-45 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyz_0 = cbuffer.data(fs_geom_20_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyz_0, g_yz_0_yz_0, g_yz_0_yz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyz_0[k] = -g_yz_0_yz_0[k] * ab_x + g_yz_0_yz_x[k];
            }

            /// Set up 45-46 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xzz_0 = cbuffer.data(fs_geom_20_off + 45 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xzz_0, g_yz_0_zz_0, g_yz_0_zz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xzz_0[k] = -g_yz_0_zz_0[k] * ab_x + g_yz_0_zz_x[k];
            }

            /// Set up 46-47 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyy_0 = cbuffer.data(fs_geom_20_off + 46 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yy_0, g_yz_0_yy_y, g_yz_0_yyy_0, g_z_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyy_0[k] = -g_z_0_yy_0[k] - g_yz_0_yy_0[k] * ab_y + g_yz_0_yy_y[k];
            }

            /// Set up 47-48 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyz_0 = cbuffer.data(fs_geom_20_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyz_0, g_yz_0_yz_0, g_yz_0_yz_y, g_z_0_yz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyz_0[k] = -g_z_0_yz_0[k] - g_yz_0_yz_0[k] * ab_y + g_yz_0_yz_y[k];
            }

            /// Set up 48-49 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yzz_0 = cbuffer.data(fs_geom_20_off + 48 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yzz_0, g_yz_0_zz_0, g_yz_0_zz_y, g_z_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yzz_0[k] = -g_z_0_zz_0[k] - g_yz_0_zz_0[k] * ab_y + g_yz_0_zz_y[k];
            }

            /// Set up 49-50 components of targeted buffer : cbuffer.data(

            auto g_yz_0_zzz_0 = cbuffer.data(fs_geom_20_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zz_0, g_yz_0_zz_0, g_yz_0_zz_z, g_yz_0_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_zzz_0[k] = -g_y_0_zz_0[k] - g_yz_0_zz_0[k] * ab_z + g_yz_0_zz_z[k];
            }

            /// Set up 50-51 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxx_0 = cbuffer.data(fs_geom_20_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xx_0, g_zz_0_xx_x, g_zz_0_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxx_0[k] = -g_zz_0_xx_0[k] * ab_x + g_zz_0_xx_x[k];
            }

            /// Set up 51-52 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxy_0 = cbuffer.data(fs_geom_20_off + 51 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxy_0, g_zz_0_xy_0, g_zz_0_xy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxy_0[k] = -g_zz_0_xy_0[k] * ab_x + g_zz_0_xy_x[k];
            }

            /// Set up 52-53 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxz_0 = cbuffer.data(fs_geom_20_off + 52 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxz_0, g_zz_0_xz_0, g_zz_0_xz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxz_0[k] = -g_zz_0_xz_0[k] * ab_x + g_zz_0_xz_x[k];
            }

            /// Set up 53-54 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyy_0 = cbuffer.data(fs_geom_20_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyy_0, g_zz_0_yy_0, g_zz_0_yy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyy_0[k] = -g_zz_0_yy_0[k] * ab_x + g_zz_0_yy_x[k];
            }

            /// Set up 54-55 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyz_0 = cbuffer.data(fs_geom_20_off + 54 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyz_0, g_zz_0_yz_0, g_zz_0_yz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyz_0[k] = -g_zz_0_yz_0[k] * ab_x + g_zz_0_yz_x[k];
            }

            /// Set up 55-56 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xzz_0 = cbuffer.data(fs_geom_20_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xzz_0, g_zz_0_zz_0, g_zz_0_zz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xzz_0[k] = -g_zz_0_zz_0[k] * ab_x + g_zz_0_zz_x[k];
            }

            /// Set up 56-57 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyy_0 = cbuffer.data(fs_geom_20_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yy_0, g_zz_0_yy_y, g_zz_0_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyy_0[k] = -g_zz_0_yy_0[k] * ab_y + g_zz_0_yy_y[k];
            }

            /// Set up 57-58 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyz_0 = cbuffer.data(fs_geom_20_off + 57 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyz_0, g_zz_0_yz_0, g_zz_0_yz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyz_0[k] = -g_zz_0_yz_0[k] * ab_y + g_zz_0_yz_y[k];
            }

            /// Set up 58-59 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yzz_0 = cbuffer.data(fs_geom_20_off + 58 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yzz_0, g_zz_0_zz_0, g_zz_0_zz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yzz_0[k] = -g_zz_0_zz_0[k] * ab_y + g_zz_0_zz_y[k];
            }

            /// Set up 59-60 components of targeted buffer : cbuffer.data(

            auto g_zz_0_zzz_0 = cbuffer.data(fs_geom_20_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zz_0, g_zz_0_zz_0, g_zz_0_zz_z, g_zz_0_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_zzz_0[k] = -2.0 * g_z_0_zz_0[k] - g_zz_0_zz_0[k] * ab_z + g_zz_0_zz_z[k];
            }
        }
    }
}

} // erirec namespace

