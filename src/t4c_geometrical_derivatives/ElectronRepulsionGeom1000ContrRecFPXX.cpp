#include "ElectronRepulsionGeom1000ContrRecFPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_fpxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_fpxx,
                                            const size_t idx_dpxx,
                                            const size_t idx_geom_10_dpxx,
                                            const size_t idx_geom_10_ddxx,
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
            /// Set up components of auxilary buffer : DPSS

            const auto dp_off = idx_dpxx + i * dcomps + j;

            auto g_xx_x = cbuffer.data(dp_off + 0 * ccomps * dcomps);

            auto g_xx_y = cbuffer.data(dp_off + 1 * ccomps * dcomps);

            auto g_xx_z = cbuffer.data(dp_off + 2 * ccomps * dcomps);

            auto g_xy_x = cbuffer.data(dp_off + 3 * ccomps * dcomps);

            auto g_xy_y = cbuffer.data(dp_off + 4 * ccomps * dcomps);

            auto g_xy_z = cbuffer.data(dp_off + 5 * ccomps * dcomps);

            auto g_xz_x = cbuffer.data(dp_off + 6 * ccomps * dcomps);

            auto g_xz_y = cbuffer.data(dp_off + 7 * ccomps * dcomps);

            auto g_xz_z = cbuffer.data(dp_off + 8 * ccomps * dcomps);

            auto g_yy_x = cbuffer.data(dp_off + 9 * ccomps * dcomps);

            auto g_yy_y = cbuffer.data(dp_off + 10 * ccomps * dcomps);

            auto g_yy_z = cbuffer.data(dp_off + 11 * ccomps * dcomps);

            auto g_yz_x = cbuffer.data(dp_off + 12 * ccomps * dcomps);

            auto g_yz_y = cbuffer.data(dp_off + 13 * ccomps * dcomps);

            auto g_yz_z = cbuffer.data(dp_off + 14 * ccomps * dcomps);

            auto g_zz_x = cbuffer.data(dp_off + 15 * ccomps * dcomps);

            auto g_zz_y = cbuffer.data(dp_off + 16 * ccomps * dcomps);

            auto g_zz_z = cbuffer.data(dp_off + 17 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DPSS

            const auto dp_geom_10_off = idx_geom_10_dpxx + i * dcomps + j;

            auto g_x_0_xx_x = cbuffer.data(dp_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xx_y = cbuffer.data(dp_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xx_z = cbuffer.data(dp_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xy_x = cbuffer.data(dp_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xy_y = cbuffer.data(dp_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xy_z = cbuffer.data(dp_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xz_x = cbuffer.data(dp_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xz_y = cbuffer.data(dp_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xz_z = cbuffer.data(dp_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_yy_x = cbuffer.data(dp_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_yy_y = cbuffer.data(dp_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_yy_z = cbuffer.data(dp_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_yz_x = cbuffer.data(dp_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_yz_y = cbuffer.data(dp_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_yz_z = cbuffer.data(dp_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_zz_x = cbuffer.data(dp_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_zz_y = cbuffer.data(dp_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_zz_z = cbuffer.data(dp_geom_10_off + 17 * ccomps * dcomps);

            auto g_y_0_xx_x = cbuffer.data(dp_geom_10_off + 18 * ccomps * dcomps);

            auto g_y_0_xx_y = cbuffer.data(dp_geom_10_off + 19 * ccomps * dcomps);

            auto g_y_0_xx_z = cbuffer.data(dp_geom_10_off + 20 * ccomps * dcomps);

            auto g_y_0_xy_x = cbuffer.data(dp_geom_10_off + 21 * ccomps * dcomps);

            auto g_y_0_xy_y = cbuffer.data(dp_geom_10_off + 22 * ccomps * dcomps);

            auto g_y_0_xy_z = cbuffer.data(dp_geom_10_off + 23 * ccomps * dcomps);

            auto g_y_0_xz_x = cbuffer.data(dp_geom_10_off + 24 * ccomps * dcomps);

            auto g_y_0_xz_y = cbuffer.data(dp_geom_10_off + 25 * ccomps * dcomps);

            auto g_y_0_xz_z = cbuffer.data(dp_geom_10_off + 26 * ccomps * dcomps);

            auto g_y_0_yy_x = cbuffer.data(dp_geom_10_off + 27 * ccomps * dcomps);

            auto g_y_0_yy_y = cbuffer.data(dp_geom_10_off + 28 * ccomps * dcomps);

            auto g_y_0_yy_z = cbuffer.data(dp_geom_10_off + 29 * ccomps * dcomps);

            auto g_y_0_yz_x = cbuffer.data(dp_geom_10_off + 30 * ccomps * dcomps);

            auto g_y_0_yz_y = cbuffer.data(dp_geom_10_off + 31 * ccomps * dcomps);

            auto g_y_0_yz_z = cbuffer.data(dp_geom_10_off + 32 * ccomps * dcomps);

            auto g_y_0_zz_x = cbuffer.data(dp_geom_10_off + 33 * ccomps * dcomps);

            auto g_y_0_zz_y = cbuffer.data(dp_geom_10_off + 34 * ccomps * dcomps);

            auto g_y_0_zz_z = cbuffer.data(dp_geom_10_off + 35 * ccomps * dcomps);

            auto g_z_0_xx_x = cbuffer.data(dp_geom_10_off + 36 * ccomps * dcomps);

            auto g_z_0_xx_y = cbuffer.data(dp_geom_10_off + 37 * ccomps * dcomps);

            auto g_z_0_xx_z = cbuffer.data(dp_geom_10_off + 38 * ccomps * dcomps);

            auto g_z_0_xy_x = cbuffer.data(dp_geom_10_off + 39 * ccomps * dcomps);

            auto g_z_0_xy_y = cbuffer.data(dp_geom_10_off + 40 * ccomps * dcomps);

            auto g_z_0_xy_z = cbuffer.data(dp_geom_10_off + 41 * ccomps * dcomps);

            auto g_z_0_xz_x = cbuffer.data(dp_geom_10_off + 42 * ccomps * dcomps);

            auto g_z_0_xz_y = cbuffer.data(dp_geom_10_off + 43 * ccomps * dcomps);

            auto g_z_0_xz_z = cbuffer.data(dp_geom_10_off + 44 * ccomps * dcomps);

            auto g_z_0_yy_x = cbuffer.data(dp_geom_10_off + 45 * ccomps * dcomps);

            auto g_z_0_yy_y = cbuffer.data(dp_geom_10_off + 46 * ccomps * dcomps);

            auto g_z_0_yy_z = cbuffer.data(dp_geom_10_off + 47 * ccomps * dcomps);

            auto g_z_0_yz_x = cbuffer.data(dp_geom_10_off + 48 * ccomps * dcomps);

            auto g_z_0_yz_y = cbuffer.data(dp_geom_10_off + 49 * ccomps * dcomps);

            auto g_z_0_yz_z = cbuffer.data(dp_geom_10_off + 50 * ccomps * dcomps);

            auto g_z_0_zz_x = cbuffer.data(dp_geom_10_off + 51 * ccomps * dcomps);

            auto g_z_0_zz_y = cbuffer.data(dp_geom_10_off + 52 * ccomps * dcomps);

            auto g_z_0_zz_z = cbuffer.data(dp_geom_10_off + 53 * ccomps * dcomps);

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

            /// set up bra offset for contr_buffer_fpxx

            const auto fp_geom_10_off = idx_geom_10_fpxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxx_x = cbuffer.data(fp_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxx_y = cbuffer.data(fp_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxx_z = cbuffer.data(fp_geom_10_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xx_x, g_x_0_xx_xx, g_x_0_xx_xy, g_x_0_xx_xz, g_x_0_xx_y, g_x_0_xx_z, g_x_0_xxx_x, g_x_0_xxx_y, g_x_0_xxx_z, g_xx_x, g_xx_y, g_xx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxx_x[k] = -g_xx_x[k] - g_x_0_xx_x[k] * ab_x + g_x_0_xx_xx[k];

                g_x_0_xxx_y[k] = -g_xx_y[k] - g_x_0_xx_y[k] * ab_x + g_x_0_xx_xy[k];

                g_x_0_xxx_z[k] = -g_xx_z[k] - g_x_0_xx_z[k] * ab_x + g_x_0_xx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxy_x = cbuffer.data(fp_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxy_y = cbuffer.data(fp_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxy_z = cbuffer.data(fp_geom_10_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xx_x, g_x_0_xx_xy, g_x_0_xx_y, g_x_0_xx_yy, g_x_0_xx_yz, g_x_0_xx_z, g_x_0_xxy_x, g_x_0_xxy_y, g_x_0_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxy_x[k] = -g_x_0_xx_x[k] * ab_y + g_x_0_xx_xy[k];

                g_x_0_xxy_y[k] = -g_x_0_xx_y[k] * ab_y + g_x_0_xx_yy[k];

                g_x_0_xxy_z[k] = -g_x_0_xx_z[k] * ab_y + g_x_0_xx_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxz_x = cbuffer.data(fp_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxz_y = cbuffer.data(fp_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxz_z = cbuffer.data(fp_geom_10_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xx_x, g_x_0_xx_xz, g_x_0_xx_y, g_x_0_xx_yz, g_x_0_xx_z, g_x_0_xx_zz, g_x_0_xxz_x, g_x_0_xxz_y, g_x_0_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxz_x[k] = -g_x_0_xx_x[k] * ab_z + g_x_0_xx_xz[k];

                g_x_0_xxz_y[k] = -g_x_0_xx_y[k] * ab_z + g_x_0_xx_yz[k];

                g_x_0_xxz_z[k] = -g_x_0_xx_z[k] * ab_z + g_x_0_xx_zz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyy_x = cbuffer.data(fp_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xyy_y = cbuffer.data(fp_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xyy_z = cbuffer.data(fp_geom_10_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xy_x, g_x_0_xy_xy, g_x_0_xy_y, g_x_0_xy_yy, g_x_0_xy_yz, g_x_0_xy_z, g_x_0_xyy_x, g_x_0_xyy_y, g_x_0_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyy_x[k] = -g_x_0_xy_x[k] * ab_y + g_x_0_xy_xy[k];

                g_x_0_xyy_y[k] = -g_x_0_xy_y[k] * ab_y + g_x_0_xy_yy[k];

                g_x_0_xyy_z[k] = -g_x_0_xy_z[k] * ab_y + g_x_0_xy_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyz_x = cbuffer.data(fp_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xyz_y = cbuffer.data(fp_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xyz_z = cbuffer.data(fp_geom_10_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyz_x, g_x_0_xyz_y, g_x_0_xyz_z, g_x_0_xz_x, g_x_0_xz_xy, g_x_0_xz_y, g_x_0_xz_yy, g_x_0_xz_yz, g_x_0_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyz_x[k] = -g_x_0_xz_x[k] * ab_y + g_x_0_xz_xy[k];

                g_x_0_xyz_y[k] = -g_x_0_xz_y[k] * ab_y + g_x_0_xz_yy[k];

                g_x_0_xyz_z[k] = -g_x_0_xz_z[k] * ab_y + g_x_0_xz_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzz_x = cbuffer.data(fp_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xzz_y = cbuffer.data(fp_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xzz_z = cbuffer.data(fp_geom_10_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xz_x, g_x_0_xz_xz, g_x_0_xz_y, g_x_0_xz_yz, g_x_0_xz_z, g_x_0_xz_zz, g_x_0_xzz_x, g_x_0_xzz_y, g_x_0_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzz_x[k] = -g_x_0_xz_x[k] * ab_z + g_x_0_xz_xz[k];

                g_x_0_xzz_y[k] = -g_x_0_xz_y[k] * ab_z + g_x_0_xz_yz[k];

                g_x_0_xzz_z[k] = -g_x_0_xz_z[k] * ab_z + g_x_0_xz_zz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyy_x = cbuffer.data(fp_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_yyy_y = cbuffer.data(fp_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_yyy_z = cbuffer.data(fp_geom_10_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yy_x, g_x_0_yy_xy, g_x_0_yy_y, g_x_0_yy_yy, g_x_0_yy_yz, g_x_0_yy_z, g_x_0_yyy_x, g_x_0_yyy_y, g_x_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyy_x[k] = -g_x_0_yy_x[k] * ab_y + g_x_0_yy_xy[k];

                g_x_0_yyy_y[k] = -g_x_0_yy_y[k] * ab_y + g_x_0_yy_yy[k];

                g_x_0_yyy_z[k] = -g_x_0_yy_z[k] * ab_y + g_x_0_yy_yz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyz_x = cbuffer.data(fp_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_yyz_y = cbuffer.data(fp_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_yyz_z = cbuffer.data(fp_geom_10_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyz_x, g_x_0_yyz_y, g_x_0_yyz_z, g_x_0_yz_x, g_x_0_yz_xy, g_x_0_yz_y, g_x_0_yz_yy, g_x_0_yz_yz, g_x_0_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyz_x[k] = -g_x_0_yz_x[k] * ab_y + g_x_0_yz_xy[k];

                g_x_0_yyz_y[k] = -g_x_0_yz_y[k] * ab_y + g_x_0_yz_yy[k];

                g_x_0_yyz_z[k] = -g_x_0_yz_z[k] * ab_y + g_x_0_yz_yz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzz_x = cbuffer.data(fp_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_yzz_y = cbuffer.data(fp_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_yzz_z = cbuffer.data(fp_geom_10_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yzz_x, g_x_0_yzz_y, g_x_0_yzz_z, g_x_0_zz_x, g_x_0_zz_xy, g_x_0_zz_y, g_x_0_zz_yy, g_x_0_zz_yz, g_x_0_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzz_x[k] = -g_x_0_zz_x[k] * ab_y + g_x_0_zz_xy[k];

                g_x_0_yzz_y[k] = -g_x_0_zz_y[k] * ab_y + g_x_0_zz_yy[k];

                g_x_0_yzz_z[k] = -g_x_0_zz_z[k] * ab_y + g_x_0_zz_yz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzz_x = cbuffer.data(fp_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_zzz_y = cbuffer.data(fp_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_zzz_z = cbuffer.data(fp_geom_10_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zz_x, g_x_0_zz_xz, g_x_0_zz_y, g_x_0_zz_yz, g_x_0_zz_z, g_x_0_zz_zz, g_x_0_zzz_x, g_x_0_zzz_y, g_x_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzz_x[k] = -g_x_0_zz_x[k] * ab_z + g_x_0_zz_xz[k];

                g_x_0_zzz_y[k] = -g_x_0_zz_y[k] * ab_z + g_x_0_zz_yz[k];

                g_x_0_zzz_z[k] = -g_x_0_zz_z[k] * ab_z + g_x_0_zz_zz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxx_x = cbuffer.data(fp_geom_10_off + 30 * ccomps * dcomps);

            auto g_y_0_xxx_y = cbuffer.data(fp_geom_10_off + 31 * ccomps * dcomps);

            auto g_y_0_xxx_z = cbuffer.data(fp_geom_10_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xx_x, g_y_0_xx_xx, g_y_0_xx_xy, g_y_0_xx_xz, g_y_0_xx_y, g_y_0_xx_z, g_y_0_xxx_x, g_y_0_xxx_y, g_y_0_xxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxx_x[k] = -g_y_0_xx_x[k] * ab_x + g_y_0_xx_xx[k];

                g_y_0_xxx_y[k] = -g_y_0_xx_y[k] * ab_x + g_y_0_xx_xy[k];

                g_y_0_xxx_z[k] = -g_y_0_xx_z[k] * ab_x + g_y_0_xx_xz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxy_x = cbuffer.data(fp_geom_10_off + 33 * ccomps * dcomps);

            auto g_y_0_xxy_y = cbuffer.data(fp_geom_10_off + 34 * ccomps * dcomps);

            auto g_y_0_xxy_z = cbuffer.data(fp_geom_10_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxy_x, g_y_0_xxy_y, g_y_0_xxy_z, g_y_0_xy_x, g_y_0_xy_xx, g_y_0_xy_xy, g_y_0_xy_xz, g_y_0_xy_y, g_y_0_xy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxy_x[k] = -g_y_0_xy_x[k] * ab_x + g_y_0_xy_xx[k];

                g_y_0_xxy_y[k] = -g_y_0_xy_y[k] * ab_x + g_y_0_xy_xy[k];

                g_y_0_xxy_z[k] = -g_y_0_xy_z[k] * ab_x + g_y_0_xy_xz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxz_x = cbuffer.data(fp_geom_10_off + 36 * ccomps * dcomps);

            auto g_y_0_xxz_y = cbuffer.data(fp_geom_10_off + 37 * ccomps * dcomps);

            auto g_y_0_xxz_z = cbuffer.data(fp_geom_10_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxz_x, g_y_0_xxz_y, g_y_0_xxz_z, g_y_0_xz_x, g_y_0_xz_xx, g_y_0_xz_xy, g_y_0_xz_xz, g_y_0_xz_y, g_y_0_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxz_x[k] = -g_y_0_xz_x[k] * ab_x + g_y_0_xz_xx[k];

                g_y_0_xxz_y[k] = -g_y_0_xz_y[k] * ab_x + g_y_0_xz_xy[k];

                g_y_0_xxz_z[k] = -g_y_0_xz_z[k] * ab_x + g_y_0_xz_xz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyy_x = cbuffer.data(fp_geom_10_off + 39 * ccomps * dcomps);

            auto g_y_0_xyy_y = cbuffer.data(fp_geom_10_off + 40 * ccomps * dcomps);

            auto g_y_0_xyy_z = cbuffer.data(fp_geom_10_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyy_x, g_y_0_xyy_y, g_y_0_xyy_z, g_y_0_yy_x, g_y_0_yy_xx, g_y_0_yy_xy, g_y_0_yy_xz, g_y_0_yy_y, g_y_0_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyy_x[k] = -g_y_0_yy_x[k] * ab_x + g_y_0_yy_xx[k];

                g_y_0_xyy_y[k] = -g_y_0_yy_y[k] * ab_x + g_y_0_yy_xy[k];

                g_y_0_xyy_z[k] = -g_y_0_yy_z[k] * ab_x + g_y_0_yy_xz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyz_x = cbuffer.data(fp_geom_10_off + 42 * ccomps * dcomps);

            auto g_y_0_xyz_y = cbuffer.data(fp_geom_10_off + 43 * ccomps * dcomps);

            auto g_y_0_xyz_z = cbuffer.data(fp_geom_10_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyz_x, g_y_0_xyz_y, g_y_0_xyz_z, g_y_0_yz_x, g_y_0_yz_xx, g_y_0_yz_xy, g_y_0_yz_xz, g_y_0_yz_y, g_y_0_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyz_x[k] = -g_y_0_yz_x[k] * ab_x + g_y_0_yz_xx[k];

                g_y_0_xyz_y[k] = -g_y_0_yz_y[k] * ab_x + g_y_0_yz_xy[k];

                g_y_0_xyz_z[k] = -g_y_0_yz_z[k] * ab_x + g_y_0_yz_xz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzz_x = cbuffer.data(fp_geom_10_off + 45 * ccomps * dcomps);

            auto g_y_0_xzz_y = cbuffer.data(fp_geom_10_off + 46 * ccomps * dcomps);

            auto g_y_0_xzz_z = cbuffer.data(fp_geom_10_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xzz_x, g_y_0_xzz_y, g_y_0_xzz_z, g_y_0_zz_x, g_y_0_zz_xx, g_y_0_zz_xy, g_y_0_zz_xz, g_y_0_zz_y, g_y_0_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzz_x[k] = -g_y_0_zz_x[k] * ab_x + g_y_0_zz_xx[k];

                g_y_0_xzz_y[k] = -g_y_0_zz_y[k] * ab_x + g_y_0_zz_xy[k];

                g_y_0_xzz_z[k] = -g_y_0_zz_z[k] * ab_x + g_y_0_zz_xz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyy_x = cbuffer.data(fp_geom_10_off + 48 * ccomps * dcomps);

            auto g_y_0_yyy_y = cbuffer.data(fp_geom_10_off + 49 * ccomps * dcomps);

            auto g_y_0_yyy_z = cbuffer.data(fp_geom_10_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yy_x, g_y_0_yy_xy, g_y_0_yy_y, g_y_0_yy_yy, g_y_0_yy_yz, g_y_0_yy_z, g_y_0_yyy_x, g_y_0_yyy_y, g_y_0_yyy_z, g_yy_x, g_yy_y, g_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyy_x[k] = -g_yy_x[k] - g_y_0_yy_x[k] * ab_y + g_y_0_yy_xy[k];

                g_y_0_yyy_y[k] = -g_yy_y[k] - g_y_0_yy_y[k] * ab_y + g_y_0_yy_yy[k];

                g_y_0_yyy_z[k] = -g_yy_z[k] - g_y_0_yy_z[k] * ab_y + g_y_0_yy_yz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyz_x = cbuffer.data(fp_geom_10_off + 51 * ccomps * dcomps);

            auto g_y_0_yyz_y = cbuffer.data(fp_geom_10_off + 52 * ccomps * dcomps);

            auto g_y_0_yyz_z = cbuffer.data(fp_geom_10_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yy_x, g_y_0_yy_xz, g_y_0_yy_y, g_y_0_yy_yz, g_y_0_yy_z, g_y_0_yy_zz, g_y_0_yyz_x, g_y_0_yyz_y, g_y_0_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyz_x[k] = -g_y_0_yy_x[k] * ab_z + g_y_0_yy_xz[k];

                g_y_0_yyz_y[k] = -g_y_0_yy_y[k] * ab_z + g_y_0_yy_yz[k];

                g_y_0_yyz_z[k] = -g_y_0_yy_z[k] * ab_z + g_y_0_yy_zz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzz_x = cbuffer.data(fp_geom_10_off + 54 * ccomps * dcomps);

            auto g_y_0_yzz_y = cbuffer.data(fp_geom_10_off + 55 * ccomps * dcomps);

            auto g_y_0_yzz_z = cbuffer.data(fp_geom_10_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yz_x, g_y_0_yz_xz, g_y_0_yz_y, g_y_0_yz_yz, g_y_0_yz_z, g_y_0_yz_zz, g_y_0_yzz_x, g_y_0_yzz_y, g_y_0_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzz_x[k] = -g_y_0_yz_x[k] * ab_z + g_y_0_yz_xz[k];

                g_y_0_yzz_y[k] = -g_y_0_yz_y[k] * ab_z + g_y_0_yz_yz[k];

                g_y_0_yzz_z[k] = -g_y_0_yz_z[k] * ab_z + g_y_0_yz_zz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzz_x = cbuffer.data(fp_geom_10_off + 57 * ccomps * dcomps);

            auto g_y_0_zzz_y = cbuffer.data(fp_geom_10_off + 58 * ccomps * dcomps);

            auto g_y_0_zzz_z = cbuffer.data(fp_geom_10_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zz_x, g_y_0_zz_xz, g_y_0_zz_y, g_y_0_zz_yz, g_y_0_zz_z, g_y_0_zz_zz, g_y_0_zzz_x, g_y_0_zzz_y, g_y_0_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzz_x[k] = -g_y_0_zz_x[k] * ab_z + g_y_0_zz_xz[k];

                g_y_0_zzz_y[k] = -g_y_0_zz_y[k] * ab_z + g_y_0_zz_yz[k];

                g_y_0_zzz_z[k] = -g_y_0_zz_z[k] * ab_z + g_y_0_zz_zz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxx_x = cbuffer.data(fp_geom_10_off + 60 * ccomps * dcomps);

            auto g_z_0_xxx_y = cbuffer.data(fp_geom_10_off + 61 * ccomps * dcomps);

            auto g_z_0_xxx_z = cbuffer.data(fp_geom_10_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xx_x, g_z_0_xx_xx, g_z_0_xx_xy, g_z_0_xx_xz, g_z_0_xx_y, g_z_0_xx_z, g_z_0_xxx_x, g_z_0_xxx_y, g_z_0_xxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxx_x[k] = -g_z_0_xx_x[k] * ab_x + g_z_0_xx_xx[k];

                g_z_0_xxx_y[k] = -g_z_0_xx_y[k] * ab_x + g_z_0_xx_xy[k];

                g_z_0_xxx_z[k] = -g_z_0_xx_z[k] * ab_x + g_z_0_xx_xz[k];
            }

            /// Set up 63-66 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxy_x = cbuffer.data(fp_geom_10_off + 63 * ccomps * dcomps);

            auto g_z_0_xxy_y = cbuffer.data(fp_geom_10_off + 64 * ccomps * dcomps);

            auto g_z_0_xxy_z = cbuffer.data(fp_geom_10_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxy_x, g_z_0_xxy_y, g_z_0_xxy_z, g_z_0_xy_x, g_z_0_xy_xx, g_z_0_xy_xy, g_z_0_xy_xz, g_z_0_xy_y, g_z_0_xy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxy_x[k] = -g_z_0_xy_x[k] * ab_x + g_z_0_xy_xx[k];

                g_z_0_xxy_y[k] = -g_z_0_xy_y[k] * ab_x + g_z_0_xy_xy[k];

                g_z_0_xxy_z[k] = -g_z_0_xy_z[k] * ab_x + g_z_0_xy_xz[k];
            }

            /// Set up 66-69 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxz_x = cbuffer.data(fp_geom_10_off + 66 * ccomps * dcomps);

            auto g_z_0_xxz_y = cbuffer.data(fp_geom_10_off + 67 * ccomps * dcomps);

            auto g_z_0_xxz_z = cbuffer.data(fp_geom_10_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxz_x, g_z_0_xxz_y, g_z_0_xxz_z, g_z_0_xz_x, g_z_0_xz_xx, g_z_0_xz_xy, g_z_0_xz_xz, g_z_0_xz_y, g_z_0_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxz_x[k] = -g_z_0_xz_x[k] * ab_x + g_z_0_xz_xx[k];

                g_z_0_xxz_y[k] = -g_z_0_xz_y[k] * ab_x + g_z_0_xz_xy[k];

                g_z_0_xxz_z[k] = -g_z_0_xz_z[k] * ab_x + g_z_0_xz_xz[k];
            }

            /// Set up 69-72 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyy_x = cbuffer.data(fp_geom_10_off + 69 * ccomps * dcomps);

            auto g_z_0_xyy_y = cbuffer.data(fp_geom_10_off + 70 * ccomps * dcomps);

            auto g_z_0_xyy_z = cbuffer.data(fp_geom_10_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyy_x, g_z_0_xyy_y, g_z_0_xyy_z, g_z_0_yy_x, g_z_0_yy_xx, g_z_0_yy_xy, g_z_0_yy_xz, g_z_0_yy_y, g_z_0_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyy_x[k] = -g_z_0_yy_x[k] * ab_x + g_z_0_yy_xx[k];

                g_z_0_xyy_y[k] = -g_z_0_yy_y[k] * ab_x + g_z_0_yy_xy[k];

                g_z_0_xyy_z[k] = -g_z_0_yy_z[k] * ab_x + g_z_0_yy_xz[k];
            }

            /// Set up 72-75 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyz_x = cbuffer.data(fp_geom_10_off + 72 * ccomps * dcomps);

            auto g_z_0_xyz_y = cbuffer.data(fp_geom_10_off + 73 * ccomps * dcomps);

            auto g_z_0_xyz_z = cbuffer.data(fp_geom_10_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyz_x, g_z_0_xyz_y, g_z_0_xyz_z, g_z_0_yz_x, g_z_0_yz_xx, g_z_0_yz_xy, g_z_0_yz_xz, g_z_0_yz_y, g_z_0_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyz_x[k] = -g_z_0_yz_x[k] * ab_x + g_z_0_yz_xx[k];

                g_z_0_xyz_y[k] = -g_z_0_yz_y[k] * ab_x + g_z_0_yz_xy[k];

                g_z_0_xyz_z[k] = -g_z_0_yz_z[k] * ab_x + g_z_0_yz_xz[k];
            }

            /// Set up 75-78 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzz_x = cbuffer.data(fp_geom_10_off + 75 * ccomps * dcomps);

            auto g_z_0_xzz_y = cbuffer.data(fp_geom_10_off + 76 * ccomps * dcomps);

            auto g_z_0_xzz_z = cbuffer.data(fp_geom_10_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xzz_x, g_z_0_xzz_y, g_z_0_xzz_z, g_z_0_zz_x, g_z_0_zz_xx, g_z_0_zz_xy, g_z_0_zz_xz, g_z_0_zz_y, g_z_0_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzz_x[k] = -g_z_0_zz_x[k] * ab_x + g_z_0_zz_xx[k];

                g_z_0_xzz_y[k] = -g_z_0_zz_y[k] * ab_x + g_z_0_zz_xy[k];

                g_z_0_xzz_z[k] = -g_z_0_zz_z[k] * ab_x + g_z_0_zz_xz[k];
            }

            /// Set up 78-81 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyy_x = cbuffer.data(fp_geom_10_off + 78 * ccomps * dcomps);

            auto g_z_0_yyy_y = cbuffer.data(fp_geom_10_off + 79 * ccomps * dcomps);

            auto g_z_0_yyy_z = cbuffer.data(fp_geom_10_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yy_x, g_z_0_yy_xy, g_z_0_yy_y, g_z_0_yy_yy, g_z_0_yy_yz, g_z_0_yy_z, g_z_0_yyy_x, g_z_0_yyy_y, g_z_0_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyy_x[k] = -g_z_0_yy_x[k] * ab_y + g_z_0_yy_xy[k];

                g_z_0_yyy_y[k] = -g_z_0_yy_y[k] * ab_y + g_z_0_yy_yy[k];

                g_z_0_yyy_z[k] = -g_z_0_yy_z[k] * ab_y + g_z_0_yy_yz[k];
            }

            /// Set up 81-84 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyz_x = cbuffer.data(fp_geom_10_off + 81 * ccomps * dcomps);

            auto g_z_0_yyz_y = cbuffer.data(fp_geom_10_off + 82 * ccomps * dcomps);

            auto g_z_0_yyz_z = cbuffer.data(fp_geom_10_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyz_x, g_z_0_yyz_y, g_z_0_yyz_z, g_z_0_yz_x, g_z_0_yz_xy, g_z_0_yz_y, g_z_0_yz_yy, g_z_0_yz_yz, g_z_0_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyz_x[k] = -g_z_0_yz_x[k] * ab_y + g_z_0_yz_xy[k];

                g_z_0_yyz_y[k] = -g_z_0_yz_y[k] * ab_y + g_z_0_yz_yy[k];

                g_z_0_yyz_z[k] = -g_z_0_yz_z[k] * ab_y + g_z_0_yz_yz[k];
            }

            /// Set up 84-87 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzz_x = cbuffer.data(fp_geom_10_off + 84 * ccomps * dcomps);

            auto g_z_0_yzz_y = cbuffer.data(fp_geom_10_off + 85 * ccomps * dcomps);

            auto g_z_0_yzz_z = cbuffer.data(fp_geom_10_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yzz_x, g_z_0_yzz_y, g_z_0_yzz_z, g_z_0_zz_x, g_z_0_zz_xy, g_z_0_zz_y, g_z_0_zz_yy, g_z_0_zz_yz, g_z_0_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzz_x[k] = -g_z_0_zz_x[k] * ab_y + g_z_0_zz_xy[k];

                g_z_0_yzz_y[k] = -g_z_0_zz_y[k] * ab_y + g_z_0_zz_yy[k];

                g_z_0_yzz_z[k] = -g_z_0_zz_z[k] * ab_y + g_z_0_zz_yz[k];
            }

            /// Set up 87-90 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzz_x = cbuffer.data(fp_geom_10_off + 87 * ccomps * dcomps);

            auto g_z_0_zzz_y = cbuffer.data(fp_geom_10_off + 88 * ccomps * dcomps);

            auto g_z_0_zzz_z = cbuffer.data(fp_geom_10_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zz_x, g_z_0_zz_xz, g_z_0_zz_y, g_z_0_zz_yz, g_z_0_zz_z, g_z_0_zz_zz, g_z_0_zzz_x, g_z_0_zzz_y, g_z_0_zzz_z, g_zz_x, g_zz_y, g_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzz_x[k] = -g_zz_x[k] - g_z_0_zz_x[k] * ab_z + g_z_0_zz_xz[k];

                g_z_0_zzz_y[k] = -g_zz_y[k] - g_z_0_zz_y[k] * ab_z + g_z_0_zz_yz[k];

                g_z_0_zzz_z[k] = -g_zz_z[k] - g_z_0_zz_z[k] * ab_z + g_z_0_zz_zz[k];
            }
        }
    }
}

} // erirec namespace

