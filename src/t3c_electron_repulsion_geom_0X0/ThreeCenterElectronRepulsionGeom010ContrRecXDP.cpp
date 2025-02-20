#include "ThreeCenterElectronRepulsionGeom010ContrRecXDP.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xdp(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xdp,
                                        const size_t idx_xpp,
                                        const size_t idx_geom_10_xpp,
                                        const size_t idx_geom_10_xpd,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        /// Set up components of auxilary buffer : SPP

        const auto pp_off = idx_xpp + i * 9;

        auto g_x_x = cbuffer.data(pp_off + 0);

        auto g_x_y = cbuffer.data(pp_off + 1);

        auto g_x_z = cbuffer.data(pp_off + 2);

        auto g_y_x = cbuffer.data(pp_off + 3);

        auto g_y_y = cbuffer.data(pp_off + 4);

        auto g_y_z = cbuffer.data(pp_off + 5);

        auto g_z_x = cbuffer.data(pp_off + 6);

        auto g_z_y = cbuffer.data(pp_off + 7);

        auto g_z_z = cbuffer.data(pp_off + 8);

        /// Set up components of auxilary buffer : SPP

        const auto pp_geom_10_off = idx_geom_10_xpp + i * 9;

        auto g_x_0_x_x = cbuffer.data(pp_geom_10_off + 0 * acomps + 0);

        auto g_x_0_x_y = cbuffer.data(pp_geom_10_off + 0 * acomps + 1);

        auto g_x_0_x_z = cbuffer.data(pp_geom_10_off + 0 * acomps + 2);

        auto g_x_0_y_x = cbuffer.data(pp_geom_10_off + 0 * acomps + 3);

        auto g_x_0_y_y = cbuffer.data(pp_geom_10_off + 0 * acomps + 4);

        auto g_x_0_y_z = cbuffer.data(pp_geom_10_off + 0 * acomps + 5);

        auto g_x_0_z_x = cbuffer.data(pp_geom_10_off + 0 * acomps + 6);

        auto g_x_0_z_y = cbuffer.data(pp_geom_10_off + 0 * acomps + 7);

        auto g_x_0_z_z = cbuffer.data(pp_geom_10_off + 0 * acomps + 8);

        auto g_y_0_x_x = cbuffer.data(pp_geom_10_off + 9 * acomps + 0);

        auto g_y_0_x_y = cbuffer.data(pp_geom_10_off + 9 * acomps + 1);

        auto g_y_0_x_z = cbuffer.data(pp_geom_10_off + 9 * acomps + 2);

        auto g_y_0_y_x = cbuffer.data(pp_geom_10_off + 9 * acomps + 3);

        auto g_y_0_y_y = cbuffer.data(pp_geom_10_off + 9 * acomps + 4);

        auto g_y_0_y_z = cbuffer.data(pp_geom_10_off + 9 * acomps + 5);

        auto g_y_0_z_x = cbuffer.data(pp_geom_10_off + 9 * acomps + 6);

        auto g_y_0_z_y = cbuffer.data(pp_geom_10_off + 9 * acomps + 7);

        auto g_y_0_z_z = cbuffer.data(pp_geom_10_off + 9 * acomps + 8);

        auto g_z_0_x_x = cbuffer.data(pp_geom_10_off + 18 * acomps + 0);

        auto g_z_0_x_y = cbuffer.data(pp_geom_10_off + 18 * acomps + 1);

        auto g_z_0_x_z = cbuffer.data(pp_geom_10_off + 18 * acomps + 2);

        auto g_z_0_y_x = cbuffer.data(pp_geom_10_off + 18 * acomps + 3);

        auto g_z_0_y_y = cbuffer.data(pp_geom_10_off + 18 * acomps + 4);

        auto g_z_0_y_z = cbuffer.data(pp_geom_10_off + 18 * acomps + 5);

        auto g_z_0_z_x = cbuffer.data(pp_geom_10_off + 18 * acomps + 6);

        auto g_z_0_z_y = cbuffer.data(pp_geom_10_off + 18 * acomps + 7);

        auto g_z_0_z_z = cbuffer.data(pp_geom_10_off + 18 * acomps + 8);

        /// Set up components of auxilary buffer : SPD

        const auto pd_geom_10_off = idx_geom_10_xpd + i * 18;

        auto g_x_0_x_xx = cbuffer.data(pd_geom_10_off + 0 * acomps + 0);

        auto g_x_0_x_xy = cbuffer.data(pd_geom_10_off + 0 * acomps + 1);

        auto g_x_0_x_xz = cbuffer.data(pd_geom_10_off + 0 * acomps + 2);

        auto g_x_0_x_yy = cbuffer.data(pd_geom_10_off + 0 * acomps + 3);

        auto g_x_0_x_yz = cbuffer.data(pd_geom_10_off + 0 * acomps + 4);

        auto g_x_0_x_zz = cbuffer.data(pd_geom_10_off + 0 * acomps + 5);

        auto g_x_0_y_xx = cbuffer.data(pd_geom_10_off + 0 * acomps + 6);

        auto g_x_0_y_xy = cbuffer.data(pd_geom_10_off + 0 * acomps + 7);

        auto g_x_0_y_xz = cbuffer.data(pd_geom_10_off + 0 * acomps + 8);

        auto g_x_0_y_yy = cbuffer.data(pd_geom_10_off + 0 * acomps + 9);

        auto g_x_0_y_yz = cbuffer.data(pd_geom_10_off + 0 * acomps + 10);

        auto g_x_0_y_zz = cbuffer.data(pd_geom_10_off + 0 * acomps + 11);

        auto g_x_0_z_xx = cbuffer.data(pd_geom_10_off + 0 * acomps + 12);

        auto g_x_0_z_xy = cbuffer.data(pd_geom_10_off + 0 * acomps + 13);

        auto g_x_0_z_xz = cbuffer.data(pd_geom_10_off + 0 * acomps + 14);

        auto g_x_0_z_yy = cbuffer.data(pd_geom_10_off + 0 * acomps + 15);

        auto g_x_0_z_yz = cbuffer.data(pd_geom_10_off + 0 * acomps + 16);

        auto g_x_0_z_zz = cbuffer.data(pd_geom_10_off + 0 * acomps + 17);

        auto g_y_0_x_xx = cbuffer.data(pd_geom_10_off + 18 * acomps + 0);

        auto g_y_0_x_xy = cbuffer.data(pd_geom_10_off + 18 * acomps + 1);

        auto g_y_0_x_xz = cbuffer.data(pd_geom_10_off + 18 * acomps + 2);

        auto g_y_0_x_yy = cbuffer.data(pd_geom_10_off + 18 * acomps + 3);

        auto g_y_0_x_yz = cbuffer.data(pd_geom_10_off + 18 * acomps + 4);

        auto g_y_0_x_zz = cbuffer.data(pd_geom_10_off + 18 * acomps + 5);

        auto g_y_0_y_xx = cbuffer.data(pd_geom_10_off + 18 * acomps + 6);

        auto g_y_0_y_xy = cbuffer.data(pd_geom_10_off + 18 * acomps + 7);

        auto g_y_0_y_xz = cbuffer.data(pd_geom_10_off + 18 * acomps + 8);

        auto g_y_0_y_yy = cbuffer.data(pd_geom_10_off + 18 * acomps + 9);

        auto g_y_0_y_yz = cbuffer.data(pd_geom_10_off + 18 * acomps + 10);

        auto g_y_0_y_zz = cbuffer.data(pd_geom_10_off + 18 * acomps + 11);

        auto g_y_0_z_xx = cbuffer.data(pd_geom_10_off + 18 * acomps + 12);

        auto g_y_0_z_xy = cbuffer.data(pd_geom_10_off + 18 * acomps + 13);

        auto g_y_0_z_xz = cbuffer.data(pd_geom_10_off + 18 * acomps + 14);

        auto g_y_0_z_yy = cbuffer.data(pd_geom_10_off + 18 * acomps + 15);

        auto g_y_0_z_yz = cbuffer.data(pd_geom_10_off + 18 * acomps + 16);

        auto g_y_0_z_zz = cbuffer.data(pd_geom_10_off + 18 * acomps + 17);

        auto g_z_0_x_xx = cbuffer.data(pd_geom_10_off + 36 * acomps + 0);

        auto g_z_0_x_xy = cbuffer.data(pd_geom_10_off + 36 * acomps + 1);

        auto g_z_0_x_xz = cbuffer.data(pd_geom_10_off + 36 * acomps + 2);

        auto g_z_0_x_yy = cbuffer.data(pd_geom_10_off + 36 * acomps + 3);

        auto g_z_0_x_yz = cbuffer.data(pd_geom_10_off + 36 * acomps + 4);

        auto g_z_0_x_zz = cbuffer.data(pd_geom_10_off + 36 * acomps + 5);

        auto g_z_0_y_xx = cbuffer.data(pd_geom_10_off + 36 * acomps + 6);

        auto g_z_0_y_xy = cbuffer.data(pd_geom_10_off + 36 * acomps + 7);

        auto g_z_0_y_xz = cbuffer.data(pd_geom_10_off + 36 * acomps + 8);

        auto g_z_0_y_yy = cbuffer.data(pd_geom_10_off + 36 * acomps + 9);

        auto g_z_0_y_yz = cbuffer.data(pd_geom_10_off + 36 * acomps + 10);

        auto g_z_0_y_zz = cbuffer.data(pd_geom_10_off + 36 * acomps + 11);

        auto g_z_0_z_xx = cbuffer.data(pd_geom_10_off + 36 * acomps + 12);

        auto g_z_0_z_xy = cbuffer.data(pd_geom_10_off + 36 * acomps + 13);

        auto g_z_0_z_xz = cbuffer.data(pd_geom_10_off + 36 * acomps + 14);

        auto g_z_0_z_yy = cbuffer.data(pd_geom_10_off + 36 * acomps + 15);

        auto g_z_0_z_yz = cbuffer.data(pd_geom_10_off + 36 * acomps + 16);

        auto g_z_0_z_zz = cbuffer.data(pd_geom_10_off + 36 * acomps + 17);

        /// set up bra offset for contr_buffer_xxdp

        const auto dp_geom_10_off = idx_geom_10_xdp + i * 18;

        /// Set up 0-3 components of targeted buffer : cbuffer.data(

        auto g_x_0_xx_x = cbuffer.data(dp_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_xx_y = cbuffer.data(dp_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_xx_z = cbuffer.data(dp_geom_10_off + 0 * acomps  + 2);

        #pragma omp simd aligned(cd_x, g_x_0_x_x, g_x_0_x_xx, g_x_0_x_xy, g_x_0_x_xz, g_x_0_x_y, g_x_0_x_z, g_x_0_xx_x, g_x_0_xx_y, g_x_0_xx_z, g_x_x, g_x_y, g_x_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xx_x[k] = -g_x_x[k] - g_x_0_x_x[k] * cd_x[k] + g_x_0_x_xx[k];

            g_x_0_xx_y[k] = -g_x_y[k] - g_x_0_x_y[k] * cd_x[k] + g_x_0_x_xy[k];

            g_x_0_xx_z[k] = -g_x_z[k] - g_x_0_x_z[k] * cd_x[k] + g_x_0_x_xz[k];
        }

        /// Set up 3-6 components of targeted buffer : cbuffer.data(

        auto g_x_0_xy_x = cbuffer.data(dp_geom_10_off + 0 * acomps  + 3);

        auto g_x_0_xy_y = cbuffer.data(dp_geom_10_off + 0 * acomps  + 4);

        auto g_x_0_xy_z = cbuffer.data(dp_geom_10_off + 0 * acomps  + 5);

        #pragma omp simd aligned(cd_y, g_x_0_x_x, g_x_0_x_xy, g_x_0_x_y, g_x_0_x_yy, g_x_0_x_yz, g_x_0_x_z, g_x_0_xy_x, g_x_0_xy_y, g_x_0_xy_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xy_x[k] = -g_x_0_x_x[k] * cd_y[k] + g_x_0_x_xy[k];

            g_x_0_xy_y[k] = -g_x_0_x_y[k] * cd_y[k] + g_x_0_x_yy[k];

            g_x_0_xy_z[k] = -g_x_0_x_z[k] * cd_y[k] + g_x_0_x_yz[k];
        }

        /// Set up 6-9 components of targeted buffer : cbuffer.data(

        auto g_x_0_xz_x = cbuffer.data(dp_geom_10_off + 0 * acomps  + 6);

        auto g_x_0_xz_y = cbuffer.data(dp_geom_10_off + 0 * acomps  + 7);

        auto g_x_0_xz_z = cbuffer.data(dp_geom_10_off + 0 * acomps  + 8);

        #pragma omp simd aligned(cd_z, g_x_0_x_x, g_x_0_x_xz, g_x_0_x_y, g_x_0_x_yz, g_x_0_x_z, g_x_0_x_zz, g_x_0_xz_x, g_x_0_xz_y, g_x_0_xz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xz_x[k] = -g_x_0_x_x[k] * cd_z[k] + g_x_0_x_xz[k];

            g_x_0_xz_y[k] = -g_x_0_x_y[k] * cd_z[k] + g_x_0_x_yz[k];

            g_x_0_xz_z[k] = -g_x_0_x_z[k] * cd_z[k] + g_x_0_x_zz[k];
        }

        /// Set up 9-12 components of targeted buffer : cbuffer.data(

        auto g_x_0_yy_x = cbuffer.data(dp_geom_10_off + 0 * acomps  + 9);

        auto g_x_0_yy_y = cbuffer.data(dp_geom_10_off + 0 * acomps  + 10);

        auto g_x_0_yy_z = cbuffer.data(dp_geom_10_off + 0 * acomps  + 11);

        #pragma omp simd aligned(cd_y, g_x_0_y_x, g_x_0_y_xy, g_x_0_y_y, g_x_0_y_yy, g_x_0_y_yz, g_x_0_y_z, g_x_0_yy_x, g_x_0_yy_y, g_x_0_yy_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yy_x[k] = -g_x_0_y_x[k] * cd_y[k] + g_x_0_y_xy[k];

            g_x_0_yy_y[k] = -g_x_0_y_y[k] * cd_y[k] + g_x_0_y_yy[k];

            g_x_0_yy_z[k] = -g_x_0_y_z[k] * cd_y[k] + g_x_0_y_yz[k];
        }

        /// Set up 12-15 components of targeted buffer : cbuffer.data(

        auto g_x_0_yz_x = cbuffer.data(dp_geom_10_off + 0 * acomps  + 12);

        auto g_x_0_yz_y = cbuffer.data(dp_geom_10_off + 0 * acomps  + 13);

        auto g_x_0_yz_z = cbuffer.data(dp_geom_10_off + 0 * acomps  + 14);

        #pragma omp simd aligned(cd_y, g_x_0_yz_x, g_x_0_yz_y, g_x_0_yz_z, g_x_0_z_x, g_x_0_z_xy, g_x_0_z_y, g_x_0_z_yy, g_x_0_z_yz, g_x_0_z_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yz_x[k] = -g_x_0_z_x[k] * cd_y[k] + g_x_0_z_xy[k];

            g_x_0_yz_y[k] = -g_x_0_z_y[k] * cd_y[k] + g_x_0_z_yy[k];

            g_x_0_yz_z[k] = -g_x_0_z_z[k] * cd_y[k] + g_x_0_z_yz[k];
        }

        /// Set up 15-18 components of targeted buffer : cbuffer.data(

        auto g_x_0_zz_x = cbuffer.data(dp_geom_10_off + 0 * acomps  + 15);

        auto g_x_0_zz_y = cbuffer.data(dp_geom_10_off + 0 * acomps  + 16);

        auto g_x_0_zz_z = cbuffer.data(dp_geom_10_off + 0 * acomps  + 17);

        #pragma omp simd aligned(cd_z, g_x_0_z_x, g_x_0_z_xz, g_x_0_z_y, g_x_0_z_yz, g_x_0_z_z, g_x_0_z_zz, g_x_0_zz_x, g_x_0_zz_y, g_x_0_zz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_zz_x[k] = -g_x_0_z_x[k] * cd_z[k] + g_x_0_z_xz[k];

            g_x_0_zz_y[k] = -g_x_0_z_y[k] * cd_z[k] + g_x_0_z_yz[k];

            g_x_0_zz_z[k] = -g_x_0_z_z[k] * cd_z[k] + g_x_0_z_zz[k];
        }
        /// Set up 0-3 components of targeted buffer : cbuffer.data(

        auto g_y_0_xx_x = cbuffer.data(dp_geom_10_off + 18 * acomps  + 0);

        auto g_y_0_xx_y = cbuffer.data(dp_geom_10_off + 18 * acomps  + 1);

        auto g_y_0_xx_z = cbuffer.data(dp_geom_10_off + 18 * acomps  + 2);

        #pragma omp simd aligned(cd_x, g_y_0_x_x, g_y_0_x_xx, g_y_0_x_xy, g_y_0_x_xz, g_y_0_x_y, g_y_0_x_z, g_y_0_xx_x, g_y_0_xx_y, g_y_0_xx_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xx_x[k] = -g_y_0_x_x[k] * cd_x[k] + g_y_0_x_xx[k];

            g_y_0_xx_y[k] = -g_y_0_x_y[k] * cd_x[k] + g_y_0_x_xy[k];

            g_y_0_xx_z[k] = -g_y_0_x_z[k] * cd_x[k] + g_y_0_x_xz[k];
        }

        /// Set up 3-6 components of targeted buffer : cbuffer.data(

        auto g_y_0_xy_x = cbuffer.data(dp_geom_10_off + 18 * acomps  + 3);

        auto g_y_0_xy_y = cbuffer.data(dp_geom_10_off + 18 * acomps  + 4);

        auto g_y_0_xy_z = cbuffer.data(dp_geom_10_off + 18 * acomps  + 5);

        #pragma omp simd aligned(cd_x, g_y_0_xy_x, g_y_0_xy_y, g_y_0_xy_z, g_y_0_y_x, g_y_0_y_xx, g_y_0_y_xy, g_y_0_y_xz, g_y_0_y_y, g_y_0_y_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xy_x[k] = -g_y_0_y_x[k] * cd_x[k] + g_y_0_y_xx[k];

            g_y_0_xy_y[k] = -g_y_0_y_y[k] * cd_x[k] + g_y_0_y_xy[k];

            g_y_0_xy_z[k] = -g_y_0_y_z[k] * cd_x[k] + g_y_0_y_xz[k];
        }

        /// Set up 6-9 components of targeted buffer : cbuffer.data(

        auto g_y_0_xz_x = cbuffer.data(dp_geom_10_off + 18 * acomps  + 6);

        auto g_y_0_xz_y = cbuffer.data(dp_geom_10_off + 18 * acomps  + 7);

        auto g_y_0_xz_z = cbuffer.data(dp_geom_10_off + 18 * acomps  + 8);

        #pragma omp simd aligned(cd_x, g_y_0_xz_x, g_y_0_xz_y, g_y_0_xz_z, g_y_0_z_x, g_y_0_z_xx, g_y_0_z_xy, g_y_0_z_xz, g_y_0_z_y, g_y_0_z_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xz_x[k] = -g_y_0_z_x[k] * cd_x[k] + g_y_0_z_xx[k];

            g_y_0_xz_y[k] = -g_y_0_z_y[k] * cd_x[k] + g_y_0_z_xy[k];

            g_y_0_xz_z[k] = -g_y_0_z_z[k] * cd_x[k] + g_y_0_z_xz[k];
        }

        /// Set up 9-12 components of targeted buffer : cbuffer.data(

        auto g_y_0_yy_x = cbuffer.data(dp_geom_10_off + 18 * acomps  + 9);

        auto g_y_0_yy_y = cbuffer.data(dp_geom_10_off + 18 * acomps  + 10);

        auto g_y_0_yy_z = cbuffer.data(dp_geom_10_off + 18 * acomps  + 11);

        #pragma omp simd aligned(cd_y, g_y_0_y_x, g_y_0_y_xy, g_y_0_y_y, g_y_0_y_yy, g_y_0_y_yz, g_y_0_y_z, g_y_0_yy_x, g_y_0_yy_y, g_y_0_yy_z, g_y_x, g_y_y, g_y_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yy_x[k] = -g_y_x[k] - g_y_0_y_x[k] * cd_y[k] + g_y_0_y_xy[k];

            g_y_0_yy_y[k] = -g_y_y[k] - g_y_0_y_y[k] * cd_y[k] + g_y_0_y_yy[k];

            g_y_0_yy_z[k] = -g_y_z[k] - g_y_0_y_z[k] * cd_y[k] + g_y_0_y_yz[k];
        }

        /// Set up 12-15 components of targeted buffer : cbuffer.data(

        auto g_y_0_yz_x = cbuffer.data(dp_geom_10_off + 18 * acomps  + 12);

        auto g_y_0_yz_y = cbuffer.data(dp_geom_10_off + 18 * acomps  + 13);

        auto g_y_0_yz_z = cbuffer.data(dp_geom_10_off + 18 * acomps  + 14);

        #pragma omp simd aligned(cd_z, g_y_0_y_x, g_y_0_y_xz, g_y_0_y_y, g_y_0_y_yz, g_y_0_y_z, g_y_0_y_zz, g_y_0_yz_x, g_y_0_yz_y, g_y_0_yz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yz_x[k] = -g_y_0_y_x[k] * cd_z[k] + g_y_0_y_xz[k];

            g_y_0_yz_y[k] = -g_y_0_y_y[k] * cd_z[k] + g_y_0_y_yz[k];

            g_y_0_yz_z[k] = -g_y_0_y_z[k] * cd_z[k] + g_y_0_y_zz[k];
        }

        /// Set up 15-18 components of targeted buffer : cbuffer.data(

        auto g_y_0_zz_x = cbuffer.data(dp_geom_10_off + 18 * acomps  + 15);

        auto g_y_0_zz_y = cbuffer.data(dp_geom_10_off + 18 * acomps  + 16);

        auto g_y_0_zz_z = cbuffer.data(dp_geom_10_off + 18 * acomps  + 17);

        #pragma omp simd aligned(cd_z, g_y_0_z_x, g_y_0_z_xz, g_y_0_z_y, g_y_0_z_yz, g_y_0_z_z, g_y_0_z_zz, g_y_0_zz_x, g_y_0_zz_y, g_y_0_zz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_zz_x[k] = -g_y_0_z_x[k] * cd_z[k] + g_y_0_z_xz[k];

            g_y_0_zz_y[k] = -g_y_0_z_y[k] * cd_z[k] + g_y_0_z_yz[k];

            g_y_0_zz_z[k] = -g_y_0_z_z[k] * cd_z[k] + g_y_0_z_zz[k];
        }
        /// Set up 0-3 components of targeted buffer : cbuffer.data(

        auto g_z_0_xx_x = cbuffer.data(dp_geom_10_off + 36 * acomps  + 0);

        auto g_z_0_xx_y = cbuffer.data(dp_geom_10_off + 36 * acomps  + 1);

        auto g_z_0_xx_z = cbuffer.data(dp_geom_10_off + 36 * acomps  + 2);

        #pragma omp simd aligned(cd_x, g_z_0_x_x, g_z_0_x_xx, g_z_0_x_xy, g_z_0_x_xz, g_z_0_x_y, g_z_0_x_z, g_z_0_xx_x, g_z_0_xx_y, g_z_0_xx_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xx_x[k] = -g_z_0_x_x[k] * cd_x[k] + g_z_0_x_xx[k];

            g_z_0_xx_y[k] = -g_z_0_x_y[k] * cd_x[k] + g_z_0_x_xy[k];

            g_z_0_xx_z[k] = -g_z_0_x_z[k] * cd_x[k] + g_z_0_x_xz[k];
        }

        /// Set up 3-6 components of targeted buffer : cbuffer.data(

        auto g_z_0_xy_x = cbuffer.data(dp_geom_10_off + 36 * acomps  + 3);

        auto g_z_0_xy_y = cbuffer.data(dp_geom_10_off + 36 * acomps  + 4);

        auto g_z_0_xy_z = cbuffer.data(dp_geom_10_off + 36 * acomps  + 5);

        #pragma omp simd aligned(cd_x, g_z_0_xy_x, g_z_0_xy_y, g_z_0_xy_z, g_z_0_y_x, g_z_0_y_xx, g_z_0_y_xy, g_z_0_y_xz, g_z_0_y_y, g_z_0_y_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xy_x[k] = -g_z_0_y_x[k] * cd_x[k] + g_z_0_y_xx[k];

            g_z_0_xy_y[k] = -g_z_0_y_y[k] * cd_x[k] + g_z_0_y_xy[k];

            g_z_0_xy_z[k] = -g_z_0_y_z[k] * cd_x[k] + g_z_0_y_xz[k];
        }

        /// Set up 6-9 components of targeted buffer : cbuffer.data(

        auto g_z_0_xz_x = cbuffer.data(dp_geom_10_off + 36 * acomps  + 6);

        auto g_z_0_xz_y = cbuffer.data(dp_geom_10_off + 36 * acomps  + 7);

        auto g_z_0_xz_z = cbuffer.data(dp_geom_10_off + 36 * acomps  + 8);

        #pragma omp simd aligned(cd_x, g_z_0_xz_x, g_z_0_xz_y, g_z_0_xz_z, g_z_0_z_x, g_z_0_z_xx, g_z_0_z_xy, g_z_0_z_xz, g_z_0_z_y, g_z_0_z_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xz_x[k] = -g_z_0_z_x[k] * cd_x[k] + g_z_0_z_xx[k];

            g_z_0_xz_y[k] = -g_z_0_z_y[k] * cd_x[k] + g_z_0_z_xy[k];

            g_z_0_xz_z[k] = -g_z_0_z_z[k] * cd_x[k] + g_z_0_z_xz[k];
        }

        /// Set up 9-12 components of targeted buffer : cbuffer.data(

        auto g_z_0_yy_x = cbuffer.data(dp_geom_10_off + 36 * acomps  + 9);

        auto g_z_0_yy_y = cbuffer.data(dp_geom_10_off + 36 * acomps  + 10);

        auto g_z_0_yy_z = cbuffer.data(dp_geom_10_off + 36 * acomps  + 11);

        #pragma omp simd aligned(cd_y, g_z_0_y_x, g_z_0_y_xy, g_z_0_y_y, g_z_0_y_yy, g_z_0_y_yz, g_z_0_y_z, g_z_0_yy_x, g_z_0_yy_y, g_z_0_yy_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yy_x[k] = -g_z_0_y_x[k] * cd_y[k] + g_z_0_y_xy[k];

            g_z_0_yy_y[k] = -g_z_0_y_y[k] * cd_y[k] + g_z_0_y_yy[k];

            g_z_0_yy_z[k] = -g_z_0_y_z[k] * cd_y[k] + g_z_0_y_yz[k];
        }

        /// Set up 12-15 components of targeted buffer : cbuffer.data(

        auto g_z_0_yz_x = cbuffer.data(dp_geom_10_off + 36 * acomps  + 12);

        auto g_z_0_yz_y = cbuffer.data(dp_geom_10_off + 36 * acomps  + 13);

        auto g_z_0_yz_z = cbuffer.data(dp_geom_10_off + 36 * acomps  + 14);

        #pragma omp simd aligned(cd_y, g_z_0_yz_x, g_z_0_yz_y, g_z_0_yz_z, g_z_0_z_x, g_z_0_z_xy, g_z_0_z_y, g_z_0_z_yy, g_z_0_z_yz, g_z_0_z_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yz_x[k] = -g_z_0_z_x[k] * cd_y[k] + g_z_0_z_xy[k];

            g_z_0_yz_y[k] = -g_z_0_z_y[k] * cd_y[k] + g_z_0_z_yy[k];

            g_z_0_yz_z[k] = -g_z_0_z_z[k] * cd_y[k] + g_z_0_z_yz[k];
        }

        /// Set up 15-18 components of targeted buffer : cbuffer.data(

        auto g_z_0_zz_x = cbuffer.data(dp_geom_10_off + 36 * acomps  + 15);

        auto g_z_0_zz_y = cbuffer.data(dp_geom_10_off + 36 * acomps  + 16);

        auto g_z_0_zz_z = cbuffer.data(dp_geom_10_off + 36 * acomps  + 17);

        #pragma omp simd aligned(cd_z, g_z_0_z_x, g_z_0_z_xz, g_z_0_z_y, g_z_0_z_yz, g_z_0_z_z, g_z_0_z_zz, g_z_0_zz_x, g_z_0_zz_y, g_z_0_zz_z, g_z_x, g_z_y, g_z_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_zz_x[k] = -g_z_x[k] - g_z_0_z_x[k] * cd_z[k] + g_z_0_z_xz[k];

            g_z_0_zz_y[k] = -g_z_y[k] - g_z_0_z_y[k] * cd_z[k] + g_z_0_z_yz[k];

            g_z_0_zz_z[k] = -g_z_z[k] - g_z_0_z_z[k] * cd_z[k] + g_z_0_z_zz[k];
        }
    }
}

} // t3ceri namespace

