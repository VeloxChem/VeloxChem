#include "ElectronRepulsionGeom1100ContrRecPPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom11_hrr_electron_repulsion_ppxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_ppxx,
                                            const size_t idx_geom_01_spxx,
                                            const size_t idx_geom_10_spxx,
                                            const size_t idx_geom_11_spxx,
                                            const size_t idx_geom_11_sdxx,
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
            /// Set up components of auxilary buffer : SPSS

            const auto sp_geom_01_off = idx_geom_01_spxx + i * dcomps + j;

            auto g_0_x_0_x = cbuffer.data(sp_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_0_y = cbuffer.data(sp_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_0_z = cbuffer.data(sp_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_y_0_x = cbuffer.data(sp_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_y_0_y = cbuffer.data(sp_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_y_0_z = cbuffer.data(sp_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_z_0_x = cbuffer.data(sp_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_z_0_y = cbuffer.data(sp_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_z_0_z = cbuffer.data(sp_geom_01_off + 8 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SPSS

            const auto sp_geom_10_off = idx_geom_10_spxx + i * dcomps + j;

            auto g_x_0_0_x = cbuffer.data(sp_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_0_y = cbuffer.data(sp_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_0_z = cbuffer.data(sp_geom_10_off + 2 * ccomps * dcomps);

            auto g_y_0_0_x = cbuffer.data(sp_geom_10_off + 3 * ccomps * dcomps);

            auto g_y_0_0_y = cbuffer.data(sp_geom_10_off + 4 * ccomps * dcomps);

            auto g_y_0_0_z = cbuffer.data(sp_geom_10_off + 5 * ccomps * dcomps);

            auto g_z_0_0_x = cbuffer.data(sp_geom_10_off + 6 * ccomps * dcomps);

            auto g_z_0_0_y = cbuffer.data(sp_geom_10_off + 7 * ccomps * dcomps);

            auto g_z_0_0_z = cbuffer.data(sp_geom_10_off + 8 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SPSS

            const auto sp_geom_11_off = idx_geom_11_spxx + i * dcomps + j;

            auto g_x_x_0_x = cbuffer.data(sp_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_0_y = cbuffer.data(sp_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_0_z = cbuffer.data(sp_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_y_0_x = cbuffer.data(sp_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_y_0_y = cbuffer.data(sp_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_y_0_z = cbuffer.data(sp_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_z_0_x = cbuffer.data(sp_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_z_0_y = cbuffer.data(sp_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_z_0_z = cbuffer.data(sp_geom_11_off + 8 * ccomps * dcomps);

            auto g_y_x_0_x = cbuffer.data(sp_geom_11_off + 9 * ccomps * dcomps);

            auto g_y_x_0_y = cbuffer.data(sp_geom_11_off + 10 * ccomps * dcomps);

            auto g_y_x_0_z = cbuffer.data(sp_geom_11_off + 11 * ccomps * dcomps);

            auto g_y_y_0_x = cbuffer.data(sp_geom_11_off + 12 * ccomps * dcomps);

            auto g_y_y_0_y = cbuffer.data(sp_geom_11_off + 13 * ccomps * dcomps);

            auto g_y_y_0_z = cbuffer.data(sp_geom_11_off + 14 * ccomps * dcomps);

            auto g_y_z_0_x = cbuffer.data(sp_geom_11_off + 15 * ccomps * dcomps);

            auto g_y_z_0_y = cbuffer.data(sp_geom_11_off + 16 * ccomps * dcomps);

            auto g_y_z_0_z = cbuffer.data(sp_geom_11_off + 17 * ccomps * dcomps);

            auto g_z_x_0_x = cbuffer.data(sp_geom_11_off + 18 * ccomps * dcomps);

            auto g_z_x_0_y = cbuffer.data(sp_geom_11_off + 19 * ccomps * dcomps);

            auto g_z_x_0_z = cbuffer.data(sp_geom_11_off + 20 * ccomps * dcomps);

            auto g_z_y_0_x = cbuffer.data(sp_geom_11_off + 21 * ccomps * dcomps);

            auto g_z_y_0_y = cbuffer.data(sp_geom_11_off + 22 * ccomps * dcomps);

            auto g_z_y_0_z = cbuffer.data(sp_geom_11_off + 23 * ccomps * dcomps);

            auto g_z_z_0_x = cbuffer.data(sp_geom_11_off + 24 * ccomps * dcomps);

            auto g_z_z_0_y = cbuffer.data(sp_geom_11_off + 25 * ccomps * dcomps);

            auto g_z_z_0_z = cbuffer.data(sp_geom_11_off + 26 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SDSS

            const auto sd_geom_11_off = idx_geom_11_sdxx + i * dcomps + j;

            auto g_x_x_0_xx = cbuffer.data(sd_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_0_xy = cbuffer.data(sd_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_0_xz = cbuffer.data(sd_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_0_yy = cbuffer.data(sd_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_0_yz = cbuffer.data(sd_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_0_zz = cbuffer.data(sd_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_y_0_xx = cbuffer.data(sd_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_y_0_xy = cbuffer.data(sd_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_y_0_xz = cbuffer.data(sd_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_y_0_yy = cbuffer.data(sd_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_y_0_yz = cbuffer.data(sd_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_y_0_zz = cbuffer.data(sd_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_z_0_xx = cbuffer.data(sd_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_z_0_xy = cbuffer.data(sd_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_z_0_xz = cbuffer.data(sd_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_z_0_yy = cbuffer.data(sd_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_z_0_yz = cbuffer.data(sd_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_z_0_zz = cbuffer.data(sd_geom_11_off + 17 * ccomps * dcomps);

            auto g_y_x_0_xx = cbuffer.data(sd_geom_11_off + 18 * ccomps * dcomps);

            auto g_y_x_0_xy = cbuffer.data(sd_geom_11_off + 19 * ccomps * dcomps);

            auto g_y_x_0_xz = cbuffer.data(sd_geom_11_off + 20 * ccomps * dcomps);

            auto g_y_x_0_yy = cbuffer.data(sd_geom_11_off + 21 * ccomps * dcomps);

            auto g_y_x_0_yz = cbuffer.data(sd_geom_11_off + 22 * ccomps * dcomps);

            auto g_y_x_0_zz = cbuffer.data(sd_geom_11_off + 23 * ccomps * dcomps);

            auto g_y_y_0_xx = cbuffer.data(sd_geom_11_off + 24 * ccomps * dcomps);

            auto g_y_y_0_xy = cbuffer.data(sd_geom_11_off + 25 * ccomps * dcomps);

            auto g_y_y_0_xz = cbuffer.data(sd_geom_11_off + 26 * ccomps * dcomps);

            auto g_y_y_0_yy = cbuffer.data(sd_geom_11_off + 27 * ccomps * dcomps);

            auto g_y_y_0_yz = cbuffer.data(sd_geom_11_off + 28 * ccomps * dcomps);

            auto g_y_y_0_zz = cbuffer.data(sd_geom_11_off + 29 * ccomps * dcomps);

            auto g_y_z_0_xx = cbuffer.data(sd_geom_11_off + 30 * ccomps * dcomps);

            auto g_y_z_0_xy = cbuffer.data(sd_geom_11_off + 31 * ccomps * dcomps);

            auto g_y_z_0_xz = cbuffer.data(sd_geom_11_off + 32 * ccomps * dcomps);

            auto g_y_z_0_yy = cbuffer.data(sd_geom_11_off + 33 * ccomps * dcomps);

            auto g_y_z_0_yz = cbuffer.data(sd_geom_11_off + 34 * ccomps * dcomps);

            auto g_y_z_0_zz = cbuffer.data(sd_geom_11_off + 35 * ccomps * dcomps);

            auto g_z_x_0_xx = cbuffer.data(sd_geom_11_off + 36 * ccomps * dcomps);

            auto g_z_x_0_xy = cbuffer.data(sd_geom_11_off + 37 * ccomps * dcomps);

            auto g_z_x_0_xz = cbuffer.data(sd_geom_11_off + 38 * ccomps * dcomps);

            auto g_z_x_0_yy = cbuffer.data(sd_geom_11_off + 39 * ccomps * dcomps);

            auto g_z_x_0_yz = cbuffer.data(sd_geom_11_off + 40 * ccomps * dcomps);

            auto g_z_x_0_zz = cbuffer.data(sd_geom_11_off + 41 * ccomps * dcomps);

            auto g_z_y_0_xx = cbuffer.data(sd_geom_11_off + 42 * ccomps * dcomps);

            auto g_z_y_0_xy = cbuffer.data(sd_geom_11_off + 43 * ccomps * dcomps);

            auto g_z_y_0_xz = cbuffer.data(sd_geom_11_off + 44 * ccomps * dcomps);

            auto g_z_y_0_yy = cbuffer.data(sd_geom_11_off + 45 * ccomps * dcomps);

            auto g_z_y_0_yz = cbuffer.data(sd_geom_11_off + 46 * ccomps * dcomps);

            auto g_z_y_0_zz = cbuffer.data(sd_geom_11_off + 47 * ccomps * dcomps);

            auto g_z_z_0_xx = cbuffer.data(sd_geom_11_off + 48 * ccomps * dcomps);

            auto g_z_z_0_xy = cbuffer.data(sd_geom_11_off + 49 * ccomps * dcomps);

            auto g_z_z_0_xz = cbuffer.data(sd_geom_11_off + 50 * ccomps * dcomps);

            auto g_z_z_0_yy = cbuffer.data(sd_geom_11_off + 51 * ccomps * dcomps);

            auto g_z_z_0_yz = cbuffer.data(sd_geom_11_off + 52 * ccomps * dcomps);

            auto g_z_z_0_zz = cbuffer.data(sd_geom_11_off + 53 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_ppxx

            const auto pp_geom_11_off = idx_geom_11_ppxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_x_x_x_x = cbuffer.data(pp_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_x_y = cbuffer.data(pp_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_x_z = cbuffer.data(pp_geom_11_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_x_0_0_x, g_x_0_0_y, g_x_0_0_z, g_x_x_0_x, g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_y, g_x_x_0_z, g_x_x_x_x, g_x_x_x_y, g_x_x_x_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_x_x[k] = -g_0_x_0_x[k] + g_x_0_0_x[k] - g_x_x_0_x[k] * ab_x + g_x_x_0_xx[k];

                g_x_x_x_y[k] = -g_0_x_0_y[k] + g_x_0_0_y[k] - g_x_x_0_y[k] * ab_x + g_x_x_0_xy[k];

                g_x_x_x_z[k] = -g_0_x_0_z[k] + g_x_0_0_z[k] - g_x_x_0_z[k] * ab_x + g_x_x_0_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_x_x_y_x = cbuffer.data(pp_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_y_y = cbuffer.data(pp_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_y_z = cbuffer.data(pp_geom_11_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_0_x, g_x_x_0_xy, g_x_x_0_y, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_z, g_x_x_y_x, g_x_x_y_y, g_x_x_y_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_y_x[k] = -g_x_x_0_x[k] * ab_y + g_x_x_0_xy[k];

                g_x_x_y_y[k] = -g_x_x_0_y[k] * ab_y + g_x_x_0_yy[k];

                g_x_x_y_z[k] = -g_x_x_0_z[k] * ab_y + g_x_x_0_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_x_x_z_x = cbuffer.data(pp_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_z_y = cbuffer.data(pp_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_z_z = cbuffer.data(pp_geom_11_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_0_x, g_x_x_0_xz, g_x_x_0_y, g_x_x_0_yz, g_x_x_0_z, g_x_x_0_zz, g_x_x_z_x, g_x_x_z_y, g_x_x_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_z_x[k] = -g_x_x_0_x[k] * ab_z + g_x_x_0_xz[k];

                g_x_x_z_y[k] = -g_x_x_0_y[k] * ab_z + g_x_x_0_yz[k];

                g_x_x_z_z[k] = -g_x_x_0_z[k] * ab_z + g_x_x_0_zz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_x_y_x_x = cbuffer.data(pp_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_y_x_y = cbuffer.data(pp_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_y_x_z = cbuffer.data(pp_geom_11_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_x_y_0_x, g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_y, g_x_y_0_z, g_x_y_x_x, g_x_y_x_y, g_x_y_x_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_x_x[k] = -g_0_y_0_x[k] - g_x_y_0_x[k] * ab_x + g_x_y_0_xx[k];

                g_x_y_x_y[k] = -g_0_y_0_y[k] - g_x_y_0_y[k] * ab_x + g_x_y_0_xy[k];

                g_x_y_x_z[k] = -g_0_y_0_z[k] - g_x_y_0_z[k] * ab_x + g_x_y_0_xz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_x_y_y_x = cbuffer.data(pp_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_y_y_y = cbuffer.data(pp_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_y_y_z = cbuffer.data(pp_geom_11_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_x, g_x_0_0_y, g_x_0_0_z, g_x_y_0_x, g_x_y_0_xy, g_x_y_0_y, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_z, g_x_y_y_x, g_x_y_y_y, g_x_y_y_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_y_x[k] = g_x_0_0_x[k] - g_x_y_0_x[k] * ab_y + g_x_y_0_xy[k];

                g_x_y_y_y[k] = g_x_0_0_y[k] - g_x_y_0_y[k] * ab_y + g_x_y_0_yy[k];

                g_x_y_y_z[k] = g_x_0_0_z[k] - g_x_y_0_z[k] * ab_y + g_x_y_0_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_x_y_z_x = cbuffer.data(pp_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_y_z_y = cbuffer.data(pp_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_y_z_z = cbuffer.data(pp_geom_11_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_0_x, g_x_y_0_xz, g_x_y_0_y, g_x_y_0_yz, g_x_y_0_z, g_x_y_0_zz, g_x_y_z_x, g_x_y_z_y, g_x_y_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_z_x[k] = -g_x_y_0_x[k] * ab_z + g_x_y_0_xz[k];

                g_x_y_z_y[k] = -g_x_y_0_y[k] * ab_z + g_x_y_0_yz[k];

                g_x_y_z_z[k] = -g_x_y_0_z[k] * ab_z + g_x_y_0_zz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_x_z_x_x = cbuffer.data(pp_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_z_x_y = cbuffer.data(pp_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_z_x_z = cbuffer.data(pp_geom_11_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_x_z_0_x, g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_y, g_x_z_0_z, g_x_z_x_x, g_x_z_x_y, g_x_z_x_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_x_x[k] = -g_0_z_0_x[k] - g_x_z_0_x[k] * ab_x + g_x_z_0_xx[k];

                g_x_z_x_y[k] = -g_0_z_0_y[k] - g_x_z_0_y[k] * ab_x + g_x_z_0_xy[k];

                g_x_z_x_z[k] = -g_0_z_0_z[k] - g_x_z_0_z[k] * ab_x + g_x_z_0_xz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_x_z_y_x = cbuffer.data(pp_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_z_y_y = cbuffer.data(pp_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_z_y_z = cbuffer.data(pp_geom_11_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_0_x, g_x_z_0_xy, g_x_z_0_y, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_z, g_x_z_y_x, g_x_z_y_y, g_x_z_y_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_y_x[k] = -g_x_z_0_x[k] * ab_y + g_x_z_0_xy[k];

                g_x_z_y_y[k] = -g_x_z_0_y[k] * ab_y + g_x_z_0_yy[k];

                g_x_z_y_z[k] = -g_x_z_0_z[k] * ab_y + g_x_z_0_yz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_x_z_z_x = cbuffer.data(pp_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_z_z_y = cbuffer.data(pp_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_z_z_z = cbuffer.data(pp_geom_11_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_x, g_x_0_0_y, g_x_0_0_z, g_x_z_0_x, g_x_z_0_xz, g_x_z_0_y, g_x_z_0_yz, g_x_z_0_z, g_x_z_0_zz, g_x_z_z_x, g_x_z_z_y, g_x_z_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_z_x[k] = g_x_0_0_x[k] - g_x_z_0_x[k] * ab_z + g_x_z_0_xz[k];

                g_x_z_z_y[k] = g_x_0_0_y[k] - g_x_z_0_y[k] * ab_z + g_x_z_0_yz[k];

                g_x_z_z_z[k] = g_x_0_0_z[k] - g_x_z_0_z[k] * ab_z + g_x_z_0_zz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_y_x_x_x = cbuffer.data(pp_geom_11_off + 27 * ccomps * dcomps);

            auto g_y_x_x_y = cbuffer.data(pp_geom_11_off + 28 * ccomps * dcomps);

            auto g_y_x_x_z = cbuffer.data(pp_geom_11_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_0_x, g_y_0_0_y, g_y_0_0_z, g_y_x_0_x, g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_y, g_y_x_0_z, g_y_x_x_x, g_y_x_x_y, g_y_x_x_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_x_x[k] = g_y_0_0_x[k] - g_y_x_0_x[k] * ab_x + g_y_x_0_xx[k];

                g_y_x_x_y[k] = g_y_0_0_y[k] - g_y_x_0_y[k] * ab_x + g_y_x_0_xy[k];

                g_y_x_x_z[k] = g_y_0_0_z[k] - g_y_x_0_z[k] * ab_x + g_y_x_0_xz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_y_x_y_x = cbuffer.data(pp_geom_11_off + 30 * ccomps * dcomps);

            auto g_y_x_y_y = cbuffer.data(pp_geom_11_off + 31 * ccomps * dcomps);

            auto g_y_x_y_z = cbuffer.data(pp_geom_11_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_y_x_0_x, g_y_x_0_xy, g_y_x_0_y, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_z, g_y_x_y_x, g_y_x_y_y, g_y_x_y_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_y_x[k] = -g_0_x_0_x[k] - g_y_x_0_x[k] * ab_y + g_y_x_0_xy[k];

                g_y_x_y_y[k] = -g_0_x_0_y[k] - g_y_x_0_y[k] * ab_y + g_y_x_0_yy[k];

                g_y_x_y_z[k] = -g_0_x_0_z[k] - g_y_x_0_z[k] * ab_y + g_y_x_0_yz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_y_x_z_x = cbuffer.data(pp_geom_11_off + 33 * ccomps * dcomps);

            auto g_y_x_z_y = cbuffer.data(pp_geom_11_off + 34 * ccomps * dcomps);

            auto g_y_x_z_z = cbuffer.data(pp_geom_11_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_0_x, g_y_x_0_xz, g_y_x_0_y, g_y_x_0_yz, g_y_x_0_z, g_y_x_0_zz, g_y_x_z_x, g_y_x_z_y, g_y_x_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_z_x[k] = -g_y_x_0_x[k] * ab_z + g_y_x_0_xz[k];

                g_y_x_z_y[k] = -g_y_x_0_y[k] * ab_z + g_y_x_0_yz[k];

                g_y_x_z_z[k] = -g_y_x_0_z[k] * ab_z + g_y_x_0_zz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_y_y_x_x = cbuffer.data(pp_geom_11_off + 36 * ccomps * dcomps);

            auto g_y_y_x_y = cbuffer.data(pp_geom_11_off + 37 * ccomps * dcomps);

            auto g_y_y_x_z = cbuffer.data(pp_geom_11_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_0_x, g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_y, g_y_y_0_z, g_y_y_x_x, g_y_y_x_y, g_y_y_x_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_x_x[k] = -g_y_y_0_x[k] * ab_x + g_y_y_0_xx[k];

                g_y_y_x_y[k] = -g_y_y_0_y[k] * ab_x + g_y_y_0_xy[k];

                g_y_y_x_z[k] = -g_y_y_0_z[k] * ab_x + g_y_y_0_xz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_y_y_y_x = cbuffer.data(pp_geom_11_off + 39 * ccomps * dcomps);

            auto g_y_y_y_y = cbuffer.data(pp_geom_11_off + 40 * ccomps * dcomps);

            auto g_y_y_y_z = cbuffer.data(pp_geom_11_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_y_0_0_x, g_y_0_0_y, g_y_0_0_z, g_y_y_0_x, g_y_y_0_xy, g_y_y_0_y, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_z, g_y_y_y_x, g_y_y_y_y, g_y_y_y_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_y_x[k] = -g_0_y_0_x[k] + g_y_0_0_x[k] - g_y_y_0_x[k] * ab_y + g_y_y_0_xy[k];

                g_y_y_y_y[k] = -g_0_y_0_y[k] + g_y_0_0_y[k] - g_y_y_0_y[k] * ab_y + g_y_y_0_yy[k];

                g_y_y_y_z[k] = -g_0_y_0_z[k] + g_y_0_0_z[k] - g_y_y_0_z[k] * ab_y + g_y_y_0_yz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_y_y_z_x = cbuffer.data(pp_geom_11_off + 42 * ccomps * dcomps);

            auto g_y_y_z_y = cbuffer.data(pp_geom_11_off + 43 * ccomps * dcomps);

            auto g_y_y_z_z = cbuffer.data(pp_geom_11_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_0_x, g_y_y_0_xz, g_y_y_0_y, g_y_y_0_yz, g_y_y_0_z, g_y_y_0_zz, g_y_y_z_x, g_y_y_z_y, g_y_y_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_z_x[k] = -g_y_y_0_x[k] * ab_z + g_y_y_0_xz[k];

                g_y_y_z_y[k] = -g_y_y_0_y[k] * ab_z + g_y_y_0_yz[k];

                g_y_y_z_z[k] = -g_y_y_0_z[k] * ab_z + g_y_y_0_zz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_y_z_x_x = cbuffer.data(pp_geom_11_off + 45 * ccomps * dcomps);

            auto g_y_z_x_y = cbuffer.data(pp_geom_11_off + 46 * ccomps * dcomps);

            auto g_y_z_x_z = cbuffer.data(pp_geom_11_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_0_x, g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_y, g_y_z_0_z, g_y_z_x_x, g_y_z_x_y, g_y_z_x_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_x_x[k] = -g_y_z_0_x[k] * ab_x + g_y_z_0_xx[k];

                g_y_z_x_y[k] = -g_y_z_0_y[k] * ab_x + g_y_z_0_xy[k];

                g_y_z_x_z[k] = -g_y_z_0_z[k] * ab_x + g_y_z_0_xz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_y_z_y_x = cbuffer.data(pp_geom_11_off + 48 * ccomps * dcomps);

            auto g_y_z_y_y = cbuffer.data(pp_geom_11_off + 49 * ccomps * dcomps);

            auto g_y_z_y_z = cbuffer.data(pp_geom_11_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_y_z_0_x, g_y_z_0_xy, g_y_z_0_y, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_z, g_y_z_y_x, g_y_z_y_y, g_y_z_y_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_y_x[k] = -g_0_z_0_x[k] - g_y_z_0_x[k] * ab_y + g_y_z_0_xy[k];

                g_y_z_y_y[k] = -g_0_z_0_y[k] - g_y_z_0_y[k] * ab_y + g_y_z_0_yy[k];

                g_y_z_y_z[k] = -g_0_z_0_z[k] - g_y_z_0_z[k] * ab_y + g_y_z_0_yz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_y_z_z_x = cbuffer.data(pp_geom_11_off + 51 * ccomps * dcomps);

            auto g_y_z_z_y = cbuffer.data(pp_geom_11_off + 52 * ccomps * dcomps);

            auto g_y_z_z_z = cbuffer.data(pp_geom_11_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_0_x, g_y_0_0_y, g_y_0_0_z, g_y_z_0_x, g_y_z_0_xz, g_y_z_0_y, g_y_z_0_yz, g_y_z_0_z, g_y_z_0_zz, g_y_z_z_x, g_y_z_z_y, g_y_z_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_z_x[k] = g_y_0_0_x[k] - g_y_z_0_x[k] * ab_z + g_y_z_0_xz[k];

                g_y_z_z_y[k] = g_y_0_0_y[k] - g_y_z_0_y[k] * ab_z + g_y_z_0_yz[k];

                g_y_z_z_z[k] = g_y_0_0_z[k] - g_y_z_0_z[k] * ab_z + g_y_z_0_zz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_z_x_x_x = cbuffer.data(pp_geom_11_off + 54 * ccomps * dcomps);

            auto g_z_x_x_y = cbuffer.data(pp_geom_11_off + 55 * ccomps * dcomps);

            auto g_z_x_x_z = cbuffer.data(pp_geom_11_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_0_x, g_z_0_0_y, g_z_0_0_z, g_z_x_0_x, g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_y, g_z_x_0_z, g_z_x_x_x, g_z_x_x_y, g_z_x_x_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_x_x[k] = g_z_0_0_x[k] - g_z_x_0_x[k] * ab_x + g_z_x_0_xx[k];

                g_z_x_x_y[k] = g_z_0_0_y[k] - g_z_x_0_y[k] * ab_x + g_z_x_0_xy[k];

                g_z_x_x_z[k] = g_z_0_0_z[k] - g_z_x_0_z[k] * ab_x + g_z_x_0_xz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_z_x_y_x = cbuffer.data(pp_geom_11_off + 57 * ccomps * dcomps);

            auto g_z_x_y_y = cbuffer.data(pp_geom_11_off + 58 * ccomps * dcomps);

            auto g_z_x_y_z = cbuffer.data(pp_geom_11_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_0_x, g_z_x_0_xy, g_z_x_0_y, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_z, g_z_x_y_x, g_z_x_y_y, g_z_x_y_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_y_x[k] = -g_z_x_0_x[k] * ab_y + g_z_x_0_xy[k];

                g_z_x_y_y[k] = -g_z_x_0_y[k] * ab_y + g_z_x_0_yy[k];

                g_z_x_y_z[k] = -g_z_x_0_z[k] * ab_y + g_z_x_0_yz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_z_x_z_x = cbuffer.data(pp_geom_11_off + 60 * ccomps * dcomps);

            auto g_z_x_z_y = cbuffer.data(pp_geom_11_off + 61 * ccomps * dcomps);

            auto g_z_x_z_z = cbuffer.data(pp_geom_11_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_z_x_0_x, g_z_x_0_xz, g_z_x_0_y, g_z_x_0_yz, g_z_x_0_z, g_z_x_0_zz, g_z_x_z_x, g_z_x_z_y, g_z_x_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_z_x[k] = -g_0_x_0_x[k] - g_z_x_0_x[k] * ab_z + g_z_x_0_xz[k];

                g_z_x_z_y[k] = -g_0_x_0_y[k] - g_z_x_0_y[k] * ab_z + g_z_x_0_yz[k];

                g_z_x_z_z[k] = -g_0_x_0_z[k] - g_z_x_0_z[k] * ab_z + g_z_x_0_zz[k];
            }

            /// Set up 63-66 components of targeted buffer : cbuffer.data(

            auto g_z_y_x_x = cbuffer.data(pp_geom_11_off + 63 * ccomps * dcomps);

            auto g_z_y_x_y = cbuffer.data(pp_geom_11_off + 64 * ccomps * dcomps);

            auto g_z_y_x_z = cbuffer.data(pp_geom_11_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_0_x, g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_y, g_z_y_0_z, g_z_y_x_x, g_z_y_x_y, g_z_y_x_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_x_x[k] = -g_z_y_0_x[k] * ab_x + g_z_y_0_xx[k];

                g_z_y_x_y[k] = -g_z_y_0_y[k] * ab_x + g_z_y_0_xy[k];

                g_z_y_x_z[k] = -g_z_y_0_z[k] * ab_x + g_z_y_0_xz[k];
            }

            /// Set up 66-69 components of targeted buffer : cbuffer.data(

            auto g_z_y_y_x = cbuffer.data(pp_geom_11_off + 66 * ccomps * dcomps);

            auto g_z_y_y_y = cbuffer.data(pp_geom_11_off + 67 * ccomps * dcomps);

            auto g_z_y_y_z = cbuffer.data(pp_geom_11_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_0_x, g_z_0_0_y, g_z_0_0_z, g_z_y_0_x, g_z_y_0_xy, g_z_y_0_y, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_z, g_z_y_y_x, g_z_y_y_y, g_z_y_y_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_y_x[k] = g_z_0_0_x[k] - g_z_y_0_x[k] * ab_y + g_z_y_0_xy[k];

                g_z_y_y_y[k] = g_z_0_0_y[k] - g_z_y_0_y[k] * ab_y + g_z_y_0_yy[k];

                g_z_y_y_z[k] = g_z_0_0_z[k] - g_z_y_0_z[k] * ab_y + g_z_y_0_yz[k];
            }

            /// Set up 69-72 components of targeted buffer : cbuffer.data(

            auto g_z_y_z_x = cbuffer.data(pp_geom_11_off + 69 * ccomps * dcomps);

            auto g_z_y_z_y = cbuffer.data(pp_geom_11_off + 70 * ccomps * dcomps);

            auto g_z_y_z_z = cbuffer.data(pp_geom_11_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_z_y_0_x, g_z_y_0_xz, g_z_y_0_y, g_z_y_0_yz, g_z_y_0_z, g_z_y_0_zz, g_z_y_z_x, g_z_y_z_y, g_z_y_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_z_x[k] = -g_0_y_0_x[k] - g_z_y_0_x[k] * ab_z + g_z_y_0_xz[k];

                g_z_y_z_y[k] = -g_0_y_0_y[k] - g_z_y_0_y[k] * ab_z + g_z_y_0_yz[k];

                g_z_y_z_z[k] = -g_0_y_0_z[k] - g_z_y_0_z[k] * ab_z + g_z_y_0_zz[k];
            }

            /// Set up 72-75 components of targeted buffer : cbuffer.data(

            auto g_z_z_x_x = cbuffer.data(pp_geom_11_off + 72 * ccomps * dcomps);

            auto g_z_z_x_y = cbuffer.data(pp_geom_11_off + 73 * ccomps * dcomps);

            auto g_z_z_x_z = cbuffer.data(pp_geom_11_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_0_x, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_y, g_z_z_0_z, g_z_z_x_x, g_z_z_x_y, g_z_z_x_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_x_x[k] = -g_z_z_0_x[k] * ab_x + g_z_z_0_xx[k];

                g_z_z_x_y[k] = -g_z_z_0_y[k] * ab_x + g_z_z_0_xy[k];

                g_z_z_x_z[k] = -g_z_z_0_z[k] * ab_x + g_z_z_0_xz[k];
            }

            /// Set up 75-78 components of targeted buffer : cbuffer.data(

            auto g_z_z_y_x = cbuffer.data(pp_geom_11_off + 75 * ccomps * dcomps);

            auto g_z_z_y_y = cbuffer.data(pp_geom_11_off + 76 * ccomps * dcomps);

            auto g_z_z_y_z = cbuffer.data(pp_geom_11_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_0_x, g_z_z_0_xy, g_z_z_0_y, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_z, g_z_z_y_x, g_z_z_y_y, g_z_z_y_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_y_x[k] = -g_z_z_0_x[k] * ab_y + g_z_z_0_xy[k];

                g_z_z_y_y[k] = -g_z_z_0_y[k] * ab_y + g_z_z_0_yy[k];

                g_z_z_y_z[k] = -g_z_z_0_z[k] * ab_y + g_z_z_0_yz[k];
            }

            /// Set up 78-81 components of targeted buffer : cbuffer.data(

            auto g_z_z_z_x = cbuffer.data(pp_geom_11_off + 78 * ccomps * dcomps);

            auto g_z_z_z_y = cbuffer.data(pp_geom_11_off + 79 * ccomps * dcomps);

            auto g_z_z_z_z = cbuffer.data(pp_geom_11_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_z_0_0_x, g_z_0_0_y, g_z_0_0_z, g_z_z_0_x, g_z_z_0_xz, g_z_z_0_y, g_z_z_0_yz, g_z_z_0_z, g_z_z_0_zz, g_z_z_z_x, g_z_z_z_y, g_z_z_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_z_x[k] = -g_0_z_0_x[k] + g_z_0_0_x[k] - g_z_z_0_x[k] * ab_z + g_z_z_0_xz[k];

                g_z_z_z_y[k] = -g_0_z_0_y[k] + g_z_0_0_y[k] - g_z_z_0_y[k] * ab_z + g_z_z_0_yz[k];

                g_z_z_z_z[k] = -g_0_z_0_z[k] + g_z_0_0_z[k] - g_z_z_0_z[k] * ab_z + g_z_z_0_zz[k];
            }
        }
    }
}

} // erirec namespace

