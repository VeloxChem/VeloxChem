#include "ElectronRepulsionGeom1010ContrRecDSXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom1010_hrr_electron_repulsion_dsxx(CSimdArray<double>& cbuffer,
                                              const size_t idx_geom_1010_dsxx,
                                              const size_t idx_geom_0010_psxx,
                                              const size_t idx_geom_1010_psxx,
                                              const size_t idx_geom_1010_ppxx,
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
            /// Set up components of auxilary buffer : PSSS

            const auto ps_geom_0010_off = idx_geom_0010_psxx + i * dcomps + j;

            auto g_0_0_x_0_x_0 = cbuffer.data(ps_geom_0010_off + 0 * ccomps * dcomps);

            auto g_0_0_x_0_y_0 = cbuffer.data(ps_geom_0010_off + 1 * ccomps * dcomps);

            auto g_0_0_x_0_z_0 = cbuffer.data(ps_geom_0010_off + 2 * ccomps * dcomps);

            auto g_0_0_y_0_x_0 = cbuffer.data(ps_geom_0010_off + 3 * ccomps * dcomps);

            auto g_0_0_y_0_y_0 = cbuffer.data(ps_geom_0010_off + 4 * ccomps * dcomps);

            auto g_0_0_y_0_z_0 = cbuffer.data(ps_geom_0010_off + 5 * ccomps * dcomps);

            auto g_0_0_z_0_x_0 = cbuffer.data(ps_geom_0010_off + 6 * ccomps * dcomps);

            auto g_0_0_z_0_y_0 = cbuffer.data(ps_geom_0010_off + 7 * ccomps * dcomps);

            auto g_0_0_z_0_z_0 = cbuffer.data(ps_geom_0010_off + 8 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PSSS

            const auto ps_geom_1010_off = idx_geom_1010_psxx + i * dcomps + j;

            auto g_x_0_x_0_x_0 = cbuffer.data(ps_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_y_0 = cbuffer.data(ps_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_z_0 = cbuffer.data(ps_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_y_0_x_0 = cbuffer.data(ps_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_y_0_y_0 = cbuffer.data(ps_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_y_0_z_0 = cbuffer.data(ps_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_z_0_x_0 = cbuffer.data(ps_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_z_0_y_0 = cbuffer.data(ps_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_z_0_z_0 = cbuffer.data(ps_geom_1010_off + 8 * ccomps * dcomps);

            auto g_y_0_x_0_x_0 = cbuffer.data(ps_geom_1010_off + 9 * ccomps * dcomps);

            auto g_y_0_x_0_y_0 = cbuffer.data(ps_geom_1010_off + 10 * ccomps * dcomps);

            auto g_y_0_x_0_z_0 = cbuffer.data(ps_geom_1010_off + 11 * ccomps * dcomps);

            auto g_y_0_y_0_x_0 = cbuffer.data(ps_geom_1010_off + 12 * ccomps * dcomps);

            auto g_y_0_y_0_y_0 = cbuffer.data(ps_geom_1010_off + 13 * ccomps * dcomps);

            auto g_y_0_y_0_z_0 = cbuffer.data(ps_geom_1010_off + 14 * ccomps * dcomps);

            auto g_y_0_z_0_x_0 = cbuffer.data(ps_geom_1010_off + 15 * ccomps * dcomps);

            auto g_y_0_z_0_y_0 = cbuffer.data(ps_geom_1010_off + 16 * ccomps * dcomps);

            auto g_y_0_z_0_z_0 = cbuffer.data(ps_geom_1010_off + 17 * ccomps * dcomps);

            auto g_z_0_x_0_x_0 = cbuffer.data(ps_geom_1010_off + 18 * ccomps * dcomps);

            auto g_z_0_x_0_y_0 = cbuffer.data(ps_geom_1010_off + 19 * ccomps * dcomps);

            auto g_z_0_x_0_z_0 = cbuffer.data(ps_geom_1010_off + 20 * ccomps * dcomps);

            auto g_z_0_y_0_x_0 = cbuffer.data(ps_geom_1010_off + 21 * ccomps * dcomps);

            auto g_z_0_y_0_y_0 = cbuffer.data(ps_geom_1010_off + 22 * ccomps * dcomps);

            auto g_z_0_y_0_z_0 = cbuffer.data(ps_geom_1010_off + 23 * ccomps * dcomps);

            auto g_z_0_z_0_x_0 = cbuffer.data(ps_geom_1010_off + 24 * ccomps * dcomps);

            auto g_z_0_z_0_y_0 = cbuffer.data(ps_geom_1010_off + 25 * ccomps * dcomps);

            auto g_z_0_z_0_z_0 = cbuffer.data(ps_geom_1010_off + 26 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PPSS

            const auto pp_geom_1010_off = idx_geom_1010_ppxx + i * dcomps + j;

            auto g_x_0_x_0_x_x = cbuffer.data(pp_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_x_y = cbuffer.data(pp_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_x_z = cbuffer.data(pp_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_y_x = cbuffer.data(pp_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_y_y = cbuffer.data(pp_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_y_z = cbuffer.data(pp_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_x_0_z_x = cbuffer.data(pp_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_z_y = cbuffer.data(pp_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_z_z = cbuffer.data(pp_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_y_0_x_x = cbuffer.data(pp_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_y_0_x_y = cbuffer.data(pp_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_y_0_x_z = cbuffer.data(pp_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_y_0_y_x = cbuffer.data(pp_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_y_0_y_y = cbuffer.data(pp_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_y_0_y_z = cbuffer.data(pp_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_y_0_z_x = cbuffer.data(pp_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_y_0_z_y = cbuffer.data(pp_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_y_0_z_z = cbuffer.data(pp_geom_1010_off + 17 * ccomps * dcomps);

            auto g_x_0_z_0_x_x = cbuffer.data(pp_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_z_0_x_y = cbuffer.data(pp_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_z_0_x_z = cbuffer.data(pp_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_z_0_y_x = cbuffer.data(pp_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_z_0_y_y = cbuffer.data(pp_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_z_0_y_z = cbuffer.data(pp_geom_1010_off + 23 * ccomps * dcomps);

            auto g_x_0_z_0_z_x = cbuffer.data(pp_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_z_0_z_y = cbuffer.data(pp_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_z_0_z_z = cbuffer.data(pp_geom_1010_off + 26 * ccomps * dcomps);

            auto g_y_0_x_0_x_x = cbuffer.data(pp_geom_1010_off + 27 * ccomps * dcomps);

            auto g_y_0_x_0_x_y = cbuffer.data(pp_geom_1010_off + 28 * ccomps * dcomps);

            auto g_y_0_x_0_x_z = cbuffer.data(pp_geom_1010_off + 29 * ccomps * dcomps);

            auto g_y_0_x_0_y_x = cbuffer.data(pp_geom_1010_off + 30 * ccomps * dcomps);

            auto g_y_0_x_0_y_y = cbuffer.data(pp_geom_1010_off + 31 * ccomps * dcomps);

            auto g_y_0_x_0_y_z = cbuffer.data(pp_geom_1010_off + 32 * ccomps * dcomps);

            auto g_y_0_x_0_z_x = cbuffer.data(pp_geom_1010_off + 33 * ccomps * dcomps);

            auto g_y_0_x_0_z_y = cbuffer.data(pp_geom_1010_off + 34 * ccomps * dcomps);

            auto g_y_0_x_0_z_z = cbuffer.data(pp_geom_1010_off + 35 * ccomps * dcomps);

            auto g_y_0_y_0_x_x = cbuffer.data(pp_geom_1010_off + 36 * ccomps * dcomps);

            auto g_y_0_y_0_x_y = cbuffer.data(pp_geom_1010_off + 37 * ccomps * dcomps);

            auto g_y_0_y_0_x_z = cbuffer.data(pp_geom_1010_off + 38 * ccomps * dcomps);

            auto g_y_0_y_0_y_x = cbuffer.data(pp_geom_1010_off + 39 * ccomps * dcomps);

            auto g_y_0_y_0_y_y = cbuffer.data(pp_geom_1010_off + 40 * ccomps * dcomps);

            auto g_y_0_y_0_y_z = cbuffer.data(pp_geom_1010_off + 41 * ccomps * dcomps);

            auto g_y_0_y_0_z_x = cbuffer.data(pp_geom_1010_off + 42 * ccomps * dcomps);

            auto g_y_0_y_0_z_y = cbuffer.data(pp_geom_1010_off + 43 * ccomps * dcomps);

            auto g_y_0_y_0_z_z = cbuffer.data(pp_geom_1010_off + 44 * ccomps * dcomps);

            auto g_y_0_z_0_x_x = cbuffer.data(pp_geom_1010_off + 45 * ccomps * dcomps);

            auto g_y_0_z_0_x_y = cbuffer.data(pp_geom_1010_off + 46 * ccomps * dcomps);

            auto g_y_0_z_0_x_z = cbuffer.data(pp_geom_1010_off + 47 * ccomps * dcomps);

            auto g_y_0_z_0_y_x = cbuffer.data(pp_geom_1010_off + 48 * ccomps * dcomps);

            auto g_y_0_z_0_y_y = cbuffer.data(pp_geom_1010_off + 49 * ccomps * dcomps);

            auto g_y_0_z_0_y_z = cbuffer.data(pp_geom_1010_off + 50 * ccomps * dcomps);

            auto g_y_0_z_0_z_x = cbuffer.data(pp_geom_1010_off + 51 * ccomps * dcomps);

            auto g_y_0_z_0_z_y = cbuffer.data(pp_geom_1010_off + 52 * ccomps * dcomps);

            auto g_y_0_z_0_z_z = cbuffer.data(pp_geom_1010_off + 53 * ccomps * dcomps);

            auto g_z_0_x_0_x_x = cbuffer.data(pp_geom_1010_off + 54 * ccomps * dcomps);

            auto g_z_0_x_0_x_y = cbuffer.data(pp_geom_1010_off + 55 * ccomps * dcomps);

            auto g_z_0_x_0_x_z = cbuffer.data(pp_geom_1010_off + 56 * ccomps * dcomps);

            auto g_z_0_x_0_y_x = cbuffer.data(pp_geom_1010_off + 57 * ccomps * dcomps);

            auto g_z_0_x_0_y_y = cbuffer.data(pp_geom_1010_off + 58 * ccomps * dcomps);

            auto g_z_0_x_0_y_z = cbuffer.data(pp_geom_1010_off + 59 * ccomps * dcomps);

            auto g_z_0_x_0_z_x = cbuffer.data(pp_geom_1010_off + 60 * ccomps * dcomps);

            auto g_z_0_x_0_z_y = cbuffer.data(pp_geom_1010_off + 61 * ccomps * dcomps);

            auto g_z_0_x_0_z_z = cbuffer.data(pp_geom_1010_off + 62 * ccomps * dcomps);

            auto g_z_0_y_0_x_x = cbuffer.data(pp_geom_1010_off + 63 * ccomps * dcomps);

            auto g_z_0_y_0_x_y = cbuffer.data(pp_geom_1010_off + 64 * ccomps * dcomps);

            auto g_z_0_y_0_x_z = cbuffer.data(pp_geom_1010_off + 65 * ccomps * dcomps);

            auto g_z_0_y_0_y_x = cbuffer.data(pp_geom_1010_off + 66 * ccomps * dcomps);

            auto g_z_0_y_0_y_y = cbuffer.data(pp_geom_1010_off + 67 * ccomps * dcomps);

            auto g_z_0_y_0_y_z = cbuffer.data(pp_geom_1010_off + 68 * ccomps * dcomps);

            auto g_z_0_y_0_z_x = cbuffer.data(pp_geom_1010_off + 69 * ccomps * dcomps);

            auto g_z_0_y_0_z_y = cbuffer.data(pp_geom_1010_off + 70 * ccomps * dcomps);

            auto g_z_0_y_0_z_z = cbuffer.data(pp_geom_1010_off + 71 * ccomps * dcomps);

            auto g_z_0_z_0_x_x = cbuffer.data(pp_geom_1010_off + 72 * ccomps * dcomps);

            auto g_z_0_z_0_x_y = cbuffer.data(pp_geom_1010_off + 73 * ccomps * dcomps);

            auto g_z_0_z_0_x_z = cbuffer.data(pp_geom_1010_off + 74 * ccomps * dcomps);

            auto g_z_0_z_0_y_x = cbuffer.data(pp_geom_1010_off + 75 * ccomps * dcomps);

            auto g_z_0_z_0_y_y = cbuffer.data(pp_geom_1010_off + 76 * ccomps * dcomps);

            auto g_z_0_z_0_y_z = cbuffer.data(pp_geom_1010_off + 77 * ccomps * dcomps);

            auto g_z_0_z_0_z_x = cbuffer.data(pp_geom_1010_off + 78 * ccomps * dcomps);

            auto g_z_0_z_0_z_y = cbuffer.data(pp_geom_1010_off + 79 * ccomps * dcomps);

            auto g_z_0_z_0_z_z = cbuffer.data(pp_geom_1010_off + 80 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_dsxx

            const auto ds_geom_1010_off = idx_geom_1010_dsxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xx_0 = cbuffer.data(ds_geom_1010_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_x_0, g_x_0_x_0_x_0, g_x_0_x_0_x_x, g_x_0_x_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xx_0[k] = -g_0_0_x_0_x_0[k] - g_x_0_x_0_x_0[k] * ab_x + g_x_0_x_0_x_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xy_0 = cbuffer.data(ds_geom_1010_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_x_0, g_x_0_x_0_x_y, g_x_0_x_0_xy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xy_0[k] = -g_x_0_x_0_x_0[k] * ab_y + g_x_0_x_0_x_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_xz_0 = cbuffer.data(ds_geom_1010_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_x_0, g_x_0_x_0_x_z, g_x_0_x_0_xz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_xz_0[k] = -g_x_0_x_0_x_0[k] * ab_z + g_x_0_x_0_x_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yy_0 = cbuffer.data(ds_geom_1010_off + 3 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_y_0, g_x_0_x_0_y_y, g_x_0_x_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yy_0[k] = -g_x_0_x_0_y_0[k] * ab_y + g_x_0_x_0_y_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_yz_0 = cbuffer.data(ds_geom_1010_off + 4 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_yz_0, g_x_0_x_0_z_0, g_x_0_x_0_z_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_yz_0[k] = -g_x_0_x_0_z_0[k] * ab_y + g_x_0_x_0_z_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_zz_0 = cbuffer.data(ds_geom_1010_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_z_0, g_x_0_x_0_z_z, g_x_0_x_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_zz_0[k] = -g_x_0_x_0_z_0[k] * ab_z + g_x_0_x_0_z_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xx_0 = cbuffer.data(ds_geom_1010_off + 6 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_x_0, g_x_0_y_0_x_0, g_x_0_y_0_x_x, g_x_0_y_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xx_0[k] = -g_0_0_y_0_x_0[k] - g_x_0_y_0_x_0[k] * ab_x + g_x_0_y_0_x_x[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xy_0 = cbuffer.data(ds_geom_1010_off + 7 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_x_0, g_x_0_y_0_x_y, g_x_0_y_0_xy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xy_0[k] = -g_x_0_y_0_x_0[k] * ab_y + g_x_0_y_0_x_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_xz_0 = cbuffer.data(ds_geom_1010_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_x_0, g_x_0_y_0_x_z, g_x_0_y_0_xz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_xz_0[k] = -g_x_0_y_0_x_0[k] * ab_z + g_x_0_y_0_x_z[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yy_0 = cbuffer.data(ds_geom_1010_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_y_0, g_x_0_y_0_y_y, g_x_0_y_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yy_0[k] = -g_x_0_y_0_y_0[k] * ab_y + g_x_0_y_0_y_y[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_yz_0 = cbuffer.data(ds_geom_1010_off + 10 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_yz_0, g_x_0_y_0_z_0, g_x_0_y_0_z_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_yz_0[k] = -g_x_0_y_0_z_0[k] * ab_y + g_x_0_y_0_z_y[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_zz_0 = cbuffer.data(ds_geom_1010_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_z_0, g_x_0_y_0_z_z, g_x_0_y_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_zz_0[k] = -g_x_0_y_0_z_0[k] * ab_z + g_x_0_y_0_z_z[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xx_0 = cbuffer.data(ds_geom_1010_off + 12 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_x_0, g_x_0_z_0_x_0, g_x_0_z_0_x_x, g_x_0_z_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xx_0[k] = -g_0_0_z_0_x_0[k] - g_x_0_z_0_x_0[k] * ab_x + g_x_0_z_0_x_x[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xy_0 = cbuffer.data(ds_geom_1010_off + 13 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_x_0, g_x_0_z_0_x_y, g_x_0_z_0_xy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xy_0[k] = -g_x_0_z_0_x_0[k] * ab_y + g_x_0_z_0_x_y[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_xz_0 = cbuffer.data(ds_geom_1010_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_x_0, g_x_0_z_0_x_z, g_x_0_z_0_xz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_xz_0[k] = -g_x_0_z_0_x_0[k] * ab_z + g_x_0_z_0_x_z[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yy_0 = cbuffer.data(ds_geom_1010_off + 15 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_y_0, g_x_0_z_0_y_y, g_x_0_z_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yy_0[k] = -g_x_0_z_0_y_0[k] * ab_y + g_x_0_z_0_y_y[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_yz_0 = cbuffer.data(ds_geom_1010_off + 16 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_yz_0, g_x_0_z_0_z_0, g_x_0_z_0_z_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_yz_0[k] = -g_x_0_z_0_z_0[k] * ab_y + g_x_0_z_0_z_y[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_zz_0 = cbuffer.data(ds_geom_1010_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_z_0, g_x_0_z_0_z_z, g_x_0_z_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_zz_0[k] = -g_x_0_z_0_z_0[k] * ab_z + g_x_0_z_0_z_z[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xx_0 = cbuffer.data(ds_geom_1010_off + 18 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_x_0, g_y_0_x_0_x_x, g_y_0_x_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xx_0[k] = -g_y_0_x_0_x_0[k] * ab_x + g_y_0_x_0_x_x[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xy_0 = cbuffer.data(ds_geom_1010_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xy_0, g_y_0_x_0_y_0, g_y_0_x_0_y_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xy_0[k] = -g_y_0_x_0_y_0[k] * ab_x + g_y_0_x_0_y_x[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_xz_0 = cbuffer.data(ds_geom_1010_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_xz_0, g_y_0_x_0_z_0, g_y_0_x_0_z_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_xz_0[k] = -g_y_0_x_0_z_0[k] * ab_x + g_y_0_x_0_z_x[k];
            }

            /// Set up 21-22 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yy_0 = cbuffer.data(ds_geom_1010_off + 21 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_y_0, g_y_0_x_0_y_0, g_y_0_x_0_y_y, g_y_0_x_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yy_0[k] = -g_0_0_x_0_y_0[k] - g_y_0_x_0_y_0[k] * ab_y + g_y_0_x_0_y_y[k];
            }

            /// Set up 22-23 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_yz_0 = cbuffer.data(ds_geom_1010_off + 22 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_y_0, g_y_0_x_0_y_z, g_y_0_x_0_yz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_yz_0[k] = -g_y_0_x_0_y_0[k] * ab_z + g_y_0_x_0_y_z[k];
            }

            /// Set up 23-24 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_zz_0 = cbuffer.data(ds_geom_1010_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_z_0, g_y_0_x_0_z_z, g_y_0_x_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_zz_0[k] = -g_y_0_x_0_z_0[k] * ab_z + g_y_0_x_0_z_z[k];
            }

            /// Set up 24-25 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xx_0 = cbuffer.data(ds_geom_1010_off + 24 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_x_0, g_y_0_y_0_x_x, g_y_0_y_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xx_0[k] = -g_y_0_y_0_x_0[k] * ab_x + g_y_0_y_0_x_x[k];
            }

            /// Set up 25-26 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xy_0 = cbuffer.data(ds_geom_1010_off + 25 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xy_0, g_y_0_y_0_y_0, g_y_0_y_0_y_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xy_0[k] = -g_y_0_y_0_y_0[k] * ab_x + g_y_0_y_0_y_x[k];
            }

            /// Set up 26-27 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_xz_0 = cbuffer.data(ds_geom_1010_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_xz_0, g_y_0_y_0_z_0, g_y_0_y_0_z_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_xz_0[k] = -g_y_0_y_0_z_0[k] * ab_x + g_y_0_y_0_z_x[k];
            }

            /// Set up 27-28 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yy_0 = cbuffer.data(ds_geom_1010_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_y_0, g_y_0_y_0_y_0, g_y_0_y_0_y_y, g_y_0_y_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yy_0[k] = -g_0_0_y_0_y_0[k] - g_y_0_y_0_y_0[k] * ab_y + g_y_0_y_0_y_y[k];
            }

            /// Set up 28-29 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_yz_0 = cbuffer.data(ds_geom_1010_off + 28 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_y_0, g_y_0_y_0_y_z, g_y_0_y_0_yz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_yz_0[k] = -g_y_0_y_0_y_0[k] * ab_z + g_y_0_y_0_y_z[k];
            }

            /// Set up 29-30 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_zz_0 = cbuffer.data(ds_geom_1010_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_z_0, g_y_0_y_0_z_z, g_y_0_y_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_zz_0[k] = -g_y_0_y_0_z_0[k] * ab_z + g_y_0_y_0_z_z[k];
            }

            /// Set up 30-31 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xx_0 = cbuffer.data(ds_geom_1010_off + 30 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_x_0, g_y_0_z_0_x_x, g_y_0_z_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xx_0[k] = -g_y_0_z_0_x_0[k] * ab_x + g_y_0_z_0_x_x[k];
            }

            /// Set up 31-32 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xy_0 = cbuffer.data(ds_geom_1010_off + 31 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xy_0, g_y_0_z_0_y_0, g_y_0_z_0_y_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xy_0[k] = -g_y_0_z_0_y_0[k] * ab_x + g_y_0_z_0_y_x[k];
            }

            /// Set up 32-33 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_xz_0 = cbuffer.data(ds_geom_1010_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_xz_0, g_y_0_z_0_z_0, g_y_0_z_0_z_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_xz_0[k] = -g_y_0_z_0_z_0[k] * ab_x + g_y_0_z_0_z_x[k];
            }

            /// Set up 33-34 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yy_0 = cbuffer.data(ds_geom_1010_off + 33 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_y_0, g_y_0_z_0_y_0, g_y_0_z_0_y_y, g_y_0_z_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yy_0[k] = -g_0_0_z_0_y_0[k] - g_y_0_z_0_y_0[k] * ab_y + g_y_0_z_0_y_y[k];
            }

            /// Set up 34-35 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_yz_0 = cbuffer.data(ds_geom_1010_off + 34 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_y_0, g_y_0_z_0_y_z, g_y_0_z_0_yz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_yz_0[k] = -g_y_0_z_0_y_0[k] * ab_z + g_y_0_z_0_y_z[k];
            }

            /// Set up 35-36 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_zz_0 = cbuffer.data(ds_geom_1010_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_z_0, g_y_0_z_0_z_z, g_y_0_z_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_zz_0[k] = -g_y_0_z_0_z_0[k] * ab_z + g_y_0_z_0_z_z[k];
            }

            /// Set up 36-37 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xx_0 = cbuffer.data(ds_geom_1010_off + 36 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_x_0, g_z_0_x_0_x_x, g_z_0_x_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xx_0[k] = -g_z_0_x_0_x_0[k] * ab_x + g_z_0_x_0_x_x[k];
            }

            /// Set up 37-38 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xy_0 = cbuffer.data(ds_geom_1010_off + 37 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xy_0, g_z_0_x_0_y_0, g_z_0_x_0_y_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xy_0[k] = -g_z_0_x_0_y_0[k] * ab_x + g_z_0_x_0_y_x[k];
            }

            /// Set up 38-39 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_xz_0 = cbuffer.data(ds_geom_1010_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_xz_0, g_z_0_x_0_z_0, g_z_0_x_0_z_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_xz_0[k] = -g_z_0_x_0_z_0[k] * ab_x + g_z_0_x_0_z_x[k];
            }

            /// Set up 39-40 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yy_0 = cbuffer.data(ds_geom_1010_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_y_0, g_z_0_x_0_y_y, g_z_0_x_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yy_0[k] = -g_z_0_x_0_y_0[k] * ab_y + g_z_0_x_0_y_y[k];
            }

            /// Set up 40-41 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_yz_0 = cbuffer.data(ds_geom_1010_off + 40 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_yz_0, g_z_0_x_0_z_0, g_z_0_x_0_z_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_yz_0[k] = -g_z_0_x_0_z_0[k] * ab_y + g_z_0_x_0_z_y[k];
            }

            /// Set up 41-42 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_zz_0 = cbuffer.data(ds_geom_1010_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_z_0, g_z_0_x_0_z_0, g_z_0_x_0_z_z, g_z_0_x_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_zz_0[k] = -g_0_0_x_0_z_0[k] - g_z_0_x_0_z_0[k] * ab_z + g_z_0_x_0_z_z[k];
            }

            /// Set up 42-43 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xx_0 = cbuffer.data(ds_geom_1010_off + 42 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_x_0, g_z_0_y_0_x_x, g_z_0_y_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xx_0[k] = -g_z_0_y_0_x_0[k] * ab_x + g_z_0_y_0_x_x[k];
            }

            /// Set up 43-44 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xy_0 = cbuffer.data(ds_geom_1010_off + 43 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xy_0, g_z_0_y_0_y_0, g_z_0_y_0_y_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xy_0[k] = -g_z_0_y_0_y_0[k] * ab_x + g_z_0_y_0_y_x[k];
            }

            /// Set up 44-45 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_xz_0 = cbuffer.data(ds_geom_1010_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_xz_0, g_z_0_y_0_z_0, g_z_0_y_0_z_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_xz_0[k] = -g_z_0_y_0_z_0[k] * ab_x + g_z_0_y_0_z_x[k];
            }

            /// Set up 45-46 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yy_0 = cbuffer.data(ds_geom_1010_off + 45 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_y_0, g_z_0_y_0_y_y, g_z_0_y_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yy_0[k] = -g_z_0_y_0_y_0[k] * ab_y + g_z_0_y_0_y_y[k];
            }

            /// Set up 46-47 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_yz_0 = cbuffer.data(ds_geom_1010_off + 46 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_yz_0, g_z_0_y_0_z_0, g_z_0_y_0_z_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_yz_0[k] = -g_z_0_y_0_z_0[k] * ab_y + g_z_0_y_0_z_y[k];
            }

            /// Set up 47-48 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_zz_0 = cbuffer.data(ds_geom_1010_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_z_0, g_z_0_y_0_z_0, g_z_0_y_0_z_z, g_z_0_y_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_zz_0[k] = -g_0_0_y_0_z_0[k] - g_z_0_y_0_z_0[k] * ab_z + g_z_0_y_0_z_z[k];
            }

            /// Set up 48-49 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xx_0 = cbuffer.data(ds_geom_1010_off + 48 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_x_0, g_z_0_z_0_x_x, g_z_0_z_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xx_0[k] = -g_z_0_z_0_x_0[k] * ab_x + g_z_0_z_0_x_x[k];
            }

            /// Set up 49-50 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xy_0 = cbuffer.data(ds_geom_1010_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xy_0, g_z_0_z_0_y_0, g_z_0_z_0_y_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xy_0[k] = -g_z_0_z_0_y_0[k] * ab_x + g_z_0_z_0_y_x[k];
            }

            /// Set up 50-51 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_xz_0 = cbuffer.data(ds_geom_1010_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_xz_0, g_z_0_z_0_z_0, g_z_0_z_0_z_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_xz_0[k] = -g_z_0_z_0_z_0[k] * ab_x + g_z_0_z_0_z_x[k];
            }

            /// Set up 51-52 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yy_0 = cbuffer.data(ds_geom_1010_off + 51 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_y_0, g_z_0_z_0_y_y, g_z_0_z_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yy_0[k] = -g_z_0_z_0_y_0[k] * ab_y + g_z_0_z_0_y_y[k];
            }

            /// Set up 52-53 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_yz_0 = cbuffer.data(ds_geom_1010_off + 52 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_yz_0, g_z_0_z_0_z_0, g_z_0_z_0_z_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_yz_0[k] = -g_z_0_z_0_z_0[k] * ab_y + g_z_0_z_0_z_y[k];
            }

            /// Set up 53-54 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_zz_0 = cbuffer.data(ds_geom_1010_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_z_0, g_z_0_z_0_z_0, g_z_0_z_0_z_z, g_z_0_z_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_zz_0[k] = -g_0_0_z_0_z_0[k] - g_z_0_z_0_z_0[k] * ab_z + g_z_0_z_0_z_z[k];
            }
        }
    }
}

} // erirec namespace

