#include "ElectronRepulsionGeom0010ContrRecXXDS.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxds(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxds,
                                            const size_t idx_xxps,
                                            const size_t idx_geom_10_xxps,
                                            const size_t idx_geom_10_xxpp,
                                            const CSimdArray<double>& factors,
                                            const size_t idx_cd,
                                            const int a_angmom,
                                            const int b_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_cartesian_components(std::array<int, 1>{a_angmom,});

    const auto bcomps = tensor::number_of_cartesian_components(std::array<int, 1>{b_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
            /// Set up components of auxilary buffer : SSPS

            const auto ps_off = idx_xxps + (i * bcomps + j) * 3;

            auto g_x_0 = cbuffer.data(ps_off + 0);

            auto g_y_0 = cbuffer.data(ps_off + 1);

            auto g_z_0 = cbuffer.data(ps_off + 2);

            /// Set up components of auxilary buffer : SSPS

            const auto ps_geom_10_off = idx_geom_10_xxps + (i * bcomps + j) * 3;

            auto g_x_0_x_0 = cbuffer.data(ps_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_y_0 = cbuffer.data(ps_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_z_0 = cbuffer.data(ps_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_y_0_x_0 = cbuffer.data(ps_geom_10_off + 3 * acomps * bcomps + 0);

            auto g_y_0_y_0 = cbuffer.data(ps_geom_10_off + 3 * acomps * bcomps + 1);

            auto g_y_0_z_0 = cbuffer.data(ps_geom_10_off + 3 * acomps * bcomps + 2);

            auto g_z_0_x_0 = cbuffer.data(ps_geom_10_off + 6 * acomps * bcomps + 0);

            auto g_z_0_y_0 = cbuffer.data(ps_geom_10_off + 6 * acomps * bcomps + 1);

            auto g_z_0_z_0 = cbuffer.data(ps_geom_10_off + 6 * acomps * bcomps + 2);

            /// Set up components of auxilary buffer : SSPP

            const auto pp_geom_10_off = idx_geom_10_xxpp + (i * bcomps + j) * 9;

            auto g_x_0_x_x = cbuffer.data(pp_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_x_y = cbuffer.data(pp_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_x_z = cbuffer.data(pp_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_y_y = cbuffer.data(pp_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_z_y = cbuffer.data(pp_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_z_z = cbuffer.data(pp_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_y_0_x_x = cbuffer.data(pp_geom_10_off + 9 * acomps * bcomps + 0);

            auto g_y_0_y_x = cbuffer.data(pp_geom_10_off + 9 * acomps * bcomps + 3);

            auto g_y_0_y_y = cbuffer.data(pp_geom_10_off + 9 * acomps * bcomps + 4);

            auto g_y_0_y_z = cbuffer.data(pp_geom_10_off + 9 * acomps * bcomps + 5);

            auto g_y_0_z_x = cbuffer.data(pp_geom_10_off + 9 * acomps * bcomps + 6);

            auto g_y_0_z_z = cbuffer.data(pp_geom_10_off + 9 * acomps * bcomps + 8);

            auto g_z_0_x_x = cbuffer.data(pp_geom_10_off + 18 * acomps * bcomps + 0);

            auto g_z_0_y_x = cbuffer.data(pp_geom_10_off + 18 * acomps * bcomps + 3);

            auto g_z_0_y_y = cbuffer.data(pp_geom_10_off + 18 * acomps * bcomps + 4);

            auto g_z_0_z_x = cbuffer.data(pp_geom_10_off + 18 * acomps * bcomps + 6);

            auto g_z_0_z_y = cbuffer.data(pp_geom_10_off + 18 * acomps * bcomps + 7);

            auto g_z_0_z_z = cbuffer.data(pp_geom_10_off + 18 * acomps * bcomps + 8);

            /// set up bra offset for contr_buffer_xxds

            const auto ds_geom_10_off = idx_geom_10_xxds + (i * bcomps + j) * 6;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_x_0_xx_0 = cbuffer.data(ds_geom_10_off + 0 * acomps * bcomps + 0);

            #pragma omp simd aligned(cd_x, g_x_0, g_x_0_x_0, g_x_0_x_x, g_x_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xx_0[k] = -g_x_0[k] - g_x_0_x_0[k] * cd_x[k] + g_x_0_x_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_x_0_xy_0 = cbuffer.data(ds_geom_10_off + 0 * acomps * bcomps + 1);

            #pragma omp simd aligned(cd_y, g_x_0_x_0, g_x_0_x_y, g_x_0_xy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xy_0[k] = -g_x_0_x_0[k] * cd_y[k] + g_x_0_x_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_xz_0 = cbuffer.data(ds_geom_10_off + 0 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_z, g_x_0_x_0, g_x_0_x_z, g_x_0_xz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xz_0[k] = -g_x_0_x_0[k] * cd_z[k] + g_x_0_x_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_x_0_yy_0 = cbuffer.data(ds_geom_10_off + 0 * acomps * bcomps + 3);

            #pragma omp simd aligned(cd_y, g_x_0_y_0, g_x_0_y_y, g_x_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yy_0[k] = -g_x_0_y_0[k] * cd_y[k] + g_x_0_y_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_x_0_yz_0 = cbuffer.data(ds_geom_10_off + 0 * acomps * bcomps + 4);

            #pragma omp simd aligned(cd_y, g_x_0_yz_0, g_x_0_z_0, g_x_0_z_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yz_0[k] = -g_x_0_z_0[k] * cd_y[k] + g_x_0_z_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_zz_0 = cbuffer.data(ds_geom_10_off + 0 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_z, g_x_0_z_0, g_x_0_z_z, g_x_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zz_0[k] = -g_x_0_z_0[k] * cd_z[k] + g_x_0_z_z[k];
            }
            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_y_0_xx_0 = cbuffer.data(ds_geom_10_off + 6 * acomps * bcomps + 0);

            #pragma omp simd aligned(cd_x, g_y_0_x_0, g_y_0_x_x, g_y_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xx_0[k] = -g_y_0_x_0[k] * cd_x[k] + g_y_0_x_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_y_0_xy_0 = cbuffer.data(ds_geom_10_off + 6 * acomps * bcomps + 1);

            #pragma omp simd aligned(cd_x, g_y_0_xy_0, g_y_0_y_0, g_y_0_y_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xy_0[k] = -g_y_0_y_0[k] * cd_x[k] + g_y_0_y_x[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_y_0_xz_0 = cbuffer.data(ds_geom_10_off + 6 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_x, g_y_0_xz_0, g_y_0_z_0, g_y_0_z_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xz_0[k] = -g_y_0_z_0[k] * cd_x[k] + g_y_0_z_x[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_y_0_yy_0 = cbuffer.data(ds_geom_10_off + 6 * acomps * bcomps + 3);

            #pragma omp simd aligned(cd_y, g_y_0, g_y_0_y_0, g_y_0_y_y, g_y_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yy_0[k] = -g_y_0[k] - g_y_0_y_0[k] * cd_y[k] + g_y_0_y_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_y_0_yz_0 = cbuffer.data(ds_geom_10_off + 6 * acomps * bcomps + 4);

            #pragma omp simd aligned(cd_z, g_y_0_y_0, g_y_0_y_z, g_y_0_yz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yz_0[k] = -g_y_0_y_0[k] * cd_z[k] + g_y_0_y_z[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_y_0_zz_0 = cbuffer.data(ds_geom_10_off + 6 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_z, g_y_0_z_0, g_y_0_z_z, g_y_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zz_0[k] = -g_y_0_z_0[k] * cd_z[k] + g_y_0_z_z[k];
            }
            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_z_0_xx_0 = cbuffer.data(ds_geom_10_off + 12 * acomps * bcomps + 0);

            #pragma omp simd aligned(cd_x, g_z_0_x_0, g_z_0_x_x, g_z_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xx_0[k] = -g_z_0_x_0[k] * cd_x[k] + g_z_0_x_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_z_0_xy_0 = cbuffer.data(ds_geom_10_off + 12 * acomps * bcomps + 1);

            #pragma omp simd aligned(cd_x, g_z_0_xy_0, g_z_0_y_0, g_z_0_y_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xy_0[k] = -g_z_0_y_0[k] * cd_x[k] + g_z_0_y_x[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_z_0_xz_0 = cbuffer.data(ds_geom_10_off + 12 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_x, g_z_0_xz_0, g_z_0_z_0, g_z_0_z_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xz_0[k] = -g_z_0_z_0[k] * cd_x[k] + g_z_0_z_x[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_z_0_yy_0 = cbuffer.data(ds_geom_10_off + 12 * acomps * bcomps + 3);

            #pragma omp simd aligned(cd_y, g_z_0_y_0, g_z_0_y_y, g_z_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yy_0[k] = -g_z_0_y_0[k] * cd_y[k] + g_z_0_y_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_z_0_yz_0 = cbuffer.data(ds_geom_10_off + 12 * acomps * bcomps + 4);

            #pragma omp simd aligned(cd_y, g_z_0_yz_0, g_z_0_z_0, g_z_0_z_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yz_0[k] = -g_z_0_z_0[k] * cd_y[k] + g_z_0_z_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_z_0_zz_0 = cbuffer.data(ds_geom_10_off + 12 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_z, g_z_0, g_z_0_z_0, g_z_0_z_z, g_z_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zz_0[k] = -g_z_0[k] - g_z_0_z_0[k] * cd_z[k] + g_z_0_z_z[k];
            }
        }
    }
}

} // erirec namespace

