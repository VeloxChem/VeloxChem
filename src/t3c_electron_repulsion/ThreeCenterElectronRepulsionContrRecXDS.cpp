#include "ThreeCenterElectronRepulsionContrRecXDS.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_hrr_electron_repulsion_xds(CSimdArray<double>& cbuffer,
                                const size_t idx_xds,
                                const size_t idx_xps,
                                const size_t idx_xpp,
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
        /// Set up components of auxilary buffer : SPS

        const auto ps_off = idx_xps + i * 3;

        auto g_x_0 = cbuffer.data(ps_off + 0);

        auto g_y_0 = cbuffer.data(ps_off + 1);

        auto g_z_0 = cbuffer.data(ps_off + 2);

        /// Set up components of auxilary buffer : SPP

        const auto pp_off = idx_xpp + i * 9;

        auto g_x_x = cbuffer.data(pp_off + 0);

        auto g_y_x = cbuffer.data(pp_off + 3);

        auto g_y_y = cbuffer.data(pp_off + 4);

        auto g_z_x = cbuffer.data(pp_off + 6);

        auto g_z_y = cbuffer.data(pp_off + 7);

        auto g_z_z = cbuffer.data(pp_off + 8);

        /// set up bra offset for contr_buffer_xds

        const auto ds_off = idx_xds + i * 6;

        /// Set up 0-1 components of targeted buffer : cbuffer.data(

        auto g_xx_0 = cbuffer.data(ds_off + 0);

        #pragma omp simd aligned(cd_x, g_x_0, g_x_x, g_xx_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xx_0[k] = -g_x_0[k] * cd_x[k] + g_x_x[k];
        }

        /// Set up 1-2 components of targeted buffer : cbuffer.data(

        auto g_xy_0 = cbuffer.data(ds_off + 1);

        #pragma omp simd aligned(cd_x, g_xy_0, g_y_0, g_y_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xy_0[k] = -g_y_0[k] * cd_x[k] + g_y_x[k];
        }

        /// Set up 2-3 components of targeted buffer : cbuffer.data(

        auto g_xz_0 = cbuffer.data(ds_off + 2);

        #pragma omp simd aligned(cd_x, g_xz_0, g_z_0, g_z_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xz_0[k] = -g_z_0[k] * cd_x[k] + g_z_x[k];
        }

        /// Set up 3-4 components of targeted buffer : cbuffer.data(

        auto g_yy_0 = cbuffer.data(ds_off + 3);

        #pragma omp simd aligned(cd_y, g_y_0, g_y_y, g_yy_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yy_0[k] = -g_y_0[k] * cd_y[k] + g_y_y[k];
        }

        /// Set up 4-5 components of targeted buffer : cbuffer.data(

        auto g_yz_0 = cbuffer.data(ds_off + 4);

        #pragma omp simd aligned(cd_y, g_yz_0, g_z_0, g_z_y  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yz_0[k] = -g_z_0[k] * cd_y[k] + g_z_y[k];
        }

        /// Set up 5-6 components of targeted buffer : cbuffer.data(

        auto g_zz_0 = cbuffer.data(ds_off + 5);

        #pragma omp simd aligned(cd_z, g_z_0, g_z_z, g_zz_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_zz_0[k] = -g_z_0[k] * cd_z[k] + g_z_z[k];
        }
    }
}

} // t3ceri namespace

