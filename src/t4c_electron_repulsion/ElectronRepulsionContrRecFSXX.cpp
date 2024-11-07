#include "ElectronRepulsionContrRecFSXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_hrr_electron_repulsion_fsxx(CSimdArray<double>& cbuffer,
                                     const size_t idx_fsxx,
                                     const size_t idx_dsxx,
                                     const size_t idx_dpxx,
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

            const auto ds_off = idx_dsxx + i * dcomps + j;

            auto g_xx_0 = cbuffer.data(ds_off + 0 * ccomps * dcomps);

            auto g_xy_0 = cbuffer.data(ds_off + 1 * ccomps * dcomps);

            auto g_xz_0 = cbuffer.data(ds_off + 2 * ccomps * dcomps);

            auto g_yy_0 = cbuffer.data(ds_off + 3 * ccomps * dcomps);

            auto g_yz_0 = cbuffer.data(ds_off + 4 * ccomps * dcomps);

            auto g_zz_0 = cbuffer.data(ds_off + 5 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DPSS

            const auto dp_off = idx_dpxx + i * dcomps + j;

            auto g_xx_x = cbuffer.data(dp_off + 0 * ccomps * dcomps);

            auto g_xy_x = cbuffer.data(dp_off + 3 * ccomps * dcomps);

            auto g_xz_x = cbuffer.data(dp_off + 6 * ccomps * dcomps);

            auto g_yy_x = cbuffer.data(dp_off + 9 * ccomps * dcomps);

            auto g_yy_y = cbuffer.data(dp_off + 10 * ccomps * dcomps);

            auto g_yz_x = cbuffer.data(dp_off + 12 * ccomps * dcomps);

            auto g_yz_y = cbuffer.data(dp_off + 13 * ccomps * dcomps);

            auto g_zz_x = cbuffer.data(dp_off + 15 * ccomps * dcomps);

            auto g_zz_y = cbuffer.data(dp_off + 16 * ccomps * dcomps);

            auto g_zz_z = cbuffer.data(dp_off + 17 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_fsxx

            const auto fs_off = idx_fsxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_xxx_0 = cbuffer.data(fs_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0, g_xx_x, g_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxx_0[k] = -g_xx_0[k] * ab_x + g_xx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_xxy_0 = cbuffer.data(fs_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxy_0, g_xy_0, g_xy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxy_0[k] = -g_xy_0[k] * ab_x + g_xy_x[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_xxz_0 = cbuffer.data(fs_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxz_0, g_xz_0, g_xz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxz_0[k] = -g_xz_0[k] * ab_x + g_xz_x[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_xyy_0 = cbuffer.data(fs_off + 3 * ccomps * dcomps);

            #pragma omp simd aligned(g_xyy_0, g_yy_0, g_yy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xyy_0[k] = -g_yy_0[k] * ab_x + g_yy_x[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_xyz_0 = cbuffer.data(fs_off + 4 * ccomps * dcomps);

            #pragma omp simd aligned(g_xyz_0, g_yz_0, g_yz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xyz_0[k] = -g_yz_0[k] * ab_x + g_yz_x[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_xzz_0 = cbuffer.data(fs_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_xzz_0, g_zz_0, g_zz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xzz_0[k] = -g_zz_0[k] * ab_x + g_zz_x[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_yyy_0 = cbuffer.data(fs_off + 6 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0, g_yy_y, g_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yyy_0[k] = -g_yy_0[k] * ab_y + g_yy_y[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_yyz_0 = cbuffer.data(fs_off + 7 * ccomps * dcomps);

            #pragma omp simd aligned(g_yyz_0, g_yz_0, g_yz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yyz_0[k] = -g_yz_0[k] * ab_y + g_yz_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_yzz_0 = cbuffer.data(fs_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_yzz_0, g_zz_0, g_zz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yzz_0[k] = -g_zz_0[k] * ab_y + g_zz_y[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_zzz_0 = cbuffer.data(fs_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0, g_zz_z, g_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zzz_0[k] = -g_zz_0[k] * ab_z + g_zz_z[k];
            }
        }
    }
}

} // erirec namespace

