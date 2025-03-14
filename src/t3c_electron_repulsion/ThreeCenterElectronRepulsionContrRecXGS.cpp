#include "ThreeCenterElectronRepulsionContrRecXGS.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_hrr_electron_repulsion_xgs(CSimdArray<double>& cbuffer,
                                const size_t idx_xgs,
                                const size_t idx_xfs,
                                const size_t idx_xfp,
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
        /// Set up components of auxilary buffer : SFS

        const auto fs_off = idx_xfs + i * 10;

        auto g_xxx_0 = cbuffer.data(fs_off + 0);

        auto g_xxy_0 = cbuffer.data(fs_off + 1);

        auto g_xxz_0 = cbuffer.data(fs_off + 2);

        auto g_xyy_0 = cbuffer.data(fs_off + 3);

        auto g_xyz_0 = cbuffer.data(fs_off + 4);

        auto g_xzz_0 = cbuffer.data(fs_off + 5);

        auto g_yyy_0 = cbuffer.data(fs_off + 6);

        auto g_yyz_0 = cbuffer.data(fs_off + 7);

        auto g_yzz_0 = cbuffer.data(fs_off + 8);

        auto g_zzz_0 = cbuffer.data(fs_off + 9);

        /// Set up components of auxilary buffer : SFP

        const auto fp_off = idx_xfp + i * 30;

        auto g_xxx_x = cbuffer.data(fp_off + 0);

        auto g_xxy_x = cbuffer.data(fp_off + 3);

        auto g_xxz_x = cbuffer.data(fp_off + 6);

        auto g_xyy_x = cbuffer.data(fp_off + 9);

        auto g_xyz_x = cbuffer.data(fp_off + 12);

        auto g_xzz_x = cbuffer.data(fp_off + 15);

        auto g_yyy_x = cbuffer.data(fp_off + 18);

        auto g_yyy_y = cbuffer.data(fp_off + 19);

        auto g_yyz_x = cbuffer.data(fp_off + 21);

        auto g_yyz_y = cbuffer.data(fp_off + 22);

        auto g_yzz_x = cbuffer.data(fp_off + 24);

        auto g_yzz_y = cbuffer.data(fp_off + 25);

        auto g_zzz_x = cbuffer.data(fp_off + 27);

        auto g_zzz_y = cbuffer.data(fp_off + 28);

        auto g_zzz_z = cbuffer.data(fp_off + 29);

        /// set up bra offset for contr_buffer_xgs

        const auto gs_off = idx_xgs + i * 15;

        /// Set up 0-1 components of targeted buffer : cbuffer.data(

        auto g_xxxx_0 = cbuffer.data(gs_off + 0);

        #pragma omp simd aligned(cd_x, g_xxx_0, g_xxx_x, g_xxxx_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxxx_0[k] = -g_xxx_0[k] * cd_x[k] + g_xxx_x[k];
        }

        /// Set up 1-2 components of targeted buffer : cbuffer.data(

        auto g_xxxy_0 = cbuffer.data(gs_off + 1);

        #pragma omp simd aligned(cd_x, g_xxxy_0, g_xxy_0, g_xxy_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxxy_0[k] = -g_xxy_0[k] * cd_x[k] + g_xxy_x[k];
        }

        /// Set up 2-3 components of targeted buffer : cbuffer.data(

        auto g_xxxz_0 = cbuffer.data(gs_off + 2);

        #pragma omp simd aligned(cd_x, g_xxxz_0, g_xxz_0, g_xxz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxxz_0[k] = -g_xxz_0[k] * cd_x[k] + g_xxz_x[k];
        }

        /// Set up 3-4 components of targeted buffer : cbuffer.data(

        auto g_xxyy_0 = cbuffer.data(gs_off + 3);

        #pragma omp simd aligned(cd_x, g_xxyy_0, g_xyy_0, g_xyy_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxyy_0[k] = -g_xyy_0[k] * cd_x[k] + g_xyy_x[k];
        }

        /// Set up 4-5 components of targeted buffer : cbuffer.data(

        auto g_xxyz_0 = cbuffer.data(gs_off + 4);

        #pragma omp simd aligned(cd_x, g_xxyz_0, g_xyz_0, g_xyz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxyz_0[k] = -g_xyz_0[k] * cd_x[k] + g_xyz_x[k];
        }

        /// Set up 5-6 components of targeted buffer : cbuffer.data(

        auto g_xxzz_0 = cbuffer.data(gs_off + 5);

        #pragma omp simd aligned(cd_x, g_xxzz_0, g_xzz_0, g_xzz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxzz_0[k] = -g_xzz_0[k] * cd_x[k] + g_xzz_x[k];
        }

        /// Set up 6-7 components of targeted buffer : cbuffer.data(

        auto g_xyyy_0 = cbuffer.data(gs_off + 6);

        #pragma omp simd aligned(cd_x, g_xyyy_0, g_yyy_0, g_yyy_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xyyy_0[k] = -g_yyy_0[k] * cd_x[k] + g_yyy_x[k];
        }

        /// Set up 7-8 components of targeted buffer : cbuffer.data(

        auto g_xyyz_0 = cbuffer.data(gs_off + 7);

        #pragma omp simd aligned(cd_x, g_xyyz_0, g_yyz_0, g_yyz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xyyz_0[k] = -g_yyz_0[k] * cd_x[k] + g_yyz_x[k];
        }

        /// Set up 8-9 components of targeted buffer : cbuffer.data(

        auto g_xyzz_0 = cbuffer.data(gs_off + 8);

        #pragma omp simd aligned(cd_x, g_xyzz_0, g_yzz_0, g_yzz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xyzz_0[k] = -g_yzz_0[k] * cd_x[k] + g_yzz_x[k];
        }

        /// Set up 9-10 components of targeted buffer : cbuffer.data(

        auto g_xzzz_0 = cbuffer.data(gs_off + 9);

        #pragma omp simd aligned(cd_x, g_xzzz_0, g_zzz_0, g_zzz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xzzz_0[k] = -g_zzz_0[k] * cd_x[k] + g_zzz_x[k];
        }

        /// Set up 10-11 components of targeted buffer : cbuffer.data(

        auto g_yyyy_0 = cbuffer.data(gs_off + 10);

        #pragma omp simd aligned(cd_y, g_yyy_0, g_yyy_y, g_yyyy_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yyyy_0[k] = -g_yyy_0[k] * cd_y[k] + g_yyy_y[k];
        }

        /// Set up 11-12 components of targeted buffer : cbuffer.data(

        auto g_yyyz_0 = cbuffer.data(gs_off + 11);

        #pragma omp simd aligned(cd_y, g_yyyz_0, g_yyz_0, g_yyz_y  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yyyz_0[k] = -g_yyz_0[k] * cd_y[k] + g_yyz_y[k];
        }

        /// Set up 12-13 components of targeted buffer : cbuffer.data(

        auto g_yyzz_0 = cbuffer.data(gs_off + 12);

        #pragma omp simd aligned(cd_y, g_yyzz_0, g_yzz_0, g_yzz_y  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yyzz_0[k] = -g_yzz_0[k] * cd_y[k] + g_yzz_y[k];
        }

        /// Set up 13-14 components of targeted buffer : cbuffer.data(

        auto g_yzzz_0 = cbuffer.data(gs_off + 13);

        #pragma omp simd aligned(cd_y, g_yzzz_0, g_zzz_0, g_zzz_y  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yzzz_0[k] = -g_zzz_0[k] * cd_y[k] + g_zzz_y[k];
        }

        /// Set up 14-15 components of targeted buffer : cbuffer.data(

        auto g_zzzz_0 = cbuffer.data(gs_off + 14);

        #pragma omp simd aligned(cd_z, g_zzz_0, g_zzz_z, g_zzzz_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_zzzz_0[k] = -g_zzz_0[k] * cd_z[k] + g_zzz_z[k];
        }
    }
}

} // t3ceri namespace

