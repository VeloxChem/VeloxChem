#include "ElectronRepulsionGeom0010ContrRecXXHS.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxhs(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxhs,
                                            const size_t idx_xxgs,
                                            const size_t idx_geom_10_xxgs,
                                            const size_t idx_geom_10_xxgp,
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
            /// Set up components of auxilary buffer : SSGS

            const auto gs_off = idx_xxgs + (i * bcomps + j) * 15;

            auto g_xxxx_0 = cbuffer.data(gs_off + 0);

            auto g_yyyy_0 = cbuffer.data(gs_off + 10);

            auto g_zzzz_0 = cbuffer.data(gs_off + 14);

            /// Set up components of auxilary buffer : SSGS

            const auto gs_geom_10_off = idx_geom_10_xxgs + (i * bcomps + j) * 15;

            auto g_x_0_xxxx_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxy_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxyy_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxyz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxzz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xyyy_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xyyz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xyzz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xzzz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_yyyy_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_yyyz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_yyzz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_yzzz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_zzzz_0 = cbuffer.data(gs_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_y_0_xxxx_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 0);

            auto g_y_0_xxxy_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 1);

            auto g_y_0_xxxz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 2);

            auto g_y_0_xxyy_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 3);

            auto g_y_0_xxyz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 4);

            auto g_y_0_xxzz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 5);

            auto g_y_0_xyyy_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 6);

            auto g_y_0_xyyz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 7);

            auto g_y_0_xyzz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 8);

            auto g_y_0_xzzz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 9);

            auto g_y_0_yyyy_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 10);

            auto g_y_0_yyyz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 11);

            auto g_y_0_yyzz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 12);

            auto g_y_0_yzzz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 13);

            auto g_y_0_zzzz_0 = cbuffer.data(gs_geom_10_off + 15 * acomps * bcomps + 14);

            auto g_z_0_xxxx_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 0);

            auto g_z_0_xxxy_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 1);

            auto g_z_0_xxxz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 2);

            auto g_z_0_xxyy_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 3);

            auto g_z_0_xxyz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 4);

            auto g_z_0_xxzz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 5);

            auto g_z_0_xyyy_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 6);

            auto g_z_0_xyyz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 7);

            auto g_z_0_xyzz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 8);

            auto g_z_0_xzzz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 9);

            auto g_z_0_yyyy_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 10);

            auto g_z_0_yyyz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 11);

            auto g_z_0_yyzz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 12);

            auto g_z_0_yzzz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 13);

            auto g_z_0_zzzz_0 = cbuffer.data(gs_geom_10_off + 30 * acomps * bcomps + 14);

            /// Set up components of auxilary buffer : SSGP

            const auto gp_geom_10_off = idx_geom_10_xxgp + (i * bcomps + j) * 45;

            auto g_x_0_xxxx_x = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxx_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxx_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxy_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxyy_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxyz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxzz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxzz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xyyy_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xyyz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xyzz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xzzz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xzzz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_yyyy_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_yyyz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_yyzz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_yzzz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_zzzz_y = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_zzzz_z = cbuffer.data(gp_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_y_0_xxxx_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 0);

            auto g_y_0_xxxy_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 3);

            auto g_y_0_xxxz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 6);

            auto g_y_0_xxyy_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 9);

            auto g_y_0_xxyz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 12);

            auto g_y_0_xxzz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 15);

            auto g_y_0_xyyy_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 18);

            auto g_y_0_xyyz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 21);

            auto g_y_0_xyzz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 24);

            auto g_y_0_xzzz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 27);

            auto g_y_0_yyyy_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 30);

            auto g_y_0_yyyy_y = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 31);

            auto g_y_0_yyyy_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 32);

            auto g_y_0_yyyz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 33);

            auto g_y_0_yyyz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 35);

            auto g_y_0_yyzz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 36);

            auto g_y_0_yyzz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 38);

            auto g_y_0_yzzz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 39);

            auto g_y_0_yzzz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 41);

            auto g_y_0_zzzz_x = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 42);

            auto g_y_0_zzzz_z = cbuffer.data(gp_geom_10_off + 45 * acomps * bcomps + 44);

            auto g_z_0_xxxx_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 0);

            auto g_z_0_xxxy_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 3);

            auto g_z_0_xxxz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 6);

            auto g_z_0_xxyy_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 9);

            auto g_z_0_xxyz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 12);

            auto g_z_0_xxzz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 15);

            auto g_z_0_xyyy_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 18);

            auto g_z_0_xyyz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 21);

            auto g_z_0_xyzz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 24);

            auto g_z_0_xzzz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 27);

            auto g_z_0_yyyy_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 30);

            auto g_z_0_yyyy_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 31);

            auto g_z_0_yyyz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 33);

            auto g_z_0_yyyz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 34);

            auto g_z_0_yyzz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 36);

            auto g_z_0_yyzz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 37);

            auto g_z_0_yzzz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 39);

            auto g_z_0_yzzz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 40);

            auto g_z_0_zzzz_x = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 42);

            auto g_z_0_zzzz_y = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 43);

            auto g_z_0_zzzz_z = cbuffer.data(gp_geom_10_off + 90 * acomps * bcomps + 44);

            /// set up bra offset for contr_buffer_xxhs

            const auto hs_geom_10_off = idx_geom_10_xxhs + (i * bcomps + j) * 21;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxx_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 0);

            #pragma omp simd aligned(cd_x, g_x_0_xxxx_0, g_x_0_xxxx_x, g_x_0_xxxxx_0, g_xxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxx_0[k] = -g_xxxx_0[k] - g_x_0_xxxx_0[k] * cd_x[k] + g_x_0_xxxx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxy_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 1);

            #pragma omp simd aligned(cd_y, g_x_0_xxxx_0, g_x_0_xxxx_y, g_x_0_xxxxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxy_0[k] = -g_x_0_xxxx_0[k] * cd_y[k] + g_x_0_xxxx_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_z, g_x_0_xxxx_0, g_x_0_xxxx_z, g_x_0_xxxxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxz_0[k] = -g_x_0_xxxx_0[k] * cd_z[k] + g_x_0_xxxx_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyy_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 3);

            #pragma omp simd aligned(cd_y, g_x_0_xxxy_0, g_x_0_xxxy_y, g_x_0_xxxyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyy_0[k] = -g_x_0_xxxy_0[k] * cd_y[k] + g_x_0_xxxy_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 4);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyz_0, g_x_0_xxxz_0, g_x_0_xxxz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyz_0[k] = -g_x_0_xxxz_0[k] * cd_y[k] + g_x_0_xxxz_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxzz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_z, g_x_0_xxxz_0, g_x_0_xxxz_z, g_x_0_xxxzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzz_0[k] = -g_x_0_xxxz_0[k] * cd_z[k] + g_x_0_xxxz_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyy_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 6);

            #pragma omp simd aligned(cd_y, g_x_0_xxyy_0, g_x_0_xxyy_y, g_x_0_xxyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyy_0[k] = -g_x_0_xxyy_0[k] * cd_y[k] + g_x_0_xxyy_y[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 7);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyz_0, g_x_0_xxyz_0, g_x_0_xxyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyz_0[k] = -g_x_0_xxyz_0[k] * cd_y[k] + g_x_0_xxyz_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyzz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_y, g_x_0_xxyzz_0, g_x_0_xxzz_0, g_x_0_xxzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzz_0[k] = -g_x_0_xxzz_0[k] * cd_y[k] + g_x_0_xxzz_y[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzzz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_z, g_x_0_xxzz_0, g_x_0_xxzz_z, g_x_0_xxzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzz_0[k] = -g_x_0_xxzz_0[k] * cd_z[k] + g_x_0_xxzz_z[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyy_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 10);

            #pragma omp simd aligned(cd_y, g_x_0_xyyy_0, g_x_0_xyyy_y, g_x_0_xyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyy_0[k] = -g_x_0_xyyy_0[k] * cd_y[k] + g_x_0_xyyy_y[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyz_0, g_x_0_xyyz_0, g_x_0_xyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyz_0[k] = -g_x_0_xyyz_0[k] * cd_y[k] + g_x_0_xyyz_y[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyzz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 12);

            #pragma omp simd aligned(cd_y, g_x_0_xyyzz_0, g_x_0_xyzz_0, g_x_0_xyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzz_0[k] = -g_x_0_xyzz_0[k] * cd_y[k] + g_x_0_xyzz_y[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzzz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 13);

            #pragma omp simd aligned(cd_y, g_x_0_xyzzz_0, g_x_0_xzzz_0, g_x_0_xzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzz_0[k] = -g_x_0_xzzz_0[k] * cd_y[k] + g_x_0_xzzz_y[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzzz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_z, g_x_0_xzzz_0, g_x_0_xzzz_z, g_x_0_xzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzz_0[k] = -g_x_0_xzzz_0[k] * cd_z[k] + g_x_0_xzzz_z[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyy_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 15);

            #pragma omp simd aligned(cd_y, g_x_0_yyyy_0, g_x_0_yyyy_y, g_x_0_yyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyy_0[k] = -g_x_0_yyyy_0[k] * cd_y[k] + g_x_0_yyyy_y[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 16);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyz_0, g_x_0_yyyz_0, g_x_0_yyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyz_0[k] = -g_x_0_yyyz_0[k] * cd_y[k] + g_x_0_yyyz_y[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyzz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_y, g_x_0_yyyzz_0, g_x_0_yyzz_0, g_x_0_yyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzz_0[k] = -g_x_0_yyzz_0[k] * cd_y[k] + g_x_0_yyzz_y[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzzz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 18);

            #pragma omp simd aligned(cd_y, g_x_0_yyzzz_0, g_x_0_yzzz_0, g_x_0_yzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzz_0[k] = -g_x_0_yzzz_0[k] * cd_y[k] + g_x_0_yzzz_y[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzzz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 19);

            #pragma omp simd aligned(cd_y, g_x_0_yzzzz_0, g_x_0_zzzz_0, g_x_0_zzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzz_0[k] = -g_x_0_zzzz_0[k] * cd_y[k] + g_x_0_zzzz_y[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzzz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_z, g_x_0_zzzz_0, g_x_0_zzzz_z, g_x_0_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzz_0[k] = -g_x_0_zzzz_0[k] * cd_z[k] + g_x_0_zzzz_z[k];
            }
            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxx_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 0);

            #pragma omp simd aligned(cd_x, g_y_0_xxxx_0, g_y_0_xxxx_x, g_y_0_xxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxx_0[k] = -g_y_0_xxxx_0[k] * cd_x[k] + g_y_0_xxxx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxy_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 1);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxy_0, g_y_0_xxxy_0, g_y_0_xxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxy_0[k] = -g_y_0_xxxy_0[k] * cd_x[k] + g_y_0_xxxy_x[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxz_0, g_y_0_xxxz_0, g_y_0_xxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxz_0[k] = -g_y_0_xxxz_0[k] * cd_x[k] + g_y_0_xxxz_x[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyy_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 3);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyy_0, g_y_0_xxyy_0, g_y_0_xxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyy_0[k] = -g_y_0_xxyy_0[k] * cd_x[k] + g_y_0_xxyy_x[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 4);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyz_0, g_y_0_xxyz_0, g_y_0_xxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyz_0[k] = -g_y_0_xxyz_0[k] * cd_x[k] + g_y_0_xxyz_x[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxzz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_y_0_xxxzz_0, g_y_0_xxzz_0, g_y_0_xxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzz_0[k] = -g_y_0_xxzz_0[k] * cd_x[k] + g_y_0_xxzz_x[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyy_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 6);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyy_0, g_y_0_xyyy_0, g_y_0_xyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyy_0[k] = -g_y_0_xyyy_0[k] * cd_x[k] + g_y_0_xyyy_x[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 7);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyz_0, g_y_0_xyyz_0, g_y_0_xyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyz_0[k] = -g_y_0_xyyz_0[k] * cd_x[k] + g_y_0_xyyz_x[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyzz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_x, g_y_0_xxyzz_0, g_y_0_xyzz_0, g_y_0_xyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzz_0[k] = -g_y_0_xyzz_0[k] * cd_x[k] + g_y_0_xyzz_x[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzzz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_x, g_y_0_xxzzz_0, g_y_0_xzzz_0, g_y_0_xzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzz_0[k] = -g_y_0_xzzz_0[k] * cd_x[k] + g_y_0_xzzz_x[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyy_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 10);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyy_0, g_y_0_yyyy_0, g_y_0_yyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyy_0[k] = -g_y_0_yyyy_0[k] * cd_x[k] + g_y_0_yyyy_x[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyz_0, g_y_0_yyyz_0, g_y_0_yyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyz_0[k] = -g_y_0_yyyz_0[k] * cd_x[k] + g_y_0_yyyz_x[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyzz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 12);

            #pragma omp simd aligned(cd_x, g_y_0_xyyzz_0, g_y_0_yyzz_0, g_y_0_yyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzz_0[k] = -g_y_0_yyzz_0[k] * cd_x[k] + g_y_0_yyzz_x[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzzz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 13);

            #pragma omp simd aligned(cd_x, g_y_0_xyzzz_0, g_y_0_yzzz_0, g_y_0_yzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzz_0[k] = -g_y_0_yzzz_0[k] * cd_x[k] + g_y_0_yzzz_x[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzzz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_x, g_y_0_xzzzz_0, g_y_0_zzzz_0, g_y_0_zzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzz_0[k] = -g_y_0_zzzz_0[k] * cd_x[k] + g_y_0_zzzz_x[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyy_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 15);

            #pragma omp simd aligned(cd_y, g_y_0_yyyy_0, g_y_0_yyyy_y, g_y_0_yyyyy_0, g_yyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyy_0[k] = -g_yyyy_0[k] - g_y_0_yyyy_0[k] * cd_y[k] + g_y_0_yyyy_y[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 16);

            #pragma omp simd aligned(cd_z, g_y_0_yyyy_0, g_y_0_yyyy_z, g_y_0_yyyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyz_0[k] = -g_y_0_yyyy_0[k] * cd_z[k] + g_y_0_yyyy_z[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyzz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_z, g_y_0_yyyz_0, g_y_0_yyyz_z, g_y_0_yyyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzz_0[k] = -g_y_0_yyyz_0[k] * cd_z[k] + g_y_0_yyyz_z[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzzz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 18);

            #pragma omp simd aligned(cd_z, g_y_0_yyzz_0, g_y_0_yyzz_z, g_y_0_yyzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzz_0[k] = -g_y_0_yyzz_0[k] * cd_z[k] + g_y_0_yyzz_z[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzzz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 19);

            #pragma omp simd aligned(cd_z, g_y_0_yzzz_0, g_y_0_yzzz_z, g_y_0_yzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzz_0[k] = -g_y_0_yzzz_0[k] * cd_z[k] + g_y_0_yzzz_z[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzzz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_z, g_y_0_zzzz_0, g_y_0_zzzz_z, g_y_0_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzz_0[k] = -g_y_0_zzzz_0[k] * cd_z[k] + g_y_0_zzzz_z[k];
            }
            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxx_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 0);

            #pragma omp simd aligned(cd_x, g_z_0_xxxx_0, g_z_0_xxxx_x, g_z_0_xxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxx_0[k] = -g_z_0_xxxx_0[k] * cd_x[k] + g_z_0_xxxx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxy_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 1);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxy_0, g_z_0_xxxy_0, g_z_0_xxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxy_0[k] = -g_z_0_xxxy_0[k] * cd_x[k] + g_z_0_xxxy_x[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxz_0, g_z_0_xxxz_0, g_z_0_xxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxz_0[k] = -g_z_0_xxxz_0[k] * cd_x[k] + g_z_0_xxxz_x[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyy_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 3);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyy_0, g_z_0_xxyy_0, g_z_0_xxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyy_0[k] = -g_z_0_xxyy_0[k] * cd_x[k] + g_z_0_xxyy_x[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 4);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyz_0, g_z_0_xxyz_0, g_z_0_xxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyz_0[k] = -g_z_0_xxyz_0[k] * cd_x[k] + g_z_0_xxyz_x[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxzz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_z_0_xxxzz_0, g_z_0_xxzz_0, g_z_0_xxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzz_0[k] = -g_z_0_xxzz_0[k] * cd_x[k] + g_z_0_xxzz_x[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyy_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 6);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyy_0, g_z_0_xyyy_0, g_z_0_xyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyy_0[k] = -g_z_0_xyyy_0[k] * cd_x[k] + g_z_0_xyyy_x[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 7);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyz_0, g_z_0_xyyz_0, g_z_0_xyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyz_0[k] = -g_z_0_xyyz_0[k] * cd_x[k] + g_z_0_xyyz_x[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyzz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_x, g_z_0_xxyzz_0, g_z_0_xyzz_0, g_z_0_xyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzz_0[k] = -g_z_0_xyzz_0[k] * cd_x[k] + g_z_0_xyzz_x[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzzz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_x, g_z_0_xxzzz_0, g_z_0_xzzz_0, g_z_0_xzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzz_0[k] = -g_z_0_xzzz_0[k] * cd_x[k] + g_z_0_xzzz_x[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyy_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 10);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyy_0, g_z_0_yyyy_0, g_z_0_yyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyy_0[k] = -g_z_0_yyyy_0[k] * cd_x[k] + g_z_0_yyyy_x[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyz_0, g_z_0_yyyz_0, g_z_0_yyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyz_0[k] = -g_z_0_yyyz_0[k] * cd_x[k] + g_z_0_yyyz_x[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyzz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 12);

            #pragma omp simd aligned(cd_x, g_z_0_xyyzz_0, g_z_0_yyzz_0, g_z_0_yyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzz_0[k] = -g_z_0_yyzz_0[k] * cd_x[k] + g_z_0_yyzz_x[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzzz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 13);

            #pragma omp simd aligned(cd_x, g_z_0_xyzzz_0, g_z_0_yzzz_0, g_z_0_yzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzz_0[k] = -g_z_0_yzzz_0[k] * cd_x[k] + g_z_0_yzzz_x[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzzz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_x, g_z_0_xzzzz_0, g_z_0_zzzz_0, g_z_0_zzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzz_0[k] = -g_z_0_zzzz_0[k] * cd_x[k] + g_z_0_zzzz_x[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyy_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 15);

            #pragma omp simd aligned(cd_y, g_z_0_yyyy_0, g_z_0_yyyy_y, g_z_0_yyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyy_0[k] = -g_z_0_yyyy_0[k] * cd_y[k] + g_z_0_yyyy_y[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 16);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyz_0, g_z_0_yyyz_0, g_z_0_yyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyz_0[k] = -g_z_0_yyyz_0[k] * cd_y[k] + g_z_0_yyyz_y[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyzz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_y, g_z_0_yyyzz_0, g_z_0_yyzz_0, g_z_0_yyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzz_0[k] = -g_z_0_yyzz_0[k] * cd_y[k] + g_z_0_yyzz_y[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzzz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 18);

            #pragma omp simd aligned(cd_y, g_z_0_yyzzz_0, g_z_0_yzzz_0, g_z_0_yzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzz_0[k] = -g_z_0_yzzz_0[k] * cd_y[k] + g_z_0_yzzz_y[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzzz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 19);

            #pragma omp simd aligned(cd_y, g_z_0_yzzzz_0, g_z_0_zzzz_0, g_z_0_zzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzz_0[k] = -g_z_0_zzzz_0[k] * cd_y[k] + g_z_0_zzzz_y[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzzz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_z, g_z_0_zzzz_0, g_z_0_zzzz_z, g_z_0_zzzzz_0, g_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzz_0[k] = -g_zzzz_0[k] - g_z_0_zzzz_0[k] * cd_z[k] + g_z_0_zzzz_z[k];
            }
        }
    }
}

} // erirec namespace

