#include "ElectronRepulsionGeom0010ContrRecXXPF.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxpf(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxpf,
                                            CSimdArray<double>& pbuffer,
                                            const size_t idx_xxsf,
                                            const size_t idx_geom_10_xxsf,
                                            const size_t idx_geom_10_xxsg,
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
            /// Set up components of auxilary buffer : SSSF

            const auto sf_off = idx_xxsf + (i * bcomps + j) * 10;

            auto g_0_xxx = pbuffer.data(sf_off + 0);

            auto g_0_xxy = pbuffer.data(sf_off + 1);

            auto g_0_xxz = pbuffer.data(sf_off + 2);

            auto g_0_xyy = pbuffer.data(sf_off + 3);

            auto g_0_xyz = pbuffer.data(sf_off + 4);

            auto g_0_xzz = pbuffer.data(sf_off + 5);

            auto g_0_yyy = pbuffer.data(sf_off + 6);

            auto g_0_yyz = pbuffer.data(sf_off + 7);

            auto g_0_yzz = pbuffer.data(sf_off + 8);

            auto g_0_zzz = pbuffer.data(sf_off + 9);

            /// Set up components of auxilary buffer : SSSF

            const auto sf_geom_10_off = idx_geom_10_xxsf + (i * bcomps + j) * 10;

            auto g_x_0_0_xxx = cbuffer.data(sf_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_0_xxy = cbuffer.data(sf_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_0_xxz = cbuffer.data(sf_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_0_xyy = cbuffer.data(sf_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_0_xyz = cbuffer.data(sf_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_0_xzz = cbuffer.data(sf_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_0_yyy = cbuffer.data(sf_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_0_yyz = cbuffer.data(sf_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_0_yzz = cbuffer.data(sf_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_0_zzz = cbuffer.data(sf_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_y_0_0_xxx = cbuffer.data(sf_geom_10_off + 10 * acomps * bcomps + 0);

            auto g_y_0_0_xxy = cbuffer.data(sf_geom_10_off + 10 * acomps * bcomps + 1);

            auto g_y_0_0_xxz = cbuffer.data(sf_geom_10_off + 10 * acomps * bcomps + 2);

            auto g_y_0_0_xyy = cbuffer.data(sf_geom_10_off + 10 * acomps * bcomps + 3);

            auto g_y_0_0_xyz = cbuffer.data(sf_geom_10_off + 10 * acomps * bcomps + 4);

            auto g_y_0_0_xzz = cbuffer.data(sf_geom_10_off + 10 * acomps * bcomps + 5);

            auto g_y_0_0_yyy = cbuffer.data(sf_geom_10_off + 10 * acomps * bcomps + 6);

            auto g_y_0_0_yyz = cbuffer.data(sf_geom_10_off + 10 * acomps * bcomps + 7);

            auto g_y_0_0_yzz = cbuffer.data(sf_geom_10_off + 10 * acomps * bcomps + 8);

            auto g_y_0_0_zzz = cbuffer.data(sf_geom_10_off + 10 * acomps * bcomps + 9);

            auto g_z_0_0_xxx = cbuffer.data(sf_geom_10_off + 20 * acomps * bcomps + 0);

            auto g_z_0_0_xxy = cbuffer.data(sf_geom_10_off + 20 * acomps * bcomps + 1);

            auto g_z_0_0_xxz = cbuffer.data(sf_geom_10_off + 20 * acomps * bcomps + 2);

            auto g_z_0_0_xyy = cbuffer.data(sf_geom_10_off + 20 * acomps * bcomps + 3);

            auto g_z_0_0_xyz = cbuffer.data(sf_geom_10_off + 20 * acomps * bcomps + 4);

            auto g_z_0_0_xzz = cbuffer.data(sf_geom_10_off + 20 * acomps * bcomps + 5);

            auto g_z_0_0_yyy = cbuffer.data(sf_geom_10_off + 20 * acomps * bcomps + 6);

            auto g_z_0_0_yyz = cbuffer.data(sf_geom_10_off + 20 * acomps * bcomps + 7);

            auto g_z_0_0_yzz = cbuffer.data(sf_geom_10_off + 20 * acomps * bcomps + 8);

            auto g_z_0_0_zzz = cbuffer.data(sf_geom_10_off + 20 * acomps * bcomps + 9);

            /// Set up components of auxilary buffer : SSSG

            const auto sg_geom_10_off = idx_geom_10_xxsg + (i * bcomps + j) * 15;

            auto g_x_0_0_xxxx = cbuffer.data(sg_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_0_xxxy = cbuffer.data(sg_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_0_xxxz = cbuffer.data(sg_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_0_xxyy = cbuffer.data(sg_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_0_xxyz = cbuffer.data(sg_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_0_xxzz = cbuffer.data(sg_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_0_xyyy = cbuffer.data(sg_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_0_xyyz = cbuffer.data(sg_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_0_xyzz = cbuffer.data(sg_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_0_xzzz = cbuffer.data(sg_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_0_yyyy = cbuffer.data(sg_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_0_yyyz = cbuffer.data(sg_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_0_yyzz = cbuffer.data(sg_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_0_yzzz = cbuffer.data(sg_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_0_zzzz = cbuffer.data(sg_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_y_0_0_xxxx = cbuffer.data(sg_geom_10_off + 15 * acomps * bcomps + 0);

            auto g_y_0_0_xxxy = cbuffer.data(sg_geom_10_off + 15 * acomps * bcomps + 1);

            auto g_y_0_0_xxxz = cbuffer.data(sg_geom_10_off + 15 * acomps * bcomps + 2);

            auto g_y_0_0_xxyy = cbuffer.data(sg_geom_10_off + 15 * acomps * bcomps + 3);

            auto g_y_0_0_xxyz = cbuffer.data(sg_geom_10_off + 15 * acomps * bcomps + 4);

            auto g_y_0_0_xxzz = cbuffer.data(sg_geom_10_off + 15 * acomps * bcomps + 5);

            auto g_y_0_0_xyyy = cbuffer.data(sg_geom_10_off + 15 * acomps * bcomps + 6);

            auto g_y_0_0_xyyz = cbuffer.data(sg_geom_10_off + 15 * acomps * bcomps + 7);

            auto g_y_0_0_xyzz = cbuffer.data(sg_geom_10_off + 15 * acomps * bcomps + 8);

            auto g_y_0_0_xzzz = cbuffer.data(sg_geom_10_off + 15 * acomps * bcomps + 9);

            auto g_y_0_0_yyyy = cbuffer.data(sg_geom_10_off + 15 * acomps * bcomps + 10);

            auto g_y_0_0_yyyz = cbuffer.data(sg_geom_10_off + 15 * acomps * bcomps + 11);

            auto g_y_0_0_yyzz = cbuffer.data(sg_geom_10_off + 15 * acomps * bcomps + 12);

            auto g_y_0_0_yzzz = cbuffer.data(sg_geom_10_off + 15 * acomps * bcomps + 13);

            auto g_y_0_0_zzzz = cbuffer.data(sg_geom_10_off + 15 * acomps * bcomps + 14);

            auto g_z_0_0_xxxx = cbuffer.data(sg_geom_10_off + 30 * acomps * bcomps + 0);

            auto g_z_0_0_xxxy = cbuffer.data(sg_geom_10_off + 30 * acomps * bcomps + 1);

            auto g_z_0_0_xxxz = cbuffer.data(sg_geom_10_off + 30 * acomps * bcomps + 2);

            auto g_z_0_0_xxyy = cbuffer.data(sg_geom_10_off + 30 * acomps * bcomps + 3);

            auto g_z_0_0_xxyz = cbuffer.data(sg_geom_10_off + 30 * acomps * bcomps + 4);

            auto g_z_0_0_xxzz = cbuffer.data(sg_geom_10_off + 30 * acomps * bcomps + 5);

            auto g_z_0_0_xyyy = cbuffer.data(sg_geom_10_off + 30 * acomps * bcomps + 6);

            auto g_z_0_0_xyyz = cbuffer.data(sg_geom_10_off + 30 * acomps * bcomps + 7);

            auto g_z_0_0_xyzz = cbuffer.data(sg_geom_10_off + 30 * acomps * bcomps + 8);

            auto g_z_0_0_xzzz = cbuffer.data(sg_geom_10_off + 30 * acomps * bcomps + 9);

            auto g_z_0_0_yyyy = cbuffer.data(sg_geom_10_off + 30 * acomps * bcomps + 10);

            auto g_z_0_0_yyyz = cbuffer.data(sg_geom_10_off + 30 * acomps * bcomps + 11);

            auto g_z_0_0_yyzz = cbuffer.data(sg_geom_10_off + 30 * acomps * bcomps + 12);

            auto g_z_0_0_yzzz = cbuffer.data(sg_geom_10_off + 30 * acomps * bcomps + 13);

            auto g_z_0_0_zzzz = cbuffer.data(sg_geom_10_off + 30 * acomps * bcomps + 14);

            /// set up bra offset for contr_buffer_xxpf

            const auto pf_geom_10_off = idx_geom_10_xxpf + (i * bcomps + j) * 30;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_xxx = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_x_xxy = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_x_xxz = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_x_xyy = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_x_xyz = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_x_xzz = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_x_yyy = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_x_yyz = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_x_yzz = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_x_zzz = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_x, g_0_xxx, g_0_xxy, g_0_xxz, g_0_xyy, g_0_xyz, g_0_xzz, g_0_yyy, g_0_yyz, g_0_yzz, g_0_zzz, g_x_0_0_xxx, g_x_0_0_xxxx, g_x_0_0_xxxy, g_x_0_0_xxxz, g_x_0_0_xxy, g_x_0_0_xxyy, g_x_0_0_xxyz, g_x_0_0_xxz, g_x_0_0_xxzz, g_x_0_0_xyy, g_x_0_0_xyyy, g_x_0_0_xyyz, g_x_0_0_xyz, g_x_0_0_xyzz, g_x_0_0_xzz, g_x_0_0_xzzz, g_x_0_0_yyy, g_x_0_0_yyz, g_x_0_0_yzz, g_x_0_0_zzz, g_x_0_x_xxx, g_x_0_x_xxy, g_x_0_x_xxz, g_x_0_x_xyy, g_x_0_x_xyz, g_x_0_x_xzz, g_x_0_x_yyy, g_x_0_x_yyz, g_x_0_x_yzz, g_x_0_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_xxx[k] = -g_0_xxx[k] - g_x_0_0_xxx[k] * cd_x[k] + g_x_0_0_xxxx[k];

                g_x_0_x_xxy[k] = -g_0_xxy[k] - g_x_0_0_xxy[k] * cd_x[k] + g_x_0_0_xxxy[k];

                g_x_0_x_xxz[k] = -g_0_xxz[k] - g_x_0_0_xxz[k] * cd_x[k] + g_x_0_0_xxxz[k];

                g_x_0_x_xyy[k] = -g_0_xyy[k] - g_x_0_0_xyy[k] * cd_x[k] + g_x_0_0_xxyy[k];

                g_x_0_x_xyz[k] = -g_0_xyz[k] - g_x_0_0_xyz[k] * cd_x[k] + g_x_0_0_xxyz[k];

                g_x_0_x_xzz[k] = -g_0_xzz[k] - g_x_0_0_xzz[k] * cd_x[k] + g_x_0_0_xxzz[k];

                g_x_0_x_yyy[k] = -g_0_yyy[k] - g_x_0_0_yyy[k] * cd_x[k] + g_x_0_0_xyyy[k];

                g_x_0_x_yyz[k] = -g_0_yyz[k] - g_x_0_0_yyz[k] * cd_x[k] + g_x_0_0_xyyz[k];

                g_x_0_x_yzz[k] = -g_0_yzz[k] - g_x_0_0_yzz[k] * cd_x[k] + g_x_0_0_xyzz[k];

                g_x_0_x_zzz[k] = -g_0_zzz[k] - g_x_0_0_zzz[k] * cd_x[k] + g_x_0_0_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_xxx = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_y_xxy = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_y_xxz = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_y_xyy = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_y_xyz = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_y_xzz = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_y_yyy = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_y_yyz = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_y_yzz = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_y_zzz = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 19);

            #pragma omp simd aligned(cd_y, g_x_0_0_xxx, g_x_0_0_xxxy, g_x_0_0_xxy, g_x_0_0_xxyy, g_x_0_0_xxyz, g_x_0_0_xxz, g_x_0_0_xyy, g_x_0_0_xyyy, g_x_0_0_xyyz, g_x_0_0_xyz, g_x_0_0_xyzz, g_x_0_0_xzz, g_x_0_0_yyy, g_x_0_0_yyyy, g_x_0_0_yyyz, g_x_0_0_yyz, g_x_0_0_yyzz, g_x_0_0_yzz, g_x_0_0_yzzz, g_x_0_0_zzz, g_x_0_y_xxx, g_x_0_y_xxy, g_x_0_y_xxz, g_x_0_y_xyy, g_x_0_y_xyz, g_x_0_y_xzz, g_x_0_y_yyy, g_x_0_y_yyz, g_x_0_y_yzz, g_x_0_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_xxx[k] = -g_x_0_0_xxx[k] * cd_y[k] + g_x_0_0_xxxy[k];

                g_x_0_y_xxy[k] = -g_x_0_0_xxy[k] * cd_y[k] + g_x_0_0_xxyy[k];

                g_x_0_y_xxz[k] = -g_x_0_0_xxz[k] * cd_y[k] + g_x_0_0_xxyz[k];

                g_x_0_y_xyy[k] = -g_x_0_0_xyy[k] * cd_y[k] + g_x_0_0_xyyy[k];

                g_x_0_y_xyz[k] = -g_x_0_0_xyz[k] * cd_y[k] + g_x_0_0_xyyz[k];

                g_x_0_y_xzz[k] = -g_x_0_0_xzz[k] * cd_y[k] + g_x_0_0_xyzz[k];

                g_x_0_y_yyy[k] = -g_x_0_0_yyy[k] * cd_y[k] + g_x_0_0_yyyy[k];

                g_x_0_y_yyz[k] = -g_x_0_0_yyz[k] * cd_y[k] + g_x_0_0_yyyz[k];

                g_x_0_y_yzz[k] = -g_x_0_0_yzz[k] * cd_y[k] + g_x_0_0_yyzz[k];

                g_x_0_y_zzz[k] = -g_x_0_0_zzz[k] * cd_y[k] + g_x_0_0_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_xxx = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_z_xxy = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_z_xxz = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_z_xyy = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_z_xyz = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_z_xzz = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_z_yyy = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_z_yyz = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_z_yzz = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_z_zzz = cbuffer.data(pf_geom_10_off + 0 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_z, g_x_0_0_xxx, g_x_0_0_xxxz, g_x_0_0_xxy, g_x_0_0_xxyz, g_x_0_0_xxz, g_x_0_0_xxzz, g_x_0_0_xyy, g_x_0_0_xyyz, g_x_0_0_xyz, g_x_0_0_xyzz, g_x_0_0_xzz, g_x_0_0_xzzz, g_x_0_0_yyy, g_x_0_0_yyyz, g_x_0_0_yyz, g_x_0_0_yyzz, g_x_0_0_yzz, g_x_0_0_yzzz, g_x_0_0_zzz, g_x_0_0_zzzz, g_x_0_z_xxx, g_x_0_z_xxy, g_x_0_z_xxz, g_x_0_z_xyy, g_x_0_z_xyz, g_x_0_z_xzz, g_x_0_z_yyy, g_x_0_z_yyz, g_x_0_z_yzz, g_x_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_xxx[k] = -g_x_0_0_xxx[k] * cd_z[k] + g_x_0_0_xxxz[k];

                g_x_0_z_xxy[k] = -g_x_0_0_xxy[k] * cd_z[k] + g_x_0_0_xxyz[k];

                g_x_0_z_xxz[k] = -g_x_0_0_xxz[k] * cd_z[k] + g_x_0_0_xxzz[k];

                g_x_0_z_xyy[k] = -g_x_0_0_xyy[k] * cd_z[k] + g_x_0_0_xyyz[k];

                g_x_0_z_xyz[k] = -g_x_0_0_xyz[k] * cd_z[k] + g_x_0_0_xyzz[k];

                g_x_0_z_xzz[k] = -g_x_0_0_xzz[k] * cd_z[k] + g_x_0_0_xzzz[k];

                g_x_0_z_yyy[k] = -g_x_0_0_yyy[k] * cd_z[k] + g_x_0_0_yyyz[k];

                g_x_0_z_yyz[k] = -g_x_0_0_yyz[k] * cd_z[k] + g_x_0_0_yyzz[k];

                g_x_0_z_yzz[k] = -g_x_0_0_yzz[k] * cd_z[k] + g_x_0_0_yzzz[k];

                g_x_0_z_zzz[k] = -g_x_0_0_zzz[k] * cd_z[k] + g_x_0_0_zzzz[k];
            }
            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_xxx = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 0);

            auto g_y_0_x_xxy = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 1);

            auto g_y_0_x_xxz = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 2);

            auto g_y_0_x_xyy = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 3);

            auto g_y_0_x_xyz = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 4);

            auto g_y_0_x_xzz = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 5);

            auto g_y_0_x_yyy = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 6);

            auto g_y_0_x_yyz = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 7);

            auto g_y_0_x_yzz = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 8);

            auto g_y_0_x_zzz = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_x, g_y_0_0_xxx, g_y_0_0_xxxx, g_y_0_0_xxxy, g_y_0_0_xxxz, g_y_0_0_xxy, g_y_0_0_xxyy, g_y_0_0_xxyz, g_y_0_0_xxz, g_y_0_0_xxzz, g_y_0_0_xyy, g_y_0_0_xyyy, g_y_0_0_xyyz, g_y_0_0_xyz, g_y_0_0_xyzz, g_y_0_0_xzz, g_y_0_0_xzzz, g_y_0_0_yyy, g_y_0_0_yyz, g_y_0_0_yzz, g_y_0_0_zzz, g_y_0_x_xxx, g_y_0_x_xxy, g_y_0_x_xxz, g_y_0_x_xyy, g_y_0_x_xyz, g_y_0_x_xzz, g_y_0_x_yyy, g_y_0_x_yyz, g_y_0_x_yzz, g_y_0_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_xxx[k] = -g_y_0_0_xxx[k] * cd_x[k] + g_y_0_0_xxxx[k];

                g_y_0_x_xxy[k] = -g_y_0_0_xxy[k] * cd_x[k] + g_y_0_0_xxxy[k];

                g_y_0_x_xxz[k] = -g_y_0_0_xxz[k] * cd_x[k] + g_y_0_0_xxxz[k];

                g_y_0_x_xyy[k] = -g_y_0_0_xyy[k] * cd_x[k] + g_y_0_0_xxyy[k];

                g_y_0_x_xyz[k] = -g_y_0_0_xyz[k] * cd_x[k] + g_y_0_0_xxyz[k];

                g_y_0_x_xzz[k] = -g_y_0_0_xzz[k] * cd_x[k] + g_y_0_0_xxzz[k];

                g_y_0_x_yyy[k] = -g_y_0_0_yyy[k] * cd_x[k] + g_y_0_0_xyyy[k];

                g_y_0_x_yyz[k] = -g_y_0_0_yyz[k] * cd_x[k] + g_y_0_0_xyyz[k];

                g_y_0_x_yzz[k] = -g_y_0_0_yzz[k] * cd_x[k] + g_y_0_0_xyzz[k];

                g_y_0_x_zzz[k] = -g_y_0_0_zzz[k] * cd_x[k] + g_y_0_0_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_xxx = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 10);

            auto g_y_0_y_xxy = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 11);

            auto g_y_0_y_xxz = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 12);

            auto g_y_0_y_xyy = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 13);

            auto g_y_0_y_xyz = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 14);

            auto g_y_0_y_xzz = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 15);

            auto g_y_0_y_yyy = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 16);

            auto g_y_0_y_yyz = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 17);

            auto g_y_0_y_yzz = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 18);

            auto g_y_0_y_zzz = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 19);

            #pragma omp simd aligned(cd_y, g_0_xxx, g_0_xxy, g_0_xxz, g_0_xyy, g_0_xyz, g_0_xzz, g_0_yyy, g_0_yyz, g_0_yzz, g_0_zzz, g_y_0_0_xxx, g_y_0_0_xxxy, g_y_0_0_xxy, g_y_0_0_xxyy, g_y_0_0_xxyz, g_y_0_0_xxz, g_y_0_0_xyy, g_y_0_0_xyyy, g_y_0_0_xyyz, g_y_0_0_xyz, g_y_0_0_xyzz, g_y_0_0_xzz, g_y_0_0_yyy, g_y_0_0_yyyy, g_y_0_0_yyyz, g_y_0_0_yyz, g_y_0_0_yyzz, g_y_0_0_yzz, g_y_0_0_yzzz, g_y_0_0_zzz, g_y_0_y_xxx, g_y_0_y_xxy, g_y_0_y_xxz, g_y_0_y_xyy, g_y_0_y_xyz, g_y_0_y_xzz, g_y_0_y_yyy, g_y_0_y_yyz, g_y_0_y_yzz, g_y_0_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_xxx[k] = -g_0_xxx[k] - g_y_0_0_xxx[k] * cd_y[k] + g_y_0_0_xxxy[k];

                g_y_0_y_xxy[k] = -g_0_xxy[k] - g_y_0_0_xxy[k] * cd_y[k] + g_y_0_0_xxyy[k];

                g_y_0_y_xxz[k] = -g_0_xxz[k] - g_y_0_0_xxz[k] * cd_y[k] + g_y_0_0_xxyz[k];

                g_y_0_y_xyy[k] = -g_0_xyy[k] - g_y_0_0_xyy[k] * cd_y[k] + g_y_0_0_xyyy[k];

                g_y_0_y_xyz[k] = -g_0_xyz[k] - g_y_0_0_xyz[k] * cd_y[k] + g_y_0_0_xyyz[k];

                g_y_0_y_xzz[k] = -g_0_xzz[k] - g_y_0_0_xzz[k] * cd_y[k] + g_y_0_0_xyzz[k];

                g_y_0_y_yyy[k] = -g_0_yyy[k] - g_y_0_0_yyy[k] * cd_y[k] + g_y_0_0_yyyy[k];

                g_y_0_y_yyz[k] = -g_0_yyz[k] - g_y_0_0_yyz[k] * cd_y[k] + g_y_0_0_yyyz[k];

                g_y_0_y_yzz[k] = -g_0_yzz[k] - g_y_0_0_yzz[k] * cd_y[k] + g_y_0_0_yyzz[k];

                g_y_0_y_zzz[k] = -g_0_zzz[k] - g_y_0_0_zzz[k] * cd_y[k] + g_y_0_0_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_xxx = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 20);

            auto g_y_0_z_xxy = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 21);

            auto g_y_0_z_xxz = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 22);

            auto g_y_0_z_xyy = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 23);

            auto g_y_0_z_xyz = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 24);

            auto g_y_0_z_xzz = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 25);

            auto g_y_0_z_yyy = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 26);

            auto g_y_0_z_yyz = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 27);

            auto g_y_0_z_yzz = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 28);

            auto g_y_0_z_zzz = cbuffer.data(pf_geom_10_off + 30 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_z, g_y_0_0_xxx, g_y_0_0_xxxz, g_y_0_0_xxy, g_y_0_0_xxyz, g_y_0_0_xxz, g_y_0_0_xxzz, g_y_0_0_xyy, g_y_0_0_xyyz, g_y_0_0_xyz, g_y_0_0_xyzz, g_y_0_0_xzz, g_y_0_0_xzzz, g_y_0_0_yyy, g_y_0_0_yyyz, g_y_0_0_yyz, g_y_0_0_yyzz, g_y_0_0_yzz, g_y_0_0_yzzz, g_y_0_0_zzz, g_y_0_0_zzzz, g_y_0_z_xxx, g_y_0_z_xxy, g_y_0_z_xxz, g_y_0_z_xyy, g_y_0_z_xyz, g_y_0_z_xzz, g_y_0_z_yyy, g_y_0_z_yyz, g_y_0_z_yzz, g_y_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_xxx[k] = -g_y_0_0_xxx[k] * cd_z[k] + g_y_0_0_xxxz[k];

                g_y_0_z_xxy[k] = -g_y_0_0_xxy[k] * cd_z[k] + g_y_0_0_xxyz[k];

                g_y_0_z_xxz[k] = -g_y_0_0_xxz[k] * cd_z[k] + g_y_0_0_xxzz[k];

                g_y_0_z_xyy[k] = -g_y_0_0_xyy[k] * cd_z[k] + g_y_0_0_xyyz[k];

                g_y_0_z_xyz[k] = -g_y_0_0_xyz[k] * cd_z[k] + g_y_0_0_xyzz[k];

                g_y_0_z_xzz[k] = -g_y_0_0_xzz[k] * cd_z[k] + g_y_0_0_xzzz[k];

                g_y_0_z_yyy[k] = -g_y_0_0_yyy[k] * cd_z[k] + g_y_0_0_yyyz[k];

                g_y_0_z_yyz[k] = -g_y_0_0_yyz[k] * cd_z[k] + g_y_0_0_yyzz[k];

                g_y_0_z_yzz[k] = -g_y_0_0_yzz[k] * cd_z[k] + g_y_0_0_yzzz[k];

                g_y_0_z_zzz[k] = -g_y_0_0_zzz[k] * cd_z[k] + g_y_0_0_zzzz[k];
            }
            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_xxx = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 0);

            auto g_z_0_x_xxy = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 1);

            auto g_z_0_x_xxz = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 2);

            auto g_z_0_x_xyy = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 3);

            auto g_z_0_x_xyz = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 4);

            auto g_z_0_x_xzz = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 5);

            auto g_z_0_x_yyy = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 6);

            auto g_z_0_x_yyz = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 7);

            auto g_z_0_x_yzz = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 8);

            auto g_z_0_x_zzz = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_x, g_z_0_0_xxx, g_z_0_0_xxxx, g_z_0_0_xxxy, g_z_0_0_xxxz, g_z_0_0_xxy, g_z_0_0_xxyy, g_z_0_0_xxyz, g_z_0_0_xxz, g_z_0_0_xxzz, g_z_0_0_xyy, g_z_0_0_xyyy, g_z_0_0_xyyz, g_z_0_0_xyz, g_z_0_0_xyzz, g_z_0_0_xzz, g_z_0_0_xzzz, g_z_0_0_yyy, g_z_0_0_yyz, g_z_0_0_yzz, g_z_0_0_zzz, g_z_0_x_xxx, g_z_0_x_xxy, g_z_0_x_xxz, g_z_0_x_xyy, g_z_0_x_xyz, g_z_0_x_xzz, g_z_0_x_yyy, g_z_0_x_yyz, g_z_0_x_yzz, g_z_0_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_xxx[k] = -g_z_0_0_xxx[k] * cd_x[k] + g_z_0_0_xxxx[k];

                g_z_0_x_xxy[k] = -g_z_0_0_xxy[k] * cd_x[k] + g_z_0_0_xxxy[k];

                g_z_0_x_xxz[k] = -g_z_0_0_xxz[k] * cd_x[k] + g_z_0_0_xxxz[k];

                g_z_0_x_xyy[k] = -g_z_0_0_xyy[k] * cd_x[k] + g_z_0_0_xxyy[k];

                g_z_0_x_xyz[k] = -g_z_0_0_xyz[k] * cd_x[k] + g_z_0_0_xxyz[k];

                g_z_0_x_xzz[k] = -g_z_0_0_xzz[k] * cd_x[k] + g_z_0_0_xxzz[k];

                g_z_0_x_yyy[k] = -g_z_0_0_yyy[k] * cd_x[k] + g_z_0_0_xyyy[k];

                g_z_0_x_yyz[k] = -g_z_0_0_yyz[k] * cd_x[k] + g_z_0_0_xyyz[k];

                g_z_0_x_yzz[k] = -g_z_0_0_yzz[k] * cd_x[k] + g_z_0_0_xyzz[k];

                g_z_0_x_zzz[k] = -g_z_0_0_zzz[k] * cd_x[k] + g_z_0_0_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_xxx = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 10);

            auto g_z_0_y_xxy = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 11);

            auto g_z_0_y_xxz = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 12);

            auto g_z_0_y_xyy = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 13);

            auto g_z_0_y_xyz = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 14);

            auto g_z_0_y_xzz = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 15);

            auto g_z_0_y_yyy = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 16);

            auto g_z_0_y_yyz = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 17);

            auto g_z_0_y_yzz = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 18);

            auto g_z_0_y_zzz = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 19);

            #pragma omp simd aligned(cd_y, g_z_0_0_xxx, g_z_0_0_xxxy, g_z_0_0_xxy, g_z_0_0_xxyy, g_z_0_0_xxyz, g_z_0_0_xxz, g_z_0_0_xyy, g_z_0_0_xyyy, g_z_0_0_xyyz, g_z_0_0_xyz, g_z_0_0_xyzz, g_z_0_0_xzz, g_z_0_0_yyy, g_z_0_0_yyyy, g_z_0_0_yyyz, g_z_0_0_yyz, g_z_0_0_yyzz, g_z_0_0_yzz, g_z_0_0_yzzz, g_z_0_0_zzz, g_z_0_y_xxx, g_z_0_y_xxy, g_z_0_y_xxz, g_z_0_y_xyy, g_z_0_y_xyz, g_z_0_y_xzz, g_z_0_y_yyy, g_z_0_y_yyz, g_z_0_y_yzz, g_z_0_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_xxx[k] = -g_z_0_0_xxx[k] * cd_y[k] + g_z_0_0_xxxy[k];

                g_z_0_y_xxy[k] = -g_z_0_0_xxy[k] * cd_y[k] + g_z_0_0_xxyy[k];

                g_z_0_y_xxz[k] = -g_z_0_0_xxz[k] * cd_y[k] + g_z_0_0_xxyz[k];

                g_z_0_y_xyy[k] = -g_z_0_0_xyy[k] * cd_y[k] + g_z_0_0_xyyy[k];

                g_z_0_y_xyz[k] = -g_z_0_0_xyz[k] * cd_y[k] + g_z_0_0_xyyz[k];

                g_z_0_y_xzz[k] = -g_z_0_0_xzz[k] * cd_y[k] + g_z_0_0_xyzz[k];

                g_z_0_y_yyy[k] = -g_z_0_0_yyy[k] * cd_y[k] + g_z_0_0_yyyy[k];

                g_z_0_y_yyz[k] = -g_z_0_0_yyz[k] * cd_y[k] + g_z_0_0_yyyz[k];

                g_z_0_y_yzz[k] = -g_z_0_0_yzz[k] * cd_y[k] + g_z_0_0_yyzz[k];

                g_z_0_y_zzz[k] = -g_z_0_0_zzz[k] * cd_y[k] + g_z_0_0_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_xxx = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 20);

            auto g_z_0_z_xxy = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 21);

            auto g_z_0_z_xxz = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 22);

            auto g_z_0_z_xyy = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 23);

            auto g_z_0_z_xyz = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 24);

            auto g_z_0_z_xzz = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 25);

            auto g_z_0_z_yyy = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 26);

            auto g_z_0_z_yyz = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 27);

            auto g_z_0_z_yzz = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 28);

            auto g_z_0_z_zzz = cbuffer.data(pf_geom_10_off + 60 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_z, g_0_xxx, g_0_xxy, g_0_xxz, g_0_xyy, g_0_xyz, g_0_xzz, g_0_yyy, g_0_yyz, g_0_yzz, g_0_zzz, g_z_0_0_xxx, g_z_0_0_xxxz, g_z_0_0_xxy, g_z_0_0_xxyz, g_z_0_0_xxz, g_z_0_0_xxzz, g_z_0_0_xyy, g_z_0_0_xyyz, g_z_0_0_xyz, g_z_0_0_xyzz, g_z_0_0_xzz, g_z_0_0_xzzz, g_z_0_0_yyy, g_z_0_0_yyyz, g_z_0_0_yyz, g_z_0_0_yyzz, g_z_0_0_yzz, g_z_0_0_yzzz, g_z_0_0_zzz, g_z_0_0_zzzz, g_z_0_z_xxx, g_z_0_z_xxy, g_z_0_z_xxz, g_z_0_z_xyy, g_z_0_z_xyz, g_z_0_z_xzz, g_z_0_z_yyy, g_z_0_z_yyz, g_z_0_z_yzz, g_z_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_xxx[k] = -g_0_xxx[k] - g_z_0_0_xxx[k] * cd_z[k] + g_z_0_0_xxxz[k];

                g_z_0_z_xxy[k] = -g_0_xxy[k] - g_z_0_0_xxy[k] * cd_z[k] + g_z_0_0_xxyz[k];

                g_z_0_z_xxz[k] = -g_0_xxz[k] - g_z_0_0_xxz[k] * cd_z[k] + g_z_0_0_xxzz[k];

                g_z_0_z_xyy[k] = -g_0_xyy[k] - g_z_0_0_xyy[k] * cd_z[k] + g_z_0_0_xyyz[k];

                g_z_0_z_xyz[k] = -g_0_xyz[k] - g_z_0_0_xyz[k] * cd_z[k] + g_z_0_0_xyzz[k];

                g_z_0_z_xzz[k] = -g_0_xzz[k] - g_z_0_0_xzz[k] * cd_z[k] + g_z_0_0_xzzz[k];

                g_z_0_z_yyy[k] = -g_0_yyy[k] - g_z_0_0_yyy[k] * cd_z[k] + g_z_0_0_yyyz[k];

                g_z_0_z_yyz[k] = -g_0_yyz[k] - g_z_0_0_yyz[k] * cd_z[k] + g_z_0_0_yyzz[k];

                g_z_0_z_yzz[k] = -g_0_yzz[k] - g_z_0_0_yzz[k] * cd_z[k] + g_z_0_0_yzzz[k];

                g_z_0_z_zzz[k] = -g_0_zzz[k] - g_z_0_0_zzz[k] * cd_z[k] + g_z_0_0_zzzz[k];
            }
        }
    }
}

} // erirec namespace

