#include "ElectronRepulsionGeom0010ContrRecXXSF.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxsf(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxsf,
                                            CSimdArray<double>& pbuffer,
                                            const size_t idx_xxsf,
                                            const size_t idx_xxsg,
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

            /// Set up components of auxilary buffer : SSSG

            const auto sg_off = idx_xxsg + (i * bcomps + j) * 15;

            auto g_0_xxxx = pbuffer.data(sg_off + 0);

            auto g_0_xxxy = pbuffer.data(sg_off + 1);

            auto g_0_xxxz = pbuffer.data(sg_off + 2);

            auto g_0_xxyy = pbuffer.data(sg_off + 3);

            auto g_0_xxyz = pbuffer.data(sg_off + 4);

            auto g_0_xxzz = pbuffer.data(sg_off + 5);

            auto g_0_xyyy = pbuffer.data(sg_off + 6);

            auto g_0_xyyz = pbuffer.data(sg_off + 7);

            auto g_0_xyzz = pbuffer.data(sg_off + 8);

            auto g_0_xzzz = pbuffer.data(sg_off + 9);

            auto g_0_yyyy = pbuffer.data(sg_off + 10);

            auto g_0_yyyz = pbuffer.data(sg_off + 11);

            auto g_0_yyzz = pbuffer.data(sg_off + 12);

            auto g_0_yzzz = pbuffer.data(sg_off + 13);

            auto g_0_zzzz = pbuffer.data(sg_off + 14);

            /// set up bra offset for contr_buffer_xxsf

            const auto sf_geom_10_off = idx_geom_10_xxsf + (i * bcomps + j) * 10;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_0_xxx, g_0_xxxx, g_0_xxxy, g_0_xxxz, g_0_xxy, g_0_xxyy, g_0_xxyz, g_0_xxz, g_0_xxzz, g_0_xyy, g_0_xyyy, g_0_xyyz, g_0_xyz, g_0_xyzz, g_0_xzz, g_0_xzzz, g_0_yyy, g_0_yyz, g_0_yzz, g_0_zzz, g_x_0_0_xxx, g_x_0_0_xxy, g_x_0_0_xxz, g_x_0_0_xyy, g_x_0_0_xyz, g_x_0_0_xzz, g_x_0_0_yyy, g_x_0_0_yyz, g_x_0_0_yzz, g_x_0_0_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_0_xxx[k] = -g_0_xxx[k] * cd_x[k] + g_0_xxxx[k];

                g_x_0_0_xxy[k] = -g_0_xxy[k] * cd_x[k] + g_0_xxxy[k];

                g_x_0_0_xxz[k] = -g_0_xxz[k] * cd_x[k] + g_0_xxxz[k];

                g_x_0_0_xyy[k] = -g_0_xyy[k] * cd_x[k] + g_0_xxyy[k];

                g_x_0_0_xyz[k] = -g_0_xyz[k] * cd_x[k] + g_0_xxyz[k];

                g_x_0_0_xzz[k] = -g_0_xzz[k] * cd_x[k] + g_0_xxzz[k];

                g_x_0_0_yyy[k] = -g_0_yyy[k] * cd_x[k] + g_0_xyyy[k];

                g_x_0_0_yyz[k] = -g_0_yyz[k] * cd_x[k] + g_0_xyyz[k];

                g_x_0_0_yzz[k] = -g_0_yzz[k] * cd_x[k] + g_0_xyzz[k];

                g_x_0_0_zzz[k] = -g_0_zzz[k] * cd_x[k] + g_0_xzzz[k];
            }
            /// Set up 0-10 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_0_xxx, g_0_xxxy, g_0_xxy, g_0_xxyy, g_0_xxyz, g_0_xxz, g_0_xyy, g_0_xyyy, g_0_xyyz, g_0_xyz, g_0_xyzz, g_0_xzz, g_0_yyy, g_0_yyyy, g_0_yyyz, g_0_yyz, g_0_yyzz, g_0_yzz, g_0_yzzz, g_0_zzz, g_y_0_0_xxx, g_y_0_0_xxy, g_y_0_0_xxz, g_y_0_0_xyy, g_y_0_0_xyz, g_y_0_0_xzz, g_y_0_0_yyy, g_y_0_0_yyz, g_y_0_0_yzz, g_y_0_0_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_0_xxx[k] = -g_0_xxx[k] * cd_y[k] + g_0_xxxy[k];

                g_y_0_0_xxy[k] = -g_0_xxy[k] * cd_y[k] + g_0_xxyy[k];

                g_y_0_0_xxz[k] = -g_0_xxz[k] * cd_y[k] + g_0_xxyz[k];

                g_y_0_0_xyy[k] = -g_0_xyy[k] * cd_y[k] + g_0_xyyy[k];

                g_y_0_0_xyz[k] = -g_0_xyz[k] * cd_y[k] + g_0_xyyz[k];

                g_y_0_0_xzz[k] = -g_0_xzz[k] * cd_y[k] + g_0_xyzz[k];

                g_y_0_0_yyy[k] = -g_0_yyy[k] * cd_y[k] + g_0_yyyy[k];

                g_y_0_0_yyz[k] = -g_0_yyz[k] * cd_y[k] + g_0_yyyz[k];

                g_y_0_0_yzz[k] = -g_0_yzz[k] * cd_y[k] + g_0_yyzz[k];

                g_y_0_0_zzz[k] = -g_0_zzz[k] * cd_y[k] + g_0_yzzz[k];
            }
            /// Set up 0-10 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_z, g_0_xxx, g_0_xxxz, g_0_xxy, g_0_xxyz, g_0_xxz, g_0_xxzz, g_0_xyy, g_0_xyyz, g_0_xyz, g_0_xyzz, g_0_xzz, g_0_xzzz, g_0_yyy, g_0_yyyz, g_0_yyz, g_0_yyzz, g_0_yzz, g_0_yzzz, g_0_zzz, g_0_zzzz, g_z_0_0_xxx, g_z_0_0_xxy, g_z_0_0_xxz, g_z_0_0_xyy, g_z_0_0_xyz, g_z_0_0_xzz, g_z_0_0_yyy, g_z_0_0_yyz, g_z_0_0_yzz, g_z_0_0_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_0_xxx[k] = -g_0_xxx[k] * cd_z[k] + g_0_xxxz[k];

                g_z_0_0_xxy[k] = -g_0_xxy[k] * cd_z[k] + g_0_xxyz[k];

                g_z_0_0_xxz[k] = -g_0_xxz[k] * cd_z[k] + g_0_xxzz[k];

                g_z_0_0_xyy[k] = -g_0_xyy[k] * cd_z[k] + g_0_xyyz[k];

                g_z_0_0_xyz[k] = -g_0_xyz[k] * cd_z[k] + g_0_xyzz[k];

                g_z_0_0_xzz[k] = -g_0_xzz[k] * cd_z[k] + g_0_xzzz[k];

                g_z_0_0_yyy[k] = -g_0_yyy[k] * cd_z[k] + g_0_yyyz[k];

                g_z_0_0_yyz[k] = -g_0_yyz[k] * cd_z[k] + g_0_yyzz[k];

                g_z_0_0_yzz[k] = -g_0_yzz[k] * cd_z[k] + g_0_yzzz[k];

                g_z_0_0_zzz[k] = -g_0_zzz[k] * cd_z[k] + g_0_zzzz[k];
            }
        }
    }
}

} // erirec namespace

