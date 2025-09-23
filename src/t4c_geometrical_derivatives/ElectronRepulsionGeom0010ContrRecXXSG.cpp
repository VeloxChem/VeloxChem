#include "ElectronRepulsionGeom0010ContrRecXXSG.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxsg(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxsg,
                                            CSimdArray<double>& pbuffer,
                                            const size_t idx_xxsg,
                                            const size_t idx_xxsh,
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

            /// Set up components of auxilary buffer : SSSH

            const auto sh_off = idx_xxsh + (i * bcomps + j) * 21;

            auto g_0_xxxxx = pbuffer.data(sh_off + 0);

            auto g_0_xxxxy = pbuffer.data(sh_off + 1);

            auto g_0_xxxxz = pbuffer.data(sh_off + 2);

            auto g_0_xxxyy = pbuffer.data(sh_off + 3);

            auto g_0_xxxyz = pbuffer.data(sh_off + 4);

            auto g_0_xxxzz = pbuffer.data(sh_off + 5);

            auto g_0_xxyyy = pbuffer.data(sh_off + 6);

            auto g_0_xxyyz = pbuffer.data(sh_off + 7);

            auto g_0_xxyzz = pbuffer.data(sh_off + 8);

            auto g_0_xxzzz = pbuffer.data(sh_off + 9);

            auto g_0_xyyyy = pbuffer.data(sh_off + 10);

            auto g_0_xyyyz = pbuffer.data(sh_off + 11);

            auto g_0_xyyzz = pbuffer.data(sh_off + 12);

            auto g_0_xyzzz = pbuffer.data(sh_off + 13);

            auto g_0_xzzzz = pbuffer.data(sh_off + 14);

            auto g_0_yyyyy = pbuffer.data(sh_off + 15);

            auto g_0_yyyyz = pbuffer.data(sh_off + 16);

            auto g_0_yyyzz = pbuffer.data(sh_off + 17);

            auto g_0_yyzzz = pbuffer.data(sh_off + 18);

            auto g_0_yzzzz = pbuffer.data(sh_off + 19);

            auto g_0_zzzzz = pbuffer.data(sh_off + 20);

            /// set up bra offset for contr_buffer_xxsg

            const auto sg_geom_10_off = idx_geom_10_xxsg + (i * bcomps + j) * 15;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_0_xxxx, g_0_xxxxx, g_0_xxxxy, g_0_xxxxz, g_0_xxxy, g_0_xxxyy, g_0_xxxyz, g_0_xxxz, g_0_xxxzz, g_0_xxyy, g_0_xxyyy, g_0_xxyyz, g_0_xxyz, g_0_xxyzz, g_0_xxzz, g_0_xxzzz, g_0_xyyy, g_0_xyyyy, g_0_xyyyz, g_0_xyyz, g_0_xyyzz, g_0_xyzz, g_0_xyzzz, g_0_xzzz, g_0_xzzzz, g_0_yyyy, g_0_yyyz, g_0_yyzz, g_0_yzzz, g_0_zzzz, g_x_0_0_xxxx, g_x_0_0_xxxy, g_x_0_0_xxxz, g_x_0_0_xxyy, g_x_0_0_xxyz, g_x_0_0_xxzz, g_x_0_0_xyyy, g_x_0_0_xyyz, g_x_0_0_xyzz, g_x_0_0_xzzz, g_x_0_0_yyyy, g_x_0_0_yyyz, g_x_0_0_yyzz, g_x_0_0_yzzz, g_x_0_0_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_0_xxxx[k] = -g_0_xxxx[k] * cd_x[k] + g_0_xxxxx[k];

                g_x_0_0_xxxy[k] = -g_0_xxxy[k] * cd_x[k] + g_0_xxxxy[k];

                g_x_0_0_xxxz[k] = -g_0_xxxz[k] * cd_x[k] + g_0_xxxxz[k];

                g_x_0_0_xxyy[k] = -g_0_xxyy[k] * cd_x[k] + g_0_xxxyy[k];

                g_x_0_0_xxyz[k] = -g_0_xxyz[k] * cd_x[k] + g_0_xxxyz[k];

                g_x_0_0_xxzz[k] = -g_0_xxzz[k] * cd_x[k] + g_0_xxxzz[k];

                g_x_0_0_xyyy[k] = -g_0_xyyy[k] * cd_x[k] + g_0_xxyyy[k];

                g_x_0_0_xyyz[k] = -g_0_xyyz[k] * cd_x[k] + g_0_xxyyz[k];

                g_x_0_0_xyzz[k] = -g_0_xyzz[k] * cd_x[k] + g_0_xxyzz[k];

                g_x_0_0_xzzz[k] = -g_0_xzzz[k] * cd_x[k] + g_0_xxzzz[k];

                g_x_0_0_yyyy[k] = -g_0_yyyy[k] * cd_x[k] + g_0_xyyyy[k];

                g_x_0_0_yyyz[k] = -g_0_yyyz[k] * cd_x[k] + g_0_xyyyz[k];

                g_x_0_0_yyzz[k] = -g_0_yyzz[k] * cd_x[k] + g_0_xyyzz[k];

                g_x_0_0_yzzz[k] = -g_0_yzzz[k] * cd_x[k] + g_0_xyzzz[k];

                g_x_0_0_zzzz[k] = -g_0_zzzz[k] * cd_x[k] + g_0_xzzzz[k];
            }
            /// Set up 0-15 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_0_xxxx, g_0_xxxxy, g_0_xxxy, g_0_xxxyy, g_0_xxxyz, g_0_xxxz, g_0_xxyy, g_0_xxyyy, g_0_xxyyz, g_0_xxyz, g_0_xxyzz, g_0_xxzz, g_0_xyyy, g_0_xyyyy, g_0_xyyyz, g_0_xyyz, g_0_xyyzz, g_0_xyzz, g_0_xyzzz, g_0_xzzz, g_0_yyyy, g_0_yyyyy, g_0_yyyyz, g_0_yyyz, g_0_yyyzz, g_0_yyzz, g_0_yyzzz, g_0_yzzz, g_0_yzzzz, g_0_zzzz, g_y_0_0_xxxx, g_y_0_0_xxxy, g_y_0_0_xxxz, g_y_0_0_xxyy, g_y_0_0_xxyz, g_y_0_0_xxzz, g_y_0_0_xyyy, g_y_0_0_xyyz, g_y_0_0_xyzz, g_y_0_0_xzzz, g_y_0_0_yyyy, g_y_0_0_yyyz, g_y_0_0_yyzz, g_y_0_0_yzzz, g_y_0_0_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_0_xxxx[k] = -g_0_xxxx[k] * cd_y[k] + g_0_xxxxy[k];

                g_y_0_0_xxxy[k] = -g_0_xxxy[k] * cd_y[k] + g_0_xxxyy[k];

                g_y_0_0_xxxz[k] = -g_0_xxxz[k] * cd_y[k] + g_0_xxxyz[k];

                g_y_0_0_xxyy[k] = -g_0_xxyy[k] * cd_y[k] + g_0_xxyyy[k];

                g_y_0_0_xxyz[k] = -g_0_xxyz[k] * cd_y[k] + g_0_xxyyz[k];

                g_y_0_0_xxzz[k] = -g_0_xxzz[k] * cd_y[k] + g_0_xxyzz[k];

                g_y_0_0_xyyy[k] = -g_0_xyyy[k] * cd_y[k] + g_0_xyyyy[k];

                g_y_0_0_xyyz[k] = -g_0_xyyz[k] * cd_y[k] + g_0_xyyyz[k];

                g_y_0_0_xyzz[k] = -g_0_xyzz[k] * cd_y[k] + g_0_xyyzz[k];

                g_y_0_0_xzzz[k] = -g_0_xzzz[k] * cd_y[k] + g_0_xyzzz[k];

                g_y_0_0_yyyy[k] = -g_0_yyyy[k] * cd_y[k] + g_0_yyyyy[k];

                g_y_0_0_yyyz[k] = -g_0_yyyz[k] * cd_y[k] + g_0_yyyyz[k];

                g_y_0_0_yyzz[k] = -g_0_yyzz[k] * cd_y[k] + g_0_yyyzz[k];

                g_y_0_0_yzzz[k] = -g_0_yzzz[k] * cd_y[k] + g_0_yyzzz[k];

                g_y_0_0_zzzz[k] = -g_0_zzzz[k] * cd_y[k] + g_0_yzzzz[k];
            }
            /// Set up 0-15 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_z, g_0_xxxx, g_0_xxxxz, g_0_xxxy, g_0_xxxyz, g_0_xxxz, g_0_xxxzz, g_0_xxyy, g_0_xxyyz, g_0_xxyz, g_0_xxyzz, g_0_xxzz, g_0_xxzzz, g_0_xyyy, g_0_xyyyz, g_0_xyyz, g_0_xyyzz, g_0_xyzz, g_0_xyzzz, g_0_xzzz, g_0_xzzzz, g_0_yyyy, g_0_yyyyz, g_0_yyyz, g_0_yyyzz, g_0_yyzz, g_0_yyzzz, g_0_yzzz, g_0_yzzzz, g_0_zzzz, g_0_zzzzz, g_z_0_0_xxxx, g_z_0_0_xxxy, g_z_0_0_xxxz, g_z_0_0_xxyy, g_z_0_0_xxyz, g_z_0_0_xxzz, g_z_0_0_xyyy, g_z_0_0_xyyz, g_z_0_0_xyzz, g_z_0_0_xzzz, g_z_0_0_yyyy, g_z_0_0_yyyz, g_z_0_0_yyzz, g_z_0_0_yzzz, g_z_0_0_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_0_xxxx[k] = -g_0_xxxx[k] * cd_z[k] + g_0_xxxxz[k];

                g_z_0_0_xxxy[k] = -g_0_xxxy[k] * cd_z[k] + g_0_xxxyz[k];

                g_z_0_0_xxxz[k] = -g_0_xxxz[k] * cd_z[k] + g_0_xxxzz[k];

                g_z_0_0_xxyy[k] = -g_0_xxyy[k] * cd_z[k] + g_0_xxyyz[k];

                g_z_0_0_xxyz[k] = -g_0_xxyz[k] * cd_z[k] + g_0_xxyzz[k];

                g_z_0_0_xxzz[k] = -g_0_xxzz[k] * cd_z[k] + g_0_xxzzz[k];

                g_z_0_0_xyyy[k] = -g_0_xyyy[k] * cd_z[k] + g_0_xyyyz[k];

                g_z_0_0_xyyz[k] = -g_0_xyyz[k] * cd_z[k] + g_0_xyyzz[k];

                g_z_0_0_xyzz[k] = -g_0_xyzz[k] * cd_z[k] + g_0_xyzzz[k];

                g_z_0_0_xzzz[k] = -g_0_xzzz[k] * cd_z[k] + g_0_xzzzz[k];

                g_z_0_0_yyyy[k] = -g_0_yyyy[k] * cd_z[k] + g_0_yyyyz[k];

                g_z_0_0_yyyz[k] = -g_0_yyyz[k] * cd_z[k] + g_0_yyyzz[k];

                g_z_0_0_yyzz[k] = -g_0_yyzz[k] * cd_z[k] + g_0_yyzzz[k];

                g_z_0_0_yzzz[k] = -g_0_yzzz[k] * cd_z[k] + g_0_yzzzz[k];

                g_z_0_0_zzzz[k] = -g_0_zzzz[k] * cd_z[k] + g_0_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

