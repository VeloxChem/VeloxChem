#include "ElectronRepulsionGeom0010ContrRecXXPG.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxpg(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxpg,
                                            CSimdArray<double>& pbuffer,
                                            const size_t idx_xxsg,
                                            const size_t idx_geom_10_xxsg,
                                            const size_t idx_geom_10_xxsh,
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

            /// Set up components of auxilary buffer : SSSH

            const auto sh_geom_10_off = idx_geom_10_xxsh + (i * bcomps + j) * 21;

            auto g_x_0_0_xxxxx = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_0_xxxxy = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_0_xxxxz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_0_xxxyy = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_0_xxxyz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_0_xxxzz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_0_xxyyy = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_0_xxyyz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_0_xxyzz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_0_xxzzz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_0_xyyyy = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_0_xyyyz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_0_xyyzz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_0_xyzzz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_0_xzzzz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_0_yyyyy = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_0_yyyyz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_0_yyyzz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_0_yyzzz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_0_yzzzz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_0_zzzzz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_y_0_0_xxxxx = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 0);

            auto g_y_0_0_xxxxy = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 1);

            auto g_y_0_0_xxxxz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 2);

            auto g_y_0_0_xxxyy = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 3);

            auto g_y_0_0_xxxyz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 4);

            auto g_y_0_0_xxxzz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 5);

            auto g_y_0_0_xxyyy = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 6);

            auto g_y_0_0_xxyyz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 7);

            auto g_y_0_0_xxyzz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 8);

            auto g_y_0_0_xxzzz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 9);

            auto g_y_0_0_xyyyy = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 10);

            auto g_y_0_0_xyyyz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 11);

            auto g_y_0_0_xyyzz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 12);

            auto g_y_0_0_xyzzz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 13);

            auto g_y_0_0_xzzzz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 14);

            auto g_y_0_0_yyyyy = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 15);

            auto g_y_0_0_yyyyz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 16);

            auto g_y_0_0_yyyzz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 17);

            auto g_y_0_0_yyzzz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 18);

            auto g_y_0_0_yzzzz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 19);

            auto g_y_0_0_zzzzz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 20);

            auto g_z_0_0_xxxxx = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 0);

            auto g_z_0_0_xxxxy = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 1);

            auto g_z_0_0_xxxxz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 2);

            auto g_z_0_0_xxxyy = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 3);

            auto g_z_0_0_xxxyz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 4);

            auto g_z_0_0_xxxzz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 5);

            auto g_z_0_0_xxyyy = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 6);

            auto g_z_0_0_xxyyz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 7);

            auto g_z_0_0_xxyzz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 8);

            auto g_z_0_0_xxzzz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 9);

            auto g_z_0_0_xyyyy = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 10);

            auto g_z_0_0_xyyyz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 11);

            auto g_z_0_0_xyyzz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 12);

            auto g_z_0_0_xyzzz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 13);

            auto g_z_0_0_xzzzz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 14);

            auto g_z_0_0_yyyyy = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 15);

            auto g_z_0_0_yyyyz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 16);

            auto g_z_0_0_yyyzz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 17);

            auto g_z_0_0_yyzzz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 18);

            auto g_z_0_0_yzzzz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 19);

            auto g_z_0_0_zzzzz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 20);

            /// set up bra offset for contr_buffer_xxpg

            const auto pg_geom_10_off = idx_geom_10_xxpg + (i * bcomps + j) * 45;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_xxxx = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_x_xxxy = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_x_xxxz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_x_xxyy = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_x_xxyz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_x_xxzz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_x_xyyy = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_x_xyyz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_x_xyzz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_x_xzzz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_x_yyyy = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_x_yyyz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_x_yyzz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_x_yzzz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_x_zzzz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_x, g_0_xxxx, g_0_xxxy, g_0_xxxz, g_0_xxyy, g_0_xxyz, g_0_xxzz, g_0_xyyy, g_0_xyyz, g_0_xyzz, g_0_xzzz, g_0_yyyy, g_0_yyyz, g_0_yyzz, g_0_yzzz, g_0_zzzz, g_x_0_0_xxxx, g_x_0_0_xxxxx, g_x_0_0_xxxxy, g_x_0_0_xxxxz, g_x_0_0_xxxy, g_x_0_0_xxxyy, g_x_0_0_xxxyz, g_x_0_0_xxxz, g_x_0_0_xxxzz, g_x_0_0_xxyy, g_x_0_0_xxyyy, g_x_0_0_xxyyz, g_x_0_0_xxyz, g_x_0_0_xxyzz, g_x_0_0_xxzz, g_x_0_0_xxzzz, g_x_0_0_xyyy, g_x_0_0_xyyyy, g_x_0_0_xyyyz, g_x_0_0_xyyz, g_x_0_0_xyyzz, g_x_0_0_xyzz, g_x_0_0_xyzzz, g_x_0_0_xzzz, g_x_0_0_xzzzz, g_x_0_0_yyyy, g_x_0_0_yyyz, g_x_0_0_yyzz, g_x_0_0_yzzz, g_x_0_0_zzzz, g_x_0_x_xxxx, g_x_0_x_xxxy, g_x_0_x_xxxz, g_x_0_x_xxyy, g_x_0_x_xxyz, g_x_0_x_xxzz, g_x_0_x_xyyy, g_x_0_x_xyyz, g_x_0_x_xyzz, g_x_0_x_xzzz, g_x_0_x_yyyy, g_x_0_x_yyyz, g_x_0_x_yyzz, g_x_0_x_yzzz, g_x_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_xxxx[k] = -g_0_xxxx[k] - g_x_0_0_xxxx[k] * cd_x[k] + g_x_0_0_xxxxx[k];

                g_x_0_x_xxxy[k] = -g_0_xxxy[k] - g_x_0_0_xxxy[k] * cd_x[k] + g_x_0_0_xxxxy[k];

                g_x_0_x_xxxz[k] = -g_0_xxxz[k] - g_x_0_0_xxxz[k] * cd_x[k] + g_x_0_0_xxxxz[k];

                g_x_0_x_xxyy[k] = -g_0_xxyy[k] - g_x_0_0_xxyy[k] * cd_x[k] + g_x_0_0_xxxyy[k];

                g_x_0_x_xxyz[k] = -g_0_xxyz[k] - g_x_0_0_xxyz[k] * cd_x[k] + g_x_0_0_xxxyz[k];

                g_x_0_x_xxzz[k] = -g_0_xxzz[k] - g_x_0_0_xxzz[k] * cd_x[k] + g_x_0_0_xxxzz[k];

                g_x_0_x_xyyy[k] = -g_0_xyyy[k] - g_x_0_0_xyyy[k] * cd_x[k] + g_x_0_0_xxyyy[k];

                g_x_0_x_xyyz[k] = -g_0_xyyz[k] - g_x_0_0_xyyz[k] * cd_x[k] + g_x_0_0_xxyyz[k];

                g_x_0_x_xyzz[k] = -g_0_xyzz[k] - g_x_0_0_xyzz[k] * cd_x[k] + g_x_0_0_xxyzz[k];

                g_x_0_x_xzzz[k] = -g_0_xzzz[k] - g_x_0_0_xzzz[k] * cd_x[k] + g_x_0_0_xxzzz[k];

                g_x_0_x_yyyy[k] = -g_0_yyyy[k] - g_x_0_0_yyyy[k] * cd_x[k] + g_x_0_0_xyyyy[k];

                g_x_0_x_yyyz[k] = -g_0_yyyz[k] - g_x_0_0_yyyz[k] * cd_x[k] + g_x_0_0_xyyyz[k];

                g_x_0_x_yyzz[k] = -g_0_yyzz[k] - g_x_0_0_yyzz[k] * cd_x[k] + g_x_0_0_xyyzz[k];

                g_x_0_x_yzzz[k] = -g_0_yzzz[k] - g_x_0_0_yzzz[k] * cd_x[k] + g_x_0_0_xyzzz[k];

                g_x_0_x_zzzz[k] = -g_0_zzzz[k] - g_x_0_0_zzzz[k] * cd_x[k] + g_x_0_0_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_xxxx = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_y_xxxy = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_y_xxxz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_y_xxyy = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_y_xxyz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_y_xxzz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_y_xyyy = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_y_xyyz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_y_xyzz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_y_xzzz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_y_yyyy = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_y_yyyz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_y_yyzz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_y_yzzz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_y_zzzz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_y, g_x_0_0_xxxx, g_x_0_0_xxxxy, g_x_0_0_xxxy, g_x_0_0_xxxyy, g_x_0_0_xxxyz, g_x_0_0_xxxz, g_x_0_0_xxyy, g_x_0_0_xxyyy, g_x_0_0_xxyyz, g_x_0_0_xxyz, g_x_0_0_xxyzz, g_x_0_0_xxzz, g_x_0_0_xyyy, g_x_0_0_xyyyy, g_x_0_0_xyyyz, g_x_0_0_xyyz, g_x_0_0_xyyzz, g_x_0_0_xyzz, g_x_0_0_xyzzz, g_x_0_0_xzzz, g_x_0_0_yyyy, g_x_0_0_yyyyy, g_x_0_0_yyyyz, g_x_0_0_yyyz, g_x_0_0_yyyzz, g_x_0_0_yyzz, g_x_0_0_yyzzz, g_x_0_0_yzzz, g_x_0_0_yzzzz, g_x_0_0_zzzz, g_x_0_y_xxxx, g_x_0_y_xxxy, g_x_0_y_xxxz, g_x_0_y_xxyy, g_x_0_y_xxyz, g_x_0_y_xxzz, g_x_0_y_xyyy, g_x_0_y_xyyz, g_x_0_y_xyzz, g_x_0_y_xzzz, g_x_0_y_yyyy, g_x_0_y_yyyz, g_x_0_y_yyzz, g_x_0_y_yzzz, g_x_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_xxxx[k] = -g_x_0_0_xxxx[k] * cd_y[k] + g_x_0_0_xxxxy[k];

                g_x_0_y_xxxy[k] = -g_x_0_0_xxxy[k] * cd_y[k] + g_x_0_0_xxxyy[k];

                g_x_0_y_xxxz[k] = -g_x_0_0_xxxz[k] * cd_y[k] + g_x_0_0_xxxyz[k];

                g_x_0_y_xxyy[k] = -g_x_0_0_xxyy[k] * cd_y[k] + g_x_0_0_xxyyy[k];

                g_x_0_y_xxyz[k] = -g_x_0_0_xxyz[k] * cd_y[k] + g_x_0_0_xxyyz[k];

                g_x_0_y_xxzz[k] = -g_x_0_0_xxzz[k] * cd_y[k] + g_x_0_0_xxyzz[k];

                g_x_0_y_xyyy[k] = -g_x_0_0_xyyy[k] * cd_y[k] + g_x_0_0_xyyyy[k];

                g_x_0_y_xyyz[k] = -g_x_0_0_xyyz[k] * cd_y[k] + g_x_0_0_xyyyz[k];

                g_x_0_y_xyzz[k] = -g_x_0_0_xyzz[k] * cd_y[k] + g_x_0_0_xyyzz[k];

                g_x_0_y_xzzz[k] = -g_x_0_0_xzzz[k] * cd_y[k] + g_x_0_0_xyzzz[k];

                g_x_0_y_yyyy[k] = -g_x_0_0_yyyy[k] * cd_y[k] + g_x_0_0_yyyyy[k];

                g_x_0_y_yyyz[k] = -g_x_0_0_yyyz[k] * cd_y[k] + g_x_0_0_yyyyz[k];

                g_x_0_y_yyzz[k] = -g_x_0_0_yyzz[k] * cd_y[k] + g_x_0_0_yyyzz[k];

                g_x_0_y_yzzz[k] = -g_x_0_0_yzzz[k] * cd_y[k] + g_x_0_0_yyzzz[k];

                g_x_0_y_zzzz[k] = -g_x_0_0_zzzz[k] * cd_y[k] + g_x_0_0_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_xxxx = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_z_xxxy = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_z_xxxz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_z_xxyy = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_z_xxyz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_z_xxzz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_z_xyyy = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_z_xyyz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_z_xyzz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_z_xzzz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_z_yyyy = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_z_yyyz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_z_yyzz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_z_yzzz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_z_zzzz = cbuffer.data(pg_geom_10_off + 0 * acomps * bcomps + 44);

            #pragma omp simd aligned(cd_z, g_x_0_0_xxxx, g_x_0_0_xxxxz, g_x_0_0_xxxy, g_x_0_0_xxxyz, g_x_0_0_xxxz, g_x_0_0_xxxzz, g_x_0_0_xxyy, g_x_0_0_xxyyz, g_x_0_0_xxyz, g_x_0_0_xxyzz, g_x_0_0_xxzz, g_x_0_0_xxzzz, g_x_0_0_xyyy, g_x_0_0_xyyyz, g_x_0_0_xyyz, g_x_0_0_xyyzz, g_x_0_0_xyzz, g_x_0_0_xyzzz, g_x_0_0_xzzz, g_x_0_0_xzzzz, g_x_0_0_yyyy, g_x_0_0_yyyyz, g_x_0_0_yyyz, g_x_0_0_yyyzz, g_x_0_0_yyzz, g_x_0_0_yyzzz, g_x_0_0_yzzz, g_x_0_0_yzzzz, g_x_0_0_zzzz, g_x_0_0_zzzzz, g_x_0_z_xxxx, g_x_0_z_xxxy, g_x_0_z_xxxz, g_x_0_z_xxyy, g_x_0_z_xxyz, g_x_0_z_xxzz, g_x_0_z_xyyy, g_x_0_z_xyyz, g_x_0_z_xyzz, g_x_0_z_xzzz, g_x_0_z_yyyy, g_x_0_z_yyyz, g_x_0_z_yyzz, g_x_0_z_yzzz, g_x_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_xxxx[k] = -g_x_0_0_xxxx[k] * cd_z[k] + g_x_0_0_xxxxz[k];

                g_x_0_z_xxxy[k] = -g_x_0_0_xxxy[k] * cd_z[k] + g_x_0_0_xxxyz[k];

                g_x_0_z_xxxz[k] = -g_x_0_0_xxxz[k] * cd_z[k] + g_x_0_0_xxxzz[k];

                g_x_0_z_xxyy[k] = -g_x_0_0_xxyy[k] * cd_z[k] + g_x_0_0_xxyyz[k];

                g_x_0_z_xxyz[k] = -g_x_0_0_xxyz[k] * cd_z[k] + g_x_0_0_xxyzz[k];

                g_x_0_z_xxzz[k] = -g_x_0_0_xxzz[k] * cd_z[k] + g_x_0_0_xxzzz[k];

                g_x_0_z_xyyy[k] = -g_x_0_0_xyyy[k] * cd_z[k] + g_x_0_0_xyyyz[k];

                g_x_0_z_xyyz[k] = -g_x_0_0_xyyz[k] * cd_z[k] + g_x_0_0_xyyzz[k];

                g_x_0_z_xyzz[k] = -g_x_0_0_xyzz[k] * cd_z[k] + g_x_0_0_xyzzz[k];

                g_x_0_z_xzzz[k] = -g_x_0_0_xzzz[k] * cd_z[k] + g_x_0_0_xzzzz[k];

                g_x_0_z_yyyy[k] = -g_x_0_0_yyyy[k] * cd_z[k] + g_x_0_0_yyyyz[k];

                g_x_0_z_yyyz[k] = -g_x_0_0_yyyz[k] * cd_z[k] + g_x_0_0_yyyzz[k];

                g_x_0_z_yyzz[k] = -g_x_0_0_yyzz[k] * cd_z[k] + g_x_0_0_yyzzz[k];

                g_x_0_z_yzzz[k] = -g_x_0_0_yzzz[k] * cd_z[k] + g_x_0_0_yzzzz[k];

                g_x_0_z_zzzz[k] = -g_x_0_0_zzzz[k] * cd_z[k] + g_x_0_0_zzzzz[k];
            }
            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_xxxx = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 0);

            auto g_y_0_x_xxxy = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 1);

            auto g_y_0_x_xxxz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 2);

            auto g_y_0_x_xxyy = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 3);

            auto g_y_0_x_xxyz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 4);

            auto g_y_0_x_xxzz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 5);

            auto g_y_0_x_xyyy = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 6);

            auto g_y_0_x_xyyz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 7);

            auto g_y_0_x_xyzz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 8);

            auto g_y_0_x_xzzz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 9);

            auto g_y_0_x_yyyy = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 10);

            auto g_y_0_x_yyyz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 11);

            auto g_y_0_x_yyzz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 12);

            auto g_y_0_x_yzzz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 13);

            auto g_y_0_x_zzzz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_x, g_y_0_0_xxxx, g_y_0_0_xxxxx, g_y_0_0_xxxxy, g_y_0_0_xxxxz, g_y_0_0_xxxy, g_y_0_0_xxxyy, g_y_0_0_xxxyz, g_y_0_0_xxxz, g_y_0_0_xxxzz, g_y_0_0_xxyy, g_y_0_0_xxyyy, g_y_0_0_xxyyz, g_y_0_0_xxyz, g_y_0_0_xxyzz, g_y_0_0_xxzz, g_y_0_0_xxzzz, g_y_0_0_xyyy, g_y_0_0_xyyyy, g_y_0_0_xyyyz, g_y_0_0_xyyz, g_y_0_0_xyyzz, g_y_0_0_xyzz, g_y_0_0_xyzzz, g_y_0_0_xzzz, g_y_0_0_xzzzz, g_y_0_0_yyyy, g_y_0_0_yyyz, g_y_0_0_yyzz, g_y_0_0_yzzz, g_y_0_0_zzzz, g_y_0_x_xxxx, g_y_0_x_xxxy, g_y_0_x_xxxz, g_y_0_x_xxyy, g_y_0_x_xxyz, g_y_0_x_xxzz, g_y_0_x_xyyy, g_y_0_x_xyyz, g_y_0_x_xyzz, g_y_0_x_xzzz, g_y_0_x_yyyy, g_y_0_x_yyyz, g_y_0_x_yyzz, g_y_0_x_yzzz, g_y_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_xxxx[k] = -g_y_0_0_xxxx[k] * cd_x[k] + g_y_0_0_xxxxx[k];

                g_y_0_x_xxxy[k] = -g_y_0_0_xxxy[k] * cd_x[k] + g_y_0_0_xxxxy[k];

                g_y_0_x_xxxz[k] = -g_y_0_0_xxxz[k] * cd_x[k] + g_y_0_0_xxxxz[k];

                g_y_0_x_xxyy[k] = -g_y_0_0_xxyy[k] * cd_x[k] + g_y_0_0_xxxyy[k];

                g_y_0_x_xxyz[k] = -g_y_0_0_xxyz[k] * cd_x[k] + g_y_0_0_xxxyz[k];

                g_y_0_x_xxzz[k] = -g_y_0_0_xxzz[k] * cd_x[k] + g_y_0_0_xxxzz[k];

                g_y_0_x_xyyy[k] = -g_y_0_0_xyyy[k] * cd_x[k] + g_y_0_0_xxyyy[k];

                g_y_0_x_xyyz[k] = -g_y_0_0_xyyz[k] * cd_x[k] + g_y_0_0_xxyyz[k];

                g_y_0_x_xyzz[k] = -g_y_0_0_xyzz[k] * cd_x[k] + g_y_0_0_xxyzz[k];

                g_y_0_x_xzzz[k] = -g_y_0_0_xzzz[k] * cd_x[k] + g_y_0_0_xxzzz[k];

                g_y_0_x_yyyy[k] = -g_y_0_0_yyyy[k] * cd_x[k] + g_y_0_0_xyyyy[k];

                g_y_0_x_yyyz[k] = -g_y_0_0_yyyz[k] * cd_x[k] + g_y_0_0_xyyyz[k];

                g_y_0_x_yyzz[k] = -g_y_0_0_yyzz[k] * cd_x[k] + g_y_0_0_xyyzz[k];

                g_y_0_x_yzzz[k] = -g_y_0_0_yzzz[k] * cd_x[k] + g_y_0_0_xyzzz[k];

                g_y_0_x_zzzz[k] = -g_y_0_0_zzzz[k] * cd_x[k] + g_y_0_0_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_xxxx = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 15);

            auto g_y_0_y_xxxy = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 16);

            auto g_y_0_y_xxxz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 17);

            auto g_y_0_y_xxyy = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 18);

            auto g_y_0_y_xxyz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 19);

            auto g_y_0_y_xxzz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 20);

            auto g_y_0_y_xyyy = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 21);

            auto g_y_0_y_xyyz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 22);

            auto g_y_0_y_xyzz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 23);

            auto g_y_0_y_xzzz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 24);

            auto g_y_0_y_yyyy = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 25);

            auto g_y_0_y_yyyz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 26);

            auto g_y_0_y_yyzz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 27);

            auto g_y_0_y_yzzz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 28);

            auto g_y_0_y_zzzz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_y, g_0_xxxx, g_0_xxxy, g_0_xxxz, g_0_xxyy, g_0_xxyz, g_0_xxzz, g_0_xyyy, g_0_xyyz, g_0_xyzz, g_0_xzzz, g_0_yyyy, g_0_yyyz, g_0_yyzz, g_0_yzzz, g_0_zzzz, g_y_0_0_xxxx, g_y_0_0_xxxxy, g_y_0_0_xxxy, g_y_0_0_xxxyy, g_y_0_0_xxxyz, g_y_0_0_xxxz, g_y_0_0_xxyy, g_y_0_0_xxyyy, g_y_0_0_xxyyz, g_y_0_0_xxyz, g_y_0_0_xxyzz, g_y_0_0_xxzz, g_y_0_0_xyyy, g_y_0_0_xyyyy, g_y_0_0_xyyyz, g_y_0_0_xyyz, g_y_0_0_xyyzz, g_y_0_0_xyzz, g_y_0_0_xyzzz, g_y_0_0_xzzz, g_y_0_0_yyyy, g_y_0_0_yyyyy, g_y_0_0_yyyyz, g_y_0_0_yyyz, g_y_0_0_yyyzz, g_y_0_0_yyzz, g_y_0_0_yyzzz, g_y_0_0_yzzz, g_y_0_0_yzzzz, g_y_0_0_zzzz, g_y_0_y_xxxx, g_y_0_y_xxxy, g_y_0_y_xxxz, g_y_0_y_xxyy, g_y_0_y_xxyz, g_y_0_y_xxzz, g_y_0_y_xyyy, g_y_0_y_xyyz, g_y_0_y_xyzz, g_y_0_y_xzzz, g_y_0_y_yyyy, g_y_0_y_yyyz, g_y_0_y_yyzz, g_y_0_y_yzzz, g_y_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_xxxx[k] = -g_0_xxxx[k] - g_y_0_0_xxxx[k] * cd_y[k] + g_y_0_0_xxxxy[k];

                g_y_0_y_xxxy[k] = -g_0_xxxy[k] - g_y_0_0_xxxy[k] * cd_y[k] + g_y_0_0_xxxyy[k];

                g_y_0_y_xxxz[k] = -g_0_xxxz[k] - g_y_0_0_xxxz[k] * cd_y[k] + g_y_0_0_xxxyz[k];

                g_y_0_y_xxyy[k] = -g_0_xxyy[k] - g_y_0_0_xxyy[k] * cd_y[k] + g_y_0_0_xxyyy[k];

                g_y_0_y_xxyz[k] = -g_0_xxyz[k] - g_y_0_0_xxyz[k] * cd_y[k] + g_y_0_0_xxyyz[k];

                g_y_0_y_xxzz[k] = -g_0_xxzz[k] - g_y_0_0_xxzz[k] * cd_y[k] + g_y_0_0_xxyzz[k];

                g_y_0_y_xyyy[k] = -g_0_xyyy[k] - g_y_0_0_xyyy[k] * cd_y[k] + g_y_0_0_xyyyy[k];

                g_y_0_y_xyyz[k] = -g_0_xyyz[k] - g_y_0_0_xyyz[k] * cd_y[k] + g_y_0_0_xyyyz[k];

                g_y_0_y_xyzz[k] = -g_0_xyzz[k] - g_y_0_0_xyzz[k] * cd_y[k] + g_y_0_0_xyyzz[k];

                g_y_0_y_xzzz[k] = -g_0_xzzz[k] - g_y_0_0_xzzz[k] * cd_y[k] + g_y_0_0_xyzzz[k];

                g_y_0_y_yyyy[k] = -g_0_yyyy[k] - g_y_0_0_yyyy[k] * cd_y[k] + g_y_0_0_yyyyy[k];

                g_y_0_y_yyyz[k] = -g_0_yyyz[k] - g_y_0_0_yyyz[k] * cd_y[k] + g_y_0_0_yyyyz[k];

                g_y_0_y_yyzz[k] = -g_0_yyzz[k] - g_y_0_0_yyzz[k] * cd_y[k] + g_y_0_0_yyyzz[k];

                g_y_0_y_yzzz[k] = -g_0_yzzz[k] - g_y_0_0_yzzz[k] * cd_y[k] + g_y_0_0_yyzzz[k];

                g_y_0_y_zzzz[k] = -g_0_zzzz[k] - g_y_0_0_zzzz[k] * cd_y[k] + g_y_0_0_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_xxxx = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 30);

            auto g_y_0_z_xxxy = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 31);

            auto g_y_0_z_xxxz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 32);

            auto g_y_0_z_xxyy = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 33);

            auto g_y_0_z_xxyz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 34);

            auto g_y_0_z_xxzz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 35);

            auto g_y_0_z_xyyy = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 36);

            auto g_y_0_z_xyyz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 37);

            auto g_y_0_z_xyzz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 38);

            auto g_y_0_z_xzzz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 39);

            auto g_y_0_z_yyyy = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 40);

            auto g_y_0_z_yyyz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 41);

            auto g_y_0_z_yyzz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 42);

            auto g_y_0_z_yzzz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 43);

            auto g_y_0_z_zzzz = cbuffer.data(pg_geom_10_off + 45 * acomps * bcomps + 44);

            #pragma omp simd aligned(cd_z, g_y_0_0_xxxx, g_y_0_0_xxxxz, g_y_0_0_xxxy, g_y_0_0_xxxyz, g_y_0_0_xxxz, g_y_0_0_xxxzz, g_y_0_0_xxyy, g_y_0_0_xxyyz, g_y_0_0_xxyz, g_y_0_0_xxyzz, g_y_0_0_xxzz, g_y_0_0_xxzzz, g_y_0_0_xyyy, g_y_0_0_xyyyz, g_y_0_0_xyyz, g_y_0_0_xyyzz, g_y_0_0_xyzz, g_y_0_0_xyzzz, g_y_0_0_xzzz, g_y_0_0_xzzzz, g_y_0_0_yyyy, g_y_0_0_yyyyz, g_y_0_0_yyyz, g_y_0_0_yyyzz, g_y_0_0_yyzz, g_y_0_0_yyzzz, g_y_0_0_yzzz, g_y_0_0_yzzzz, g_y_0_0_zzzz, g_y_0_0_zzzzz, g_y_0_z_xxxx, g_y_0_z_xxxy, g_y_0_z_xxxz, g_y_0_z_xxyy, g_y_0_z_xxyz, g_y_0_z_xxzz, g_y_0_z_xyyy, g_y_0_z_xyyz, g_y_0_z_xyzz, g_y_0_z_xzzz, g_y_0_z_yyyy, g_y_0_z_yyyz, g_y_0_z_yyzz, g_y_0_z_yzzz, g_y_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_xxxx[k] = -g_y_0_0_xxxx[k] * cd_z[k] + g_y_0_0_xxxxz[k];

                g_y_0_z_xxxy[k] = -g_y_0_0_xxxy[k] * cd_z[k] + g_y_0_0_xxxyz[k];

                g_y_0_z_xxxz[k] = -g_y_0_0_xxxz[k] * cd_z[k] + g_y_0_0_xxxzz[k];

                g_y_0_z_xxyy[k] = -g_y_0_0_xxyy[k] * cd_z[k] + g_y_0_0_xxyyz[k];

                g_y_0_z_xxyz[k] = -g_y_0_0_xxyz[k] * cd_z[k] + g_y_0_0_xxyzz[k];

                g_y_0_z_xxzz[k] = -g_y_0_0_xxzz[k] * cd_z[k] + g_y_0_0_xxzzz[k];

                g_y_0_z_xyyy[k] = -g_y_0_0_xyyy[k] * cd_z[k] + g_y_0_0_xyyyz[k];

                g_y_0_z_xyyz[k] = -g_y_0_0_xyyz[k] * cd_z[k] + g_y_0_0_xyyzz[k];

                g_y_0_z_xyzz[k] = -g_y_0_0_xyzz[k] * cd_z[k] + g_y_0_0_xyzzz[k];

                g_y_0_z_xzzz[k] = -g_y_0_0_xzzz[k] * cd_z[k] + g_y_0_0_xzzzz[k];

                g_y_0_z_yyyy[k] = -g_y_0_0_yyyy[k] * cd_z[k] + g_y_0_0_yyyyz[k];

                g_y_0_z_yyyz[k] = -g_y_0_0_yyyz[k] * cd_z[k] + g_y_0_0_yyyzz[k];

                g_y_0_z_yyzz[k] = -g_y_0_0_yyzz[k] * cd_z[k] + g_y_0_0_yyzzz[k];

                g_y_0_z_yzzz[k] = -g_y_0_0_yzzz[k] * cd_z[k] + g_y_0_0_yzzzz[k];

                g_y_0_z_zzzz[k] = -g_y_0_0_zzzz[k] * cd_z[k] + g_y_0_0_zzzzz[k];
            }
            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_xxxx = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 0);

            auto g_z_0_x_xxxy = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 1);

            auto g_z_0_x_xxxz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 2);

            auto g_z_0_x_xxyy = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 3);

            auto g_z_0_x_xxyz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 4);

            auto g_z_0_x_xxzz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 5);

            auto g_z_0_x_xyyy = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 6);

            auto g_z_0_x_xyyz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 7);

            auto g_z_0_x_xyzz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 8);

            auto g_z_0_x_xzzz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 9);

            auto g_z_0_x_yyyy = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 10);

            auto g_z_0_x_yyyz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 11);

            auto g_z_0_x_yyzz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 12);

            auto g_z_0_x_yzzz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 13);

            auto g_z_0_x_zzzz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_x, g_z_0_0_xxxx, g_z_0_0_xxxxx, g_z_0_0_xxxxy, g_z_0_0_xxxxz, g_z_0_0_xxxy, g_z_0_0_xxxyy, g_z_0_0_xxxyz, g_z_0_0_xxxz, g_z_0_0_xxxzz, g_z_0_0_xxyy, g_z_0_0_xxyyy, g_z_0_0_xxyyz, g_z_0_0_xxyz, g_z_0_0_xxyzz, g_z_0_0_xxzz, g_z_0_0_xxzzz, g_z_0_0_xyyy, g_z_0_0_xyyyy, g_z_0_0_xyyyz, g_z_0_0_xyyz, g_z_0_0_xyyzz, g_z_0_0_xyzz, g_z_0_0_xyzzz, g_z_0_0_xzzz, g_z_0_0_xzzzz, g_z_0_0_yyyy, g_z_0_0_yyyz, g_z_0_0_yyzz, g_z_0_0_yzzz, g_z_0_0_zzzz, g_z_0_x_xxxx, g_z_0_x_xxxy, g_z_0_x_xxxz, g_z_0_x_xxyy, g_z_0_x_xxyz, g_z_0_x_xxzz, g_z_0_x_xyyy, g_z_0_x_xyyz, g_z_0_x_xyzz, g_z_0_x_xzzz, g_z_0_x_yyyy, g_z_0_x_yyyz, g_z_0_x_yyzz, g_z_0_x_yzzz, g_z_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_xxxx[k] = -g_z_0_0_xxxx[k] * cd_x[k] + g_z_0_0_xxxxx[k];

                g_z_0_x_xxxy[k] = -g_z_0_0_xxxy[k] * cd_x[k] + g_z_0_0_xxxxy[k];

                g_z_0_x_xxxz[k] = -g_z_0_0_xxxz[k] * cd_x[k] + g_z_0_0_xxxxz[k];

                g_z_0_x_xxyy[k] = -g_z_0_0_xxyy[k] * cd_x[k] + g_z_0_0_xxxyy[k];

                g_z_0_x_xxyz[k] = -g_z_0_0_xxyz[k] * cd_x[k] + g_z_0_0_xxxyz[k];

                g_z_0_x_xxzz[k] = -g_z_0_0_xxzz[k] * cd_x[k] + g_z_0_0_xxxzz[k];

                g_z_0_x_xyyy[k] = -g_z_0_0_xyyy[k] * cd_x[k] + g_z_0_0_xxyyy[k];

                g_z_0_x_xyyz[k] = -g_z_0_0_xyyz[k] * cd_x[k] + g_z_0_0_xxyyz[k];

                g_z_0_x_xyzz[k] = -g_z_0_0_xyzz[k] * cd_x[k] + g_z_0_0_xxyzz[k];

                g_z_0_x_xzzz[k] = -g_z_0_0_xzzz[k] * cd_x[k] + g_z_0_0_xxzzz[k];

                g_z_0_x_yyyy[k] = -g_z_0_0_yyyy[k] * cd_x[k] + g_z_0_0_xyyyy[k];

                g_z_0_x_yyyz[k] = -g_z_0_0_yyyz[k] * cd_x[k] + g_z_0_0_xyyyz[k];

                g_z_0_x_yyzz[k] = -g_z_0_0_yyzz[k] * cd_x[k] + g_z_0_0_xyyzz[k];

                g_z_0_x_yzzz[k] = -g_z_0_0_yzzz[k] * cd_x[k] + g_z_0_0_xyzzz[k];

                g_z_0_x_zzzz[k] = -g_z_0_0_zzzz[k] * cd_x[k] + g_z_0_0_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_xxxx = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 15);

            auto g_z_0_y_xxxy = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 16);

            auto g_z_0_y_xxxz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 17);

            auto g_z_0_y_xxyy = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 18);

            auto g_z_0_y_xxyz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 19);

            auto g_z_0_y_xxzz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 20);

            auto g_z_0_y_xyyy = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 21);

            auto g_z_0_y_xyyz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 22);

            auto g_z_0_y_xyzz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 23);

            auto g_z_0_y_xzzz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 24);

            auto g_z_0_y_yyyy = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 25);

            auto g_z_0_y_yyyz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 26);

            auto g_z_0_y_yyzz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 27);

            auto g_z_0_y_yzzz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 28);

            auto g_z_0_y_zzzz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_y, g_z_0_0_xxxx, g_z_0_0_xxxxy, g_z_0_0_xxxy, g_z_0_0_xxxyy, g_z_0_0_xxxyz, g_z_0_0_xxxz, g_z_0_0_xxyy, g_z_0_0_xxyyy, g_z_0_0_xxyyz, g_z_0_0_xxyz, g_z_0_0_xxyzz, g_z_0_0_xxzz, g_z_0_0_xyyy, g_z_0_0_xyyyy, g_z_0_0_xyyyz, g_z_0_0_xyyz, g_z_0_0_xyyzz, g_z_0_0_xyzz, g_z_0_0_xyzzz, g_z_0_0_xzzz, g_z_0_0_yyyy, g_z_0_0_yyyyy, g_z_0_0_yyyyz, g_z_0_0_yyyz, g_z_0_0_yyyzz, g_z_0_0_yyzz, g_z_0_0_yyzzz, g_z_0_0_yzzz, g_z_0_0_yzzzz, g_z_0_0_zzzz, g_z_0_y_xxxx, g_z_0_y_xxxy, g_z_0_y_xxxz, g_z_0_y_xxyy, g_z_0_y_xxyz, g_z_0_y_xxzz, g_z_0_y_xyyy, g_z_0_y_xyyz, g_z_0_y_xyzz, g_z_0_y_xzzz, g_z_0_y_yyyy, g_z_0_y_yyyz, g_z_0_y_yyzz, g_z_0_y_yzzz, g_z_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_xxxx[k] = -g_z_0_0_xxxx[k] * cd_y[k] + g_z_0_0_xxxxy[k];

                g_z_0_y_xxxy[k] = -g_z_0_0_xxxy[k] * cd_y[k] + g_z_0_0_xxxyy[k];

                g_z_0_y_xxxz[k] = -g_z_0_0_xxxz[k] * cd_y[k] + g_z_0_0_xxxyz[k];

                g_z_0_y_xxyy[k] = -g_z_0_0_xxyy[k] * cd_y[k] + g_z_0_0_xxyyy[k];

                g_z_0_y_xxyz[k] = -g_z_0_0_xxyz[k] * cd_y[k] + g_z_0_0_xxyyz[k];

                g_z_0_y_xxzz[k] = -g_z_0_0_xxzz[k] * cd_y[k] + g_z_0_0_xxyzz[k];

                g_z_0_y_xyyy[k] = -g_z_0_0_xyyy[k] * cd_y[k] + g_z_0_0_xyyyy[k];

                g_z_0_y_xyyz[k] = -g_z_0_0_xyyz[k] * cd_y[k] + g_z_0_0_xyyyz[k];

                g_z_0_y_xyzz[k] = -g_z_0_0_xyzz[k] * cd_y[k] + g_z_0_0_xyyzz[k];

                g_z_0_y_xzzz[k] = -g_z_0_0_xzzz[k] * cd_y[k] + g_z_0_0_xyzzz[k];

                g_z_0_y_yyyy[k] = -g_z_0_0_yyyy[k] * cd_y[k] + g_z_0_0_yyyyy[k];

                g_z_0_y_yyyz[k] = -g_z_0_0_yyyz[k] * cd_y[k] + g_z_0_0_yyyyz[k];

                g_z_0_y_yyzz[k] = -g_z_0_0_yyzz[k] * cd_y[k] + g_z_0_0_yyyzz[k];

                g_z_0_y_yzzz[k] = -g_z_0_0_yzzz[k] * cd_y[k] + g_z_0_0_yyzzz[k];

                g_z_0_y_zzzz[k] = -g_z_0_0_zzzz[k] * cd_y[k] + g_z_0_0_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_xxxx = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 30);

            auto g_z_0_z_xxxy = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 31);

            auto g_z_0_z_xxxz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 32);

            auto g_z_0_z_xxyy = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 33);

            auto g_z_0_z_xxyz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 34);

            auto g_z_0_z_xxzz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 35);

            auto g_z_0_z_xyyy = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 36);

            auto g_z_0_z_xyyz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 37);

            auto g_z_0_z_xyzz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 38);

            auto g_z_0_z_xzzz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 39);

            auto g_z_0_z_yyyy = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 40);

            auto g_z_0_z_yyyz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 41);

            auto g_z_0_z_yyzz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 42);

            auto g_z_0_z_yzzz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 43);

            auto g_z_0_z_zzzz = cbuffer.data(pg_geom_10_off + 90 * acomps * bcomps + 44);

            #pragma omp simd aligned(cd_z, g_0_xxxx, g_0_xxxy, g_0_xxxz, g_0_xxyy, g_0_xxyz, g_0_xxzz, g_0_xyyy, g_0_xyyz, g_0_xyzz, g_0_xzzz, g_0_yyyy, g_0_yyyz, g_0_yyzz, g_0_yzzz, g_0_zzzz, g_z_0_0_xxxx, g_z_0_0_xxxxz, g_z_0_0_xxxy, g_z_0_0_xxxyz, g_z_0_0_xxxz, g_z_0_0_xxxzz, g_z_0_0_xxyy, g_z_0_0_xxyyz, g_z_0_0_xxyz, g_z_0_0_xxyzz, g_z_0_0_xxzz, g_z_0_0_xxzzz, g_z_0_0_xyyy, g_z_0_0_xyyyz, g_z_0_0_xyyz, g_z_0_0_xyyzz, g_z_0_0_xyzz, g_z_0_0_xyzzz, g_z_0_0_xzzz, g_z_0_0_xzzzz, g_z_0_0_yyyy, g_z_0_0_yyyyz, g_z_0_0_yyyz, g_z_0_0_yyyzz, g_z_0_0_yyzz, g_z_0_0_yyzzz, g_z_0_0_yzzz, g_z_0_0_yzzzz, g_z_0_0_zzzz, g_z_0_0_zzzzz, g_z_0_z_xxxx, g_z_0_z_xxxy, g_z_0_z_xxxz, g_z_0_z_xxyy, g_z_0_z_xxyz, g_z_0_z_xxzz, g_z_0_z_xyyy, g_z_0_z_xyyz, g_z_0_z_xyzz, g_z_0_z_xzzz, g_z_0_z_yyyy, g_z_0_z_yyyz, g_z_0_z_yyzz, g_z_0_z_yzzz, g_z_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_xxxx[k] = -g_0_xxxx[k] - g_z_0_0_xxxx[k] * cd_z[k] + g_z_0_0_xxxxz[k];

                g_z_0_z_xxxy[k] = -g_0_xxxy[k] - g_z_0_0_xxxy[k] * cd_z[k] + g_z_0_0_xxxyz[k];

                g_z_0_z_xxxz[k] = -g_0_xxxz[k] - g_z_0_0_xxxz[k] * cd_z[k] + g_z_0_0_xxxzz[k];

                g_z_0_z_xxyy[k] = -g_0_xxyy[k] - g_z_0_0_xxyy[k] * cd_z[k] + g_z_0_0_xxyyz[k];

                g_z_0_z_xxyz[k] = -g_0_xxyz[k] - g_z_0_0_xxyz[k] * cd_z[k] + g_z_0_0_xxyzz[k];

                g_z_0_z_xxzz[k] = -g_0_xxzz[k] - g_z_0_0_xxzz[k] * cd_z[k] + g_z_0_0_xxzzz[k];

                g_z_0_z_xyyy[k] = -g_0_xyyy[k] - g_z_0_0_xyyy[k] * cd_z[k] + g_z_0_0_xyyyz[k];

                g_z_0_z_xyyz[k] = -g_0_xyyz[k] - g_z_0_0_xyyz[k] * cd_z[k] + g_z_0_0_xyyzz[k];

                g_z_0_z_xyzz[k] = -g_0_xyzz[k] - g_z_0_0_xyzz[k] * cd_z[k] + g_z_0_0_xyzzz[k];

                g_z_0_z_xzzz[k] = -g_0_xzzz[k] - g_z_0_0_xzzz[k] * cd_z[k] + g_z_0_0_xzzzz[k];

                g_z_0_z_yyyy[k] = -g_0_yyyy[k] - g_z_0_0_yyyy[k] * cd_z[k] + g_z_0_0_yyyyz[k];

                g_z_0_z_yyyz[k] = -g_0_yyyz[k] - g_z_0_0_yyyz[k] * cd_z[k] + g_z_0_0_yyyzz[k];

                g_z_0_z_yyzz[k] = -g_0_yyzz[k] - g_z_0_0_yyzz[k] * cd_z[k] + g_z_0_0_yyzzz[k];

                g_z_0_z_yzzz[k] = -g_0_yzzz[k] - g_z_0_0_yzzz[k] * cd_z[k] + g_z_0_0_yzzzz[k];

                g_z_0_z_zzzz[k] = -g_0_zzzz[k] - g_z_0_0_zzzz[k] * cd_z[k] + g_z_0_0_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

