#include "ElectronRepulsionContrRecXXDG.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_hrr_electron_repulsion_xxdg(CSimdArray<double>& cbuffer,
                                     const size_t idx_xxdg,
                                     const size_t idx_xxpg,
                                     const size_t idx_xxph,
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
            /// Set up components of auxilary buffer : SSPG

            const auto pg_off = idx_xxpg + (i * bcomps + j) * 45;

            auto g_x_xxxx = cbuffer.data(pg_off + 0);

            auto g_x_xxxy = cbuffer.data(pg_off + 1);

            auto g_x_xxxz = cbuffer.data(pg_off + 2);

            auto g_x_xxyy = cbuffer.data(pg_off + 3);

            auto g_x_xxyz = cbuffer.data(pg_off + 4);

            auto g_x_xxzz = cbuffer.data(pg_off + 5);

            auto g_x_xyyy = cbuffer.data(pg_off + 6);

            auto g_x_xyyz = cbuffer.data(pg_off + 7);

            auto g_x_xyzz = cbuffer.data(pg_off + 8);

            auto g_x_xzzz = cbuffer.data(pg_off + 9);

            auto g_x_yyyy = cbuffer.data(pg_off + 10);

            auto g_x_yyyz = cbuffer.data(pg_off + 11);

            auto g_x_yyzz = cbuffer.data(pg_off + 12);

            auto g_x_yzzz = cbuffer.data(pg_off + 13);

            auto g_x_zzzz = cbuffer.data(pg_off + 14);

            auto g_y_xxxx = cbuffer.data(pg_off + 15);

            auto g_y_xxxy = cbuffer.data(pg_off + 16);

            auto g_y_xxxz = cbuffer.data(pg_off + 17);

            auto g_y_xxyy = cbuffer.data(pg_off + 18);

            auto g_y_xxyz = cbuffer.data(pg_off + 19);

            auto g_y_xxzz = cbuffer.data(pg_off + 20);

            auto g_y_xyyy = cbuffer.data(pg_off + 21);

            auto g_y_xyyz = cbuffer.data(pg_off + 22);

            auto g_y_xyzz = cbuffer.data(pg_off + 23);

            auto g_y_xzzz = cbuffer.data(pg_off + 24);

            auto g_y_yyyy = cbuffer.data(pg_off + 25);

            auto g_y_yyyz = cbuffer.data(pg_off + 26);

            auto g_y_yyzz = cbuffer.data(pg_off + 27);

            auto g_y_yzzz = cbuffer.data(pg_off + 28);

            auto g_y_zzzz = cbuffer.data(pg_off + 29);

            auto g_z_xxxx = cbuffer.data(pg_off + 30);

            auto g_z_xxxy = cbuffer.data(pg_off + 31);

            auto g_z_xxxz = cbuffer.data(pg_off + 32);

            auto g_z_xxyy = cbuffer.data(pg_off + 33);

            auto g_z_xxyz = cbuffer.data(pg_off + 34);

            auto g_z_xxzz = cbuffer.data(pg_off + 35);

            auto g_z_xyyy = cbuffer.data(pg_off + 36);

            auto g_z_xyyz = cbuffer.data(pg_off + 37);

            auto g_z_xyzz = cbuffer.data(pg_off + 38);

            auto g_z_xzzz = cbuffer.data(pg_off + 39);

            auto g_z_yyyy = cbuffer.data(pg_off + 40);

            auto g_z_yyyz = cbuffer.data(pg_off + 41);

            auto g_z_yyzz = cbuffer.data(pg_off + 42);

            auto g_z_yzzz = cbuffer.data(pg_off + 43);

            auto g_z_zzzz = cbuffer.data(pg_off + 44);

            /// Set up components of auxilary buffer : SSPH

            const auto ph_off = idx_xxph + (i * bcomps + j) * 63;

            auto g_x_xxxxx = cbuffer.data(ph_off + 0);

            auto g_x_xxxxy = cbuffer.data(ph_off + 1);

            auto g_x_xxxxz = cbuffer.data(ph_off + 2);

            auto g_x_xxxyy = cbuffer.data(ph_off + 3);

            auto g_x_xxxyz = cbuffer.data(ph_off + 4);

            auto g_x_xxxzz = cbuffer.data(ph_off + 5);

            auto g_x_xxyyy = cbuffer.data(ph_off + 6);

            auto g_x_xxyyz = cbuffer.data(ph_off + 7);

            auto g_x_xxyzz = cbuffer.data(ph_off + 8);

            auto g_x_xxzzz = cbuffer.data(ph_off + 9);

            auto g_x_xyyyy = cbuffer.data(ph_off + 10);

            auto g_x_xyyyz = cbuffer.data(ph_off + 11);

            auto g_x_xyyzz = cbuffer.data(ph_off + 12);

            auto g_x_xyzzz = cbuffer.data(ph_off + 13);

            auto g_x_xzzzz = cbuffer.data(ph_off + 14);

            auto g_y_xxxxx = cbuffer.data(ph_off + 21);

            auto g_y_xxxxy = cbuffer.data(ph_off + 22);

            auto g_y_xxxxz = cbuffer.data(ph_off + 23);

            auto g_y_xxxyy = cbuffer.data(ph_off + 24);

            auto g_y_xxxyz = cbuffer.data(ph_off + 25);

            auto g_y_xxxzz = cbuffer.data(ph_off + 26);

            auto g_y_xxyyy = cbuffer.data(ph_off + 27);

            auto g_y_xxyyz = cbuffer.data(ph_off + 28);

            auto g_y_xxyzz = cbuffer.data(ph_off + 29);

            auto g_y_xxzzz = cbuffer.data(ph_off + 30);

            auto g_y_xyyyy = cbuffer.data(ph_off + 31);

            auto g_y_xyyyz = cbuffer.data(ph_off + 32);

            auto g_y_xyyzz = cbuffer.data(ph_off + 33);

            auto g_y_xyzzz = cbuffer.data(ph_off + 34);

            auto g_y_xzzzz = cbuffer.data(ph_off + 35);

            auto g_y_yyyyy = cbuffer.data(ph_off + 36);

            auto g_y_yyyyz = cbuffer.data(ph_off + 37);

            auto g_y_yyyzz = cbuffer.data(ph_off + 38);

            auto g_y_yyzzz = cbuffer.data(ph_off + 39);

            auto g_y_yzzzz = cbuffer.data(ph_off + 40);

            auto g_z_xxxxx = cbuffer.data(ph_off + 42);

            auto g_z_xxxxy = cbuffer.data(ph_off + 43);

            auto g_z_xxxxz = cbuffer.data(ph_off + 44);

            auto g_z_xxxyy = cbuffer.data(ph_off + 45);

            auto g_z_xxxyz = cbuffer.data(ph_off + 46);

            auto g_z_xxxzz = cbuffer.data(ph_off + 47);

            auto g_z_xxyyy = cbuffer.data(ph_off + 48);

            auto g_z_xxyyz = cbuffer.data(ph_off + 49);

            auto g_z_xxyzz = cbuffer.data(ph_off + 50);

            auto g_z_xxzzz = cbuffer.data(ph_off + 51);

            auto g_z_xyyyy = cbuffer.data(ph_off + 52);

            auto g_z_xyyyz = cbuffer.data(ph_off + 53);

            auto g_z_xyyzz = cbuffer.data(ph_off + 54);

            auto g_z_xyzzz = cbuffer.data(ph_off + 55);

            auto g_z_xzzzz = cbuffer.data(ph_off + 56);

            auto g_z_yyyyy = cbuffer.data(ph_off + 57);

            auto g_z_yyyyz = cbuffer.data(ph_off + 58);

            auto g_z_yyyzz = cbuffer.data(ph_off + 59);

            auto g_z_yyzzz = cbuffer.data(ph_off + 60);

            auto g_z_yzzzz = cbuffer.data(ph_off + 61);

            auto g_z_zzzzz = cbuffer.data(ph_off + 62);

            /// set up bra offset for contr_buffer_xxdg

            const auto dg_off = idx_xxdg + (i * bcomps + j) * 90;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_xx_xxxx = cbuffer.data(dg_off + 0);

            auto g_xx_xxxy = cbuffer.data(dg_off + 1);

            auto g_xx_xxxz = cbuffer.data(dg_off + 2);

            auto g_xx_xxyy = cbuffer.data(dg_off + 3);

            auto g_xx_xxyz = cbuffer.data(dg_off + 4);

            auto g_xx_xxzz = cbuffer.data(dg_off + 5);

            auto g_xx_xyyy = cbuffer.data(dg_off + 6);

            auto g_xx_xyyz = cbuffer.data(dg_off + 7);

            auto g_xx_xyzz = cbuffer.data(dg_off + 8);

            auto g_xx_xzzz = cbuffer.data(dg_off + 9);

            auto g_xx_yyyy = cbuffer.data(dg_off + 10);

            auto g_xx_yyyz = cbuffer.data(dg_off + 11);

            auto g_xx_yyzz = cbuffer.data(dg_off + 12);

            auto g_xx_yzzz = cbuffer.data(dg_off + 13);

            auto g_xx_zzzz = cbuffer.data(dg_off + 14);

            #pragma omp simd aligned(cd_x, g_x_xxxx, g_x_xxxxx, g_x_xxxxy, g_x_xxxxz, g_x_xxxy, g_x_xxxyy, g_x_xxxyz, g_x_xxxz, g_x_xxxzz, g_x_xxyy, g_x_xxyyy, g_x_xxyyz, g_x_xxyz, g_x_xxyzz, g_x_xxzz, g_x_xxzzz, g_x_xyyy, g_x_xyyyy, g_x_xyyyz, g_x_xyyz, g_x_xyyzz, g_x_xyzz, g_x_xyzzz, g_x_xzzz, g_x_xzzzz, g_x_yyyy, g_x_yyyz, g_x_yyzz, g_x_yzzz, g_x_zzzz, g_xx_xxxx, g_xx_xxxy, g_xx_xxxz, g_xx_xxyy, g_xx_xxyz, g_xx_xxzz, g_xx_xyyy, g_xx_xyyz, g_xx_xyzz, g_xx_xzzz, g_xx_yyyy, g_xx_yyyz, g_xx_yyzz, g_xx_yzzz, g_xx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_xxxx[k] = -g_x_xxxx[k] * cd_x[k] + g_x_xxxxx[k];

                g_xx_xxxy[k] = -g_x_xxxy[k] * cd_x[k] + g_x_xxxxy[k];

                g_xx_xxxz[k] = -g_x_xxxz[k] * cd_x[k] + g_x_xxxxz[k];

                g_xx_xxyy[k] = -g_x_xxyy[k] * cd_x[k] + g_x_xxxyy[k];

                g_xx_xxyz[k] = -g_x_xxyz[k] * cd_x[k] + g_x_xxxyz[k];

                g_xx_xxzz[k] = -g_x_xxzz[k] * cd_x[k] + g_x_xxxzz[k];

                g_xx_xyyy[k] = -g_x_xyyy[k] * cd_x[k] + g_x_xxyyy[k];

                g_xx_xyyz[k] = -g_x_xyyz[k] * cd_x[k] + g_x_xxyyz[k];

                g_xx_xyzz[k] = -g_x_xyzz[k] * cd_x[k] + g_x_xxyzz[k];

                g_xx_xzzz[k] = -g_x_xzzz[k] * cd_x[k] + g_x_xxzzz[k];

                g_xx_yyyy[k] = -g_x_yyyy[k] * cd_x[k] + g_x_xyyyy[k];

                g_xx_yyyz[k] = -g_x_yyyz[k] * cd_x[k] + g_x_xyyyz[k];

                g_xx_yyzz[k] = -g_x_yyzz[k] * cd_x[k] + g_x_xyyzz[k];

                g_xx_yzzz[k] = -g_x_yzzz[k] * cd_x[k] + g_x_xyzzz[k];

                g_xx_zzzz[k] = -g_x_zzzz[k] * cd_x[k] + g_x_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_xy_xxxx = cbuffer.data(dg_off + 15);

            auto g_xy_xxxy = cbuffer.data(dg_off + 16);

            auto g_xy_xxxz = cbuffer.data(dg_off + 17);

            auto g_xy_xxyy = cbuffer.data(dg_off + 18);

            auto g_xy_xxyz = cbuffer.data(dg_off + 19);

            auto g_xy_xxzz = cbuffer.data(dg_off + 20);

            auto g_xy_xyyy = cbuffer.data(dg_off + 21);

            auto g_xy_xyyz = cbuffer.data(dg_off + 22);

            auto g_xy_xyzz = cbuffer.data(dg_off + 23);

            auto g_xy_xzzz = cbuffer.data(dg_off + 24);

            auto g_xy_yyyy = cbuffer.data(dg_off + 25);

            auto g_xy_yyyz = cbuffer.data(dg_off + 26);

            auto g_xy_yyzz = cbuffer.data(dg_off + 27);

            auto g_xy_yzzz = cbuffer.data(dg_off + 28);

            auto g_xy_zzzz = cbuffer.data(dg_off + 29);

            #pragma omp simd aligned(cd_x, g_xy_xxxx, g_xy_xxxy, g_xy_xxxz, g_xy_xxyy, g_xy_xxyz, g_xy_xxzz, g_xy_xyyy, g_xy_xyyz, g_xy_xyzz, g_xy_xzzz, g_xy_yyyy, g_xy_yyyz, g_xy_yyzz, g_xy_yzzz, g_xy_zzzz, g_y_xxxx, g_y_xxxxx, g_y_xxxxy, g_y_xxxxz, g_y_xxxy, g_y_xxxyy, g_y_xxxyz, g_y_xxxz, g_y_xxxzz, g_y_xxyy, g_y_xxyyy, g_y_xxyyz, g_y_xxyz, g_y_xxyzz, g_y_xxzz, g_y_xxzzz, g_y_xyyy, g_y_xyyyy, g_y_xyyyz, g_y_xyyz, g_y_xyyzz, g_y_xyzz, g_y_xyzzz, g_y_xzzz, g_y_xzzzz, g_y_yyyy, g_y_yyyz, g_y_yyzz, g_y_yzzz, g_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_xxxx[k] = -g_y_xxxx[k] * cd_x[k] + g_y_xxxxx[k];

                g_xy_xxxy[k] = -g_y_xxxy[k] * cd_x[k] + g_y_xxxxy[k];

                g_xy_xxxz[k] = -g_y_xxxz[k] * cd_x[k] + g_y_xxxxz[k];

                g_xy_xxyy[k] = -g_y_xxyy[k] * cd_x[k] + g_y_xxxyy[k];

                g_xy_xxyz[k] = -g_y_xxyz[k] * cd_x[k] + g_y_xxxyz[k];

                g_xy_xxzz[k] = -g_y_xxzz[k] * cd_x[k] + g_y_xxxzz[k];

                g_xy_xyyy[k] = -g_y_xyyy[k] * cd_x[k] + g_y_xxyyy[k];

                g_xy_xyyz[k] = -g_y_xyyz[k] * cd_x[k] + g_y_xxyyz[k];

                g_xy_xyzz[k] = -g_y_xyzz[k] * cd_x[k] + g_y_xxyzz[k];

                g_xy_xzzz[k] = -g_y_xzzz[k] * cd_x[k] + g_y_xxzzz[k];

                g_xy_yyyy[k] = -g_y_yyyy[k] * cd_x[k] + g_y_xyyyy[k];

                g_xy_yyyz[k] = -g_y_yyyz[k] * cd_x[k] + g_y_xyyyz[k];

                g_xy_yyzz[k] = -g_y_yyzz[k] * cd_x[k] + g_y_xyyzz[k];

                g_xy_yzzz[k] = -g_y_yzzz[k] * cd_x[k] + g_y_xyzzz[k];

                g_xy_zzzz[k] = -g_y_zzzz[k] * cd_x[k] + g_y_xzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_xz_xxxx = cbuffer.data(dg_off + 30);

            auto g_xz_xxxy = cbuffer.data(dg_off + 31);

            auto g_xz_xxxz = cbuffer.data(dg_off + 32);

            auto g_xz_xxyy = cbuffer.data(dg_off + 33);

            auto g_xz_xxyz = cbuffer.data(dg_off + 34);

            auto g_xz_xxzz = cbuffer.data(dg_off + 35);

            auto g_xz_xyyy = cbuffer.data(dg_off + 36);

            auto g_xz_xyyz = cbuffer.data(dg_off + 37);

            auto g_xz_xyzz = cbuffer.data(dg_off + 38);

            auto g_xz_xzzz = cbuffer.data(dg_off + 39);

            auto g_xz_yyyy = cbuffer.data(dg_off + 40);

            auto g_xz_yyyz = cbuffer.data(dg_off + 41);

            auto g_xz_yyzz = cbuffer.data(dg_off + 42);

            auto g_xz_yzzz = cbuffer.data(dg_off + 43);

            auto g_xz_zzzz = cbuffer.data(dg_off + 44);

            #pragma omp simd aligned(cd_x, g_xz_xxxx, g_xz_xxxy, g_xz_xxxz, g_xz_xxyy, g_xz_xxyz, g_xz_xxzz, g_xz_xyyy, g_xz_xyyz, g_xz_xyzz, g_xz_xzzz, g_xz_yyyy, g_xz_yyyz, g_xz_yyzz, g_xz_yzzz, g_xz_zzzz, g_z_xxxx, g_z_xxxxx, g_z_xxxxy, g_z_xxxxz, g_z_xxxy, g_z_xxxyy, g_z_xxxyz, g_z_xxxz, g_z_xxxzz, g_z_xxyy, g_z_xxyyy, g_z_xxyyz, g_z_xxyz, g_z_xxyzz, g_z_xxzz, g_z_xxzzz, g_z_xyyy, g_z_xyyyy, g_z_xyyyz, g_z_xyyz, g_z_xyyzz, g_z_xyzz, g_z_xyzzz, g_z_xzzz, g_z_xzzzz, g_z_yyyy, g_z_yyyz, g_z_yyzz, g_z_yzzz, g_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_xxxx[k] = -g_z_xxxx[k] * cd_x[k] + g_z_xxxxx[k];

                g_xz_xxxy[k] = -g_z_xxxy[k] * cd_x[k] + g_z_xxxxy[k];

                g_xz_xxxz[k] = -g_z_xxxz[k] * cd_x[k] + g_z_xxxxz[k];

                g_xz_xxyy[k] = -g_z_xxyy[k] * cd_x[k] + g_z_xxxyy[k];

                g_xz_xxyz[k] = -g_z_xxyz[k] * cd_x[k] + g_z_xxxyz[k];

                g_xz_xxzz[k] = -g_z_xxzz[k] * cd_x[k] + g_z_xxxzz[k];

                g_xz_xyyy[k] = -g_z_xyyy[k] * cd_x[k] + g_z_xxyyy[k];

                g_xz_xyyz[k] = -g_z_xyyz[k] * cd_x[k] + g_z_xxyyz[k];

                g_xz_xyzz[k] = -g_z_xyzz[k] * cd_x[k] + g_z_xxyzz[k];

                g_xz_xzzz[k] = -g_z_xzzz[k] * cd_x[k] + g_z_xxzzz[k];

                g_xz_yyyy[k] = -g_z_yyyy[k] * cd_x[k] + g_z_xyyyy[k];

                g_xz_yyyz[k] = -g_z_yyyz[k] * cd_x[k] + g_z_xyyyz[k];

                g_xz_yyzz[k] = -g_z_yyzz[k] * cd_x[k] + g_z_xyyzz[k];

                g_xz_yzzz[k] = -g_z_yzzz[k] * cd_x[k] + g_z_xyzzz[k];

                g_xz_zzzz[k] = -g_z_zzzz[k] * cd_x[k] + g_z_xzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_yy_xxxx = cbuffer.data(dg_off + 45);

            auto g_yy_xxxy = cbuffer.data(dg_off + 46);

            auto g_yy_xxxz = cbuffer.data(dg_off + 47);

            auto g_yy_xxyy = cbuffer.data(dg_off + 48);

            auto g_yy_xxyz = cbuffer.data(dg_off + 49);

            auto g_yy_xxzz = cbuffer.data(dg_off + 50);

            auto g_yy_xyyy = cbuffer.data(dg_off + 51);

            auto g_yy_xyyz = cbuffer.data(dg_off + 52);

            auto g_yy_xyzz = cbuffer.data(dg_off + 53);

            auto g_yy_xzzz = cbuffer.data(dg_off + 54);

            auto g_yy_yyyy = cbuffer.data(dg_off + 55);

            auto g_yy_yyyz = cbuffer.data(dg_off + 56);

            auto g_yy_yyzz = cbuffer.data(dg_off + 57);

            auto g_yy_yzzz = cbuffer.data(dg_off + 58);

            auto g_yy_zzzz = cbuffer.data(dg_off + 59);

            #pragma omp simd aligned(cd_y, g_y_xxxx, g_y_xxxxy, g_y_xxxy, g_y_xxxyy, g_y_xxxyz, g_y_xxxz, g_y_xxyy, g_y_xxyyy, g_y_xxyyz, g_y_xxyz, g_y_xxyzz, g_y_xxzz, g_y_xyyy, g_y_xyyyy, g_y_xyyyz, g_y_xyyz, g_y_xyyzz, g_y_xyzz, g_y_xyzzz, g_y_xzzz, g_y_yyyy, g_y_yyyyy, g_y_yyyyz, g_y_yyyz, g_y_yyyzz, g_y_yyzz, g_y_yyzzz, g_y_yzzz, g_y_yzzzz, g_y_zzzz, g_yy_xxxx, g_yy_xxxy, g_yy_xxxz, g_yy_xxyy, g_yy_xxyz, g_yy_xxzz, g_yy_xyyy, g_yy_xyyz, g_yy_xyzz, g_yy_xzzz, g_yy_yyyy, g_yy_yyyz, g_yy_yyzz, g_yy_yzzz, g_yy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_xxxx[k] = -g_y_xxxx[k] * cd_y[k] + g_y_xxxxy[k];

                g_yy_xxxy[k] = -g_y_xxxy[k] * cd_y[k] + g_y_xxxyy[k];

                g_yy_xxxz[k] = -g_y_xxxz[k] * cd_y[k] + g_y_xxxyz[k];

                g_yy_xxyy[k] = -g_y_xxyy[k] * cd_y[k] + g_y_xxyyy[k];

                g_yy_xxyz[k] = -g_y_xxyz[k] * cd_y[k] + g_y_xxyyz[k];

                g_yy_xxzz[k] = -g_y_xxzz[k] * cd_y[k] + g_y_xxyzz[k];

                g_yy_xyyy[k] = -g_y_xyyy[k] * cd_y[k] + g_y_xyyyy[k];

                g_yy_xyyz[k] = -g_y_xyyz[k] * cd_y[k] + g_y_xyyyz[k];

                g_yy_xyzz[k] = -g_y_xyzz[k] * cd_y[k] + g_y_xyyzz[k];

                g_yy_xzzz[k] = -g_y_xzzz[k] * cd_y[k] + g_y_xyzzz[k];

                g_yy_yyyy[k] = -g_y_yyyy[k] * cd_y[k] + g_y_yyyyy[k];

                g_yy_yyyz[k] = -g_y_yyyz[k] * cd_y[k] + g_y_yyyyz[k];

                g_yy_yyzz[k] = -g_y_yyzz[k] * cd_y[k] + g_y_yyyzz[k];

                g_yy_yzzz[k] = -g_y_yzzz[k] * cd_y[k] + g_y_yyzzz[k];

                g_yy_zzzz[k] = -g_y_zzzz[k] * cd_y[k] + g_y_yzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_yz_xxxx = cbuffer.data(dg_off + 60);

            auto g_yz_xxxy = cbuffer.data(dg_off + 61);

            auto g_yz_xxxz = cbuffer.data(dg_off + 62);

            auto g_yz_xxyy = cbuffer.data(dg_off + 63);

            auto g_yz_xxyz = cbuffer.data(dg_off + 64);

            auto g_yz_xxzz = cbuffer.data(dg_off + 65);

            auto g_yz_xyyy = cbuffer.data(dg_off + 66);

            auto g_yz_xyyz = cbuffer.data(dg_off + 67);

            auto g_yz_xyzz = cbuffer.data(dg_off + 68);

            auto g_yz_xzzz = cbuffer.data(dg_off + 69);

            auto g_yz_yyyy = cbuffer.data(dg_off + 70);

            auto g_yz_yyyz = cbuffer.data(dg_off + 71);

            auto g_yz_yyzz = cbuffer.data(dg_off + 72);

            auto g_yz_yzzz = cbuffer.data(dg_off + 73);

            auto g_yz_zzzz = cbuffer.data(dg_off + 74);

            #pragma omp simd aligned(cd_y, g_yz_xxxx, g_yz_xxxy, g_yz_xxxz, g_yz_xxyy, g_yz_xxyz, g_yz_xxzz, g_yz_xyyy, g_yz_xyyz, g_yz_xyzz, g_yz_xzzz, g_yz_yyyy, g_yz_yyyz, g_yz_yyzz, g_yz_yzzz, g_yz_zzzz, g_z_xxxx, g_z_xxxxy, g_z_xxxy, g_z_xxxyy, g_z_xxxyz, g_z_xxxz, g_z_xxyy, g_z_xxyyy, g_z_xxyyz, g_z_xxyz, g_z_xxyzz, g_z_xxzz, g_z_xyyy, g_z_xyyyy, g_z_xyyyz, g_z_xyyz, g_z_xyyzz, g_z_xyzz, g_z_xyzzz, g_z_xzzz, g_z_yyyy, g_z_yyyyy, g_z_yyyyz, g_z_yyyz, g_z_yyyzz, g_z_yyzz, g_z_yyzzz, g_z_yzzz, g_z_yzzzz, g_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_xxxx[k] = -g_z_xxxx[k] * cd_y[k] + g_z_xxxxy[k];

                g_yz_xxxy[k] = -g_z_xxxy[k] * cd_y[k] + g_z_xxxyy[k];

                g_yz_xxxz[k] = -g_z_xxxz[k] * cd_y[k] + g_z_xxxyz[k];

                g_yz_xxyy[k] = -g_z_xxyy[k] * cd_y[k] + g_z_xxyyy[k];

                g_yz_xxyz[k] = -g_z_xxyz[k] * cd_y[k] + g_z_xxyyz[k];

                g_yz_xxzz[k] = -g_z_xxzz[k] * cd_y[k] + g_z_xxyzz[k];

                g_yz_xyyy[k] = -g_z_xyyy[k] * cd_y[k] + g_z_xyyyy[k];

                g_yz_xyyz[k] = -g_z_xyyz[k] * cd_y[k] + g_z_xyyyz[k];

                g_yz_xyzz[k] = -g_z_xyzz[k] * cd_y[k] + g_z_xyyzz[k];

                g_yz_xzzz[k] = -g_z_xzzz[k] * cd_y[k] + g_z_xyzzz[k];

                g_yz_yyyy[k] = -g_z_yyyy[k] * cd_y[k] + g_z_yyyyy[k];

                g_yz_yyyz[k] = -g_z_yyyz[k] * cd_y[k] + g_z_yyyyz[k];

                g_yz_yyzz[k] = -g_z_yyzz[k] * cd_y[k] + g_z_yyyzz[k];

                g_yz_yzzz[k] = -g_z_yzzz[k] * cd_y[k] + g_z_yyzzz[k];

                g_yz_zzzz[k] = -g_z_zzzz[k] * cd_y[k] + g_z_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_zz_xxxx = cbuffer.data(dg_off + 75);

            auto g_zz_xxxy = cbuffer.data(dg_off + 76);

            auto g_zz_xxxz = cbuffer.data(dg_off + 77);

            auto g_zz_xxyy = cbuffer.data(dg_off + 78);

            auto g_zz_xxyz = cbuffer.data(dg_off + 79);

            auto g_zz_xxzz = cbuffer.data(dg_off + 80);

            auto g_zz_xyyy = cbuffer.data(dg_off + 81);

            auto g_zz_xyyz = cbuffer.data(dg_off + 82);

            auto g_zz_xyzz = cbuffer.data(dg_off + 83);

            auto g_zz_xzzz = cbuffer.data(dg_off + 84);

            auto g_zz_yyyy = cbuffer.data(dg_off + 85);

            auto g_zz_yyyz = cbuffer.data(dg_off + 86);

            auto g_zz_yyzz = cbuffer.data(dg_off + 87);

            auto g_zz_yzzz = cbuffer.data(dg_off + 88);

            auto g_zz_zzzz = cbuffer.data(dg_off + 89);

            #pragma omp simd aligned(cd_z, g_z_xxxx, g_z_xxxxz, g_z_xxxy, g_z_xxxyz, g_z_xxxz, g_z_xxxzz, g_z_xxyy, g_z_xxyyz, g_z_xxyz, g_z_xxyzz, g_z_xxzz, g_z_xxzzz, g_z_xyyy, g_z_xyyyz, g_z_xyyz, g_z_xyyzz, g_z_xyzz, g_z_xyzzz, g_z_xzzz, g_z_xzzzz, g_z_yyyy, g_z_yyyyz, g_z_yyyz, g_z_yyyzz, g_z_yyzz, g_z_yyzzz, g_z_yzzz, g_z_yzzzz, g_z_zzzz, g_z_zzzzz, g_zz_xxxx, g_zz_xxxy, g_zz_xxxz, g_zz_xxyy, g_zz_xxyz, g_zz_xxzz, g_zz_xyyy, g_zz_xyyz, g_zz_xyzz, g_zz_xzzz, g_zz_yyyy, g_zz_yyyz, g_zz_yyzz, g_zz_yzzz, g_zz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_xxxx[k] = -g_z_xxxx[k] * cd_z[k] + g_z_xxxxz[k];

                g_zz_xxxy[k] = -g_z_xxxy[k] * cd_z[k] + g_z_xxxyz[k];

                g_zz_xxxz[k] = -g_z_xxxz[k] * cd_z[k] + g_z_xxxzz[k];

                g_zz_xxyy[k] = -g_z_xxyy[k] * cd_z[k] + g_z_xxyyz[k];

                g_zz_xxyz[k] = -g_z_xxyz[k] * cd_z[k] + g_z_xxyzz[k];

                g_zz_xxzz[k] = -g_z_xxzz[k] * cd_z[k] + g_z_xxzzz[k];

                g_zz_xyyy[k] = -g_z_xyyy[k] * cd_z[k] + g_z_xyyyz[k];

                g_zz_xyyz[k] = -g_z_xyyz[k] * cd_z[k] + g_z_xyyzz[k];

                g_zz_xyzz[k] = -g_z_xyzz[k] * cd_z[k] + g_z_xyzzz[k];

                g_zz_xzzz[k] = -g_z_xzzz[k] * cd_z[k] + g_z_xzzzz[k];

                g_zz_yyyy[k] = -g_z_yyyy[k] * cd_z[k] + g_z_yyyyz[k];

                g_zz_yyyz[k] = -g_z_yyyz[k] * cd_z[k] + g_z_yyyzz[k];

                g_zz_yyzz[k] = -g_z_yyzz[k] * cd_z[k] + g_z_yyzzz[k];

                g_zz_yzzz[k] = -g_z_yzzz[k] * cd_z[k] + g_z_yzzzz[k];

                g_zz_zzzz[k] = -g_z_zzzz[k] * cd_z[k] + g_z_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

