#include "ElectronRepulsionContrRecXXPG.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_hrr_electron_repulsion_xxpg(CSimdArray<double>& contr_buffer_xxpg,
                                     const CSimdArray<double>& contr_buffer_xxsg,
                                     const CSimdArray<double>& contr_buffer_xxsh,
                                     const double* cd_x,
                                     const double* cd_y,
                                     const double* cd_z,
                                     const int a_angmom,
                                     const int b_angmom) -> void
{
    const auto ndims = contr_buffer_xxpg.number_of_columns();

    const auto acomps = tensor::number_of_cartesian_components(a_angmom);

    const auto bcomps = tensor::number_of_cartesian_components(b_angmom);

    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
            /// Set up components of auxilary buffer : contr_buffer_xxsg

            const auto sg_off = (i * bcomps + j) * 15;

            auto g_0_xxxx = contr_buffer_xxsg[sg_off + 0];

            auto g_0_xxxy = contr_buffer_xxsg[sg_off + 1];

            auto g_0_xxxz = contr_buffer_xxsg[sg_off + 2];

            auto g_0_xxyy = contr_buffer_xxsg[sg_off + 3];

            auto g_0_xxyz = contr_buffer_xxsg[sg_off + 4];

            auto g_0_xxzz = contr_buffer_xxsg[sg_off + 5];

            auto g_0_xyyy = contr_buffer_xxsg[sg_off + 6];

            auto g_0_xyyz = contr_buffer_xxsg[sg_off + 7];

            auto g_0_xyzz = contr_buffer_xxsg[sg_off + 8];

            auto g_0_xzzz = contr_buffer_xxsg[sg_off + 9];

            auto g_0_yyyy = contr_buffer_xxsg[sg_off + 10];

            auto g_0_yyyz = contr_buffer_xxsg[sg_off + 11];

            auto g_0_yyzz = contr_buffer_xxsg[sg_off + 12];

            auto g_0_yzzz = contr_buffer_xxsg[sg_off + 13];

            auto g_0_zzzz = contr_buffer_xxsg[sg_off + 14];

            /// Set up components of auxilary buffer : contr_buffer_xxsh

            const auto sh_off = (i * bcomps + j) * 21;

            auto g_0_xxxxx = contr_buffer_xxsh[sh_off + 0];

            auto g_0_xxxxy = contr_buffer_xxsh[sh_off + 1];

            auto g_0_xxxxz = contr_buffer_xxsh[sh_off + 2];

            auto g_0_xxxyy = contr_buffer_xxsh[sh_off + 3];

            auto g_0_xxxyz = contr_buffer_xxsh[sh_off + 4];

            auto g_0_xxxzz = contr_buffer_xxsh[sh_off + 5];

            auto g_0_xxyyy = contr_buffer_xxsh[sh_off + 6];

            auto g_0_xxyyz = contr_buffer_xxsh[sh_off + 7];

            auto g_0_xxyzz = contr_buffer_xxsh[sh_off + 8];

            auto g_0_xxzzz = contr_buffer_xxsh[sh_off + 9];

            auto g_0_xyyyy = contr_buffer_xxsh[sh_off + 10];

            auto g_0_xyyyz = contr_buffer_xxsh[sh_off + 11];

            auto g_0_xyyzz = contr_buffer_xxsh[sh_off + 12];

            auto g_0_xyzzz = contr_buffer_xxsh[sh_off + 13];

            auto g_0_xzzzz = contr_buffer_xxsh[sh_off + 14];

            auto g_0_yyyyy = contr_buffer_xxsh[sh_off + 15];

            auto g_0_yyyyz = contr_buffer_xxsh[sh_off + 16];

            auto g_0_yyyzz = contr_buffer_xxsh[sh_off + 17];

            auto g_0_yyzzz = contr_buffer_xxsh[sh_off + 18];

            auto g_0_yzzzz = contr_buffer_xxsh[sh_off + 19];

            auto g_0_zzzzz = contr_buffer_xxsh[sh_off + 20];

            /// set up bra offset for contr_buffer_xxpg

            const auto pg_off = (i * bcomps + j) * 45;

            /// Set up 0-15 components of targeted buffer : contr_buffer_xxpg

            auto g_x_xxxx = contr_buffer_xxpg[pg_off + 0];

            auto g_x_xxxy = contr_buffer_xxpg[pg_off + 1];

            auto g_x_xxxz = contr_buffer_xxpg[pg_off + 2];

            auto g_x_xxyy = contr_buffer_xxpg[pg_off + 3];

            auto g_x_xxyz = contr_buffer_xxpg[pg_off + 4];

            auto g_x_xxzz = contr_buffer_xxpg[pg_off + 5];

            auto g_x_xyyy = contr_buffer_xxpg[pg_off + 6];

            auto g_x_xyyz = contr_buffer_xxpg[pg_off + 7];

            auto g_x_xyzz = contr_buffer_xxpg[pg_off + 8];

            auto g_x_xzzz = contr_buffer_xxpg[pg_off + 9];

            auto g_x_yyyy = contr_buffer_xxpg[pg_off + 10];

            auto g_x_yyyz = contr_buffer_xxpg[pg_off + 11];

            auto g_x_yyzz = contr_buffer_xxpg[pg_off + 12];

            auto g_x_yzzz = contr_buffer_xxpg[pg_off + 13];

            auto g_x_zzzz = contr_buffer_xxpg[pg_off + 14];

            #pragma omp simd aligned(cd_x, g_0_xxxx, g_0_xxxxx, g_0_xxxxy, g_0_xxxxz, g_0_xxxy, g_0_xxxyy, g_0_xxxyz, g_0_xxxz, g_0_xxxzz, g_0_xxyy, g_0_xxyyy, g_0_xxyyz, g_0_xxyz, g_0_xxyzz, g_0_xxzz, g_0_xxzzz, g_0_xyyy, g_0_xyyyy, g_0_xyyyz, g_0_xyyz, g_0_xyyzz, g_0_xyzz, g_0_xyzzz, g_0_xzzz, g_0_xzzzz, g_0_yyyy, g_0_yyyz, g_0_yyzz, g_0_yzzz, g_0_zzzz, g_x_xxxx, g_x_xxxy, g_x_xxxz, g_x_xxyy, g_x_xxyz, g_x_xxzz, g_x_xyyy, g_x_xyyz, g_x_xyzz, g_x_xzzz, g_x_yyyy, g_x_yyyz, g_x_yyzz, g_x_yzzz, g_x_zzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_x_xxxx[k] = -g_0_xxxx[k] * cd_x[k] + g_0_xxxxx[k];

                g_x_xxxy[k] = -g_0_xxxy[k] * cd_x[k] + g_0_xxxxy[k];

                g_x_xxxz[k] = -g_0_xxxz[k] * cd_x[k] + g_0_xxxxz[k];

                g_x_xxyy[k] = -g_0_xxyy[k] * cd_x[k] + g_0_xxxyy[k];

                g_x_xxyz[k] = -g_0_xxyz[k] * cd_x[k] + g_0_xxxyz[k];

                g_x_xxzz[k] = -g_0_xxzz[k] * cd_x[k] + g_0_xxxzz[k];

                g_x_xyyy[k] = -g_0_xyyy[k] * cd_x[k] + g_0_xxyyy[k];

                g_x_xyyz[k] = -g_0_xyyz[k] * cd_x[k] + g_0_xxyyz[k];

                g_x_xyzz[k] = -g_0_xyzz[k] * cd_x[k] + g_0_xxyzz[k];

                g_x_xzzz[k] = -g_0_xzzz[k] * cd_x[k] + g_0_xxzzz[k];

                g_x_yyyy[k] = -g_0_yyyy[k] * cd_x[k] + g_0_xyyyy[k];

                g_x_yyyz[k] = -g_0_yyyz[k] * cd_x[k] + g_0_xyyyz[k];

                g_x_yyzz[k] = -g_0_yyzz[k] * cd_x[k] + g_0_xyyzz[k];

                g_x_yzzz[k] = -g_0_yzzz[k] * cd_x[k] + g_0_xyzzz[k];

                g_x_zzzz[k] = -g_0_zzzz[k] * cd_x[k] + g_0_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : contr_buffer_xxpg

            auto g_y_xxxx = contr_buffer_xxpg[pg_off + 15];

            auto g_y_xxxy = contr_buffer_xxpg[pg_off + 16];

            auto g_y_xxxz = contr_buffer_xxpg[pg_off + 17];

            auto g_y_xxyy = contr_buffer_xxpg[pg_off + 18];

            auto g_y_xxyz = contr_buffer_xxpg[pg_off + 19];

            auto g_y_xxzz = contr_buffer_xxpg[pg_off + 20];

            auto g_y_xyyy = contr_buffer_xxpg[pg_off + 21];

            auto g_y_xyyz = contr_buffer_xxpg[pg_off + 22];

            auto g_y_xyzz = contr_buffer_xxpg[pg_off + 23];

            auto g_y_xzzz = contr_buffer_xxpg[pg_off + 24];

            auto g_y_yyyy = contr_buffer_xxpg[pg_off + 25];

            auto g_y_yyyz = contr_buffer_xxpg[pg_off + 26];

            auto g_y_yyzz = contr_buffer_xxpg[pg_off + 27];

            auto g_y_yzzz = contr_buffer_xxpg[pg_off + 28];

            auto g_y_zzzz = contr_buffer_xxpg[pg_off + 29];

            #pragma omp simd aligned(cd_y, g_0_xxxx, g_0_xxxxy, g_0_xxxy, g_0_xxxyy, g_0_xxxyz, g_0_xxxz, g_0_xxyy, g_0_xxyyy, g_0_xxyyz, g_0_xxyz, g_0_xxyzz, g_0_xxzz, g_0_xyyy, g_0_xyyyy, g_0_xyyyz, g_0_xyyz, g_0_xyyzz, g_0_xyzz, g_0_xyzzz, g_0_xzzz, g_0_yyyy, g_0_yyyyy, g_0_yyyyz, g_0_yyyz, g_0_yyyzz, g_0_yyzz, g_0_yyzzz, g_0_yzzz, g_0_yzzzz, g_0_zzzz, g_y_xxxx, g_y_xxxy, g_y_xxxz, g_y_xxyy, g_y_xxyz, g_y_xxzz, g_y_xyyy, g_y_xyyz, g_y_xyzz, g_y_xzzz, g_y_yyyy, g_y_yyyz, g_y_yyzz, g_y_yzzz, g_y_zzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_y_xxxx[k] = -g_0_xxxx[k] * cd_y[k] + g_0_xxxxy[k];

                g_y_xxxy[k] = -g_0_xxxy[k] * cd_y[k] + g_0_xxxyy[k];

                g_y_xxxz[k] = -g_0_xxxz[k] * cd_y[k] + g_0_xxxyz[k];

                g_y_xxyy[k] = -g_0_xxyy[k] * cd_y[k] + g_0_xxyyy[k];

                g_y_xxyz[k] = -g_0_xxyz[k] * cd_y[k] + g_0_xxyyz[k];

                g_y_xxzz[k] = -g_0_xxzz[k] * cd_y[k] + g_0_xxyzz[k];

                g_y_xyyy[k] = -g_0_xyyy[k] * cd_y[k] + g_0_xyyyy[k];

                g_y_xyyz[k] = -g_0_xyyz[k] * cd_y[k] + g_0_xyyyz[k];

                g_y_xyzz[k] = -g_0_xyzz[k] * cd_y[k] + g_0_xyyzz[k];

                g_y_xzzz[k] = -g_0_xzzz[k] * cd_y[k] + g_0_xyzzz[k];

                g_y_yyyy[k] = -g_0_yyyy[k] * cd_y[k] + g_0_yyyyy[k];

                g_y_yyyz[k] = -g_0_yyyz[k] * cd_y[k] + g_0_yyyyz[k];

                g_y_yyzz[k] = -g_0_yyzz[k] * cd_y[k] + g_0_yyyzz[k];

                g_y_yzzz[k] = -g_0_yzzz[k] * cd_y[k] + g_0_yyzzz[k];

                g_y_zzzz[k] = -g_0_zzzz[k] * cd_y[k] + g_0_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : contr_buffer_xxpg

            auto g_z_xxxx = contr_buffer_xxpg[pg_off + 30];

            auto g_z_xxxy = contr_buffer_xxpg[pg_off + 31];

            auto g_z_xxxz = contr_buffer_xxpg[pg_off + 32];

            auto g_z_xxyy = contr_buffer_xxpg[pg_off + 33];

            auto g_z_xxyz = contr_buffer_xxpg[pg_off + 34];

            auto g_z_xxzz = contr_buffer_xxpg[pg_off + 35];

            auto g_z_xyyy = contr_buffer_xxpg[pg_off + 36];

            auto g_z_xyyz = contr_buffer_xxpg[pg_off + 37];

            auto g_z_xyzz = contr_buffer_xxpg[pg_off + 38];

            auto g_z_xzzz = contr_buffer_xxpg[pg_off + 39];

            auto g_z_yyyy = contr_buffer_xxpg[pg_off + 40];

            auto g_z_yyyz = contr_buffer_xxpg[pg_off + 41];

            auto g_z_yyzz = contr_buffer_xxpg[pg_off + 42];

            auto g_z_yzzz = contr_buffer_xxpg[pg_off + 43];

            auto g_z_zzzz = contr_buffer_xxpg[pg_off + 44];

            #pragma omp simd aligned(cd_z, g_0_xxxx, g_0_xxxxz, g_0_xxxy, g_0_xxxyz, g_0_xxxz, g_0_xxxzz, g_0_xxyy, g_0_xxyyz, g_0_xxyz, g_0_xxyzz, g_0_xxzz, g_0_xxzzz, g_0_xyyy, g_0_xyyyz, g_0_xyyz, g_0_xyyzz, g_0_xyzz, g_0_xyzzz, g_0_xzzz, g_0_xzzzz, g_0_yyyy, g_0_yyyyz, g_0_yyyz, g_0_yyyzz, g_0_yyzz, g_0_yyzzz, g_0_yzzz, g_0_yzzzz, g_0_zzzz, g_0_zzzzz, g_z_xxxx, g_z_xxxy, g_z_xxxz, g_z_xxyy, g_z_xxyz, g_z_xxzz, g_z_xyyy, g_z_xyyz, g_z_xyzz, g_z_xzzz, g_z_yyyy, g_z_yyyz, g_z_yyzz, g_z_yzzz, g_z_zzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_z_xxxx[k] = -g_0_xxxx[k] * cd_z[k] + g_0_xxxxz[k];

                g_z_xxxy[k] = -g_0_xxxy[k] * cd_z[k] + g_0_xxxyz[k];

                g_z_xxxz[k] = -g_0_xxxz[k] * cd_z[k] + g_0_xxxzz[k];

                g_z_xxyy[k] = -g_0_xxyy[k] * cd_z[k] + g_0_xxyyz[k];

                g_z_xxyz[k] = -g_0_xxyz[k] * cd_z[k] + g_0_xxyzz[k];

                g_z_xxzz[k] = -g_0_xxzz[k] * cd_z[k] + g_0_xxzzz[k];

                g_z_xyyy[k] = -g_0_xyyy[k] * cd_z[k] + g_0_xyyyz[k];

                g_z_xyyz[k] = -g_0_xyyz[k] * cd_z[k] + g_0_xyyzz[k];

                g_z_xyzz[k] = -g_0_xyzz[k] * cd_z[k] + g_0_xyzzz[k];

                g_z_xzzz[k] = -g_0_xzzz[k] * cd_z[k] + g_0_xzzzz[k];

                g_z_yyyy[k] = -g_0_yyyy[k] * cd_z[k] + g_0_yyyyz[k];

                g_z_yyyz[k] = -g_0_yyyz[k] * cd_z[k] + g_0_yyyzz[k];

                g_z_yyzz[k] = -g_0_yyzz[k] * cd_z[k] + g_0_yyzzz[k];

                g_z_yzzz[k] = -g_0_yzzz[k] * cd_z[k] + g_0_yzzzz[k];

                g_z_zzzz[k] = -g_0_zzzz[k] * cd_z[k] + g_0_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

