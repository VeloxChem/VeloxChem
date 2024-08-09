#include "ElectronRepulsionContrRecXXPH.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_hrr_electron_repulsion_xxph(CSimdArray<double>& contr_buffer_xxph,
                                     const CSimdArray<double>& contr_buffer_xxsh,
                                     const CSimdArray<double>& contr_buffer_xxsi,
                                     const double* cd_x,
                                     const double* cd_y,
                                     const double* cd_z,
                                     const int a_angmom,
                                     const int b_angmom) -> void
{
    const auto ndims = contr_buffer_xxph.number_of_columns();

    const auto acomps = tensor::number_of_cartesian_components(a_angmom);

    const auto bcomps = tensor::number_of_cartesian_components(b_angmom);

    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
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

            /// Set up components of auxilary buffer : contr_buffer_xxsi

            const auto si_off = (i * bcomps + j) * 28;

            auto g_0_xxxxxx = contr_buffer_xxsi[si_off + 0];

            auto g_0_xxxxxy = contr_buffer_xxsi[si_off + 1];

            auto g_0_xxxxxz = contr_buffer_xxsi[si_off + 2];

            auto g_0_xxxxyy = contr_buffer_xxsi[si_off + 3];

            auto g_0_xxxxyz = contr_buffer_xxsi[si_off + 4];

            auto g_0_xxxxzz = contr_buffer_xxsi[si_off + 5];

            auto g_0_xxxyyy = contr_buffer_xxsi[si_off + 6];

            auto g_0_xxxyyz = contr_buffer_xxsi[si_off + 7];

            auto g_0_xxxyzz = contr_buffer_xxsi[si_off + 8];

            auto g_0_xxxzzz = contr_buffer_xxsi[si_off + 9];

            auto g_0_xxyyyy = contr_buffer_xxsi[si_off + 10];

            auto g_0_xxyyyz = contr_buffer_xxsi[si_off + 11];

            auto g_0_xxyyzz = contr_buffer_xxsi[si_off + 12];

            auto g_0_xxyzzz = contr_buffer_xxsi[si_off + 13];

            auto g_0_xxzzzz = contr_buffer_xxsi[si_off + 14];

            auto g_0_xyyyyy = contr_buffer_xxsi[si_off + 15];

            auto g_0_xyyyyz = contr_buffer_xxsi[si_off + 16];

            auto g_0_xyyyzz = contr_buffer_xxsi[si_off + 17];

            auto g_0_xyyzzz = contr_buffer_xxsi[si_off + 18];

            auto g_0_xyzzzz = contr_buffer_xxsi[si_off + 19];

            auto g_0_xzzzzz = contr_buffer_xxsi[si_off + 20];

            auto g_0_yyyyyy = contr_buffer_xxsi[si_off + 21];

            auto g_0_yyyyyz = contr_buffer_xxsi[si_off + 22];

            auto g_0_yyyyzz = contr_buffer_xxsi[si_off + 23];

            auto g_0_yyyzzz = contr_buffer_xxsi[si_off + 24];

            auto g_0_yyzzzz = contr_buffer_xxsi[si_off + 25];

            auto g_0_yzzzzz = contr_buffer_xxsi[si_off + 26];

            auto g_0_zzzzzz = contr_buffer_xxsi[si_off + 27];

            /// set up bra offset for contr_buffer_xxph

            const auto ph_off = (i * bcomps + j) * 63;

            /// Set up 0-21 components of targeted buffer : contr_buffer_xxph

            auto g_x_xxxxx = contr_buffer_xxph[ph_off + 0];

            auto g_x_xxxxy = contr_buffer_xxph[ph_off + 1];

            auto g_x_xxxxz = contr_buffer_xxph[ph_off + 2];

            auto g_x_xxxyy = contr_buffer_xxph[ph_off + 3];

            auto g_x_xxxyz = contr_buffer_xxph[ph_off + 4];

            auto g_x_xxxzz = contr_buffer_xxph[ph_off + 5];

            auto g_x_xxyyy = contr_buffer_xxph[ph_off + 6];

            auto g_x_xxyyz = contr_buffer_xxph[ph_off + 7];

            auto g_x_xxyzz = contr_buffer_xxph[ph_off + 8];

            auto g_x_xxzzz = contr_buffer_xxph[ph_off + 9];

            auto g_x_xyyyy = contr_buffer_xxph[ph_off + 10];

            auto g_x_xyyyz = contr_buffer_xxph[ph_off + 11];

            auto g_x_xyyzz = contr_buffer_xxph[ph_off + 12];

            auto g_x_xyzzz = contr_buffer_xxph[ph_off + 13];

            auto g_x_xzzzz = contr_buffer_xxph[ph_off + 14];

            auto g_x_yyyyy = contr_buffer_xxph[ph_off + 15];

            auto g_x_yyyyz = contr_buffer_xxph[ph_off + 16];

            auto g_x_yyyzz = contr_buffer_xxph[ph_off + 17];

            auto g_x_yyzzz = contr_buffer_xxph[ph_off + 18];

            auto g_x_yzzzz = contr_buffer_xxph[ph_off + 19];

            auto g_x_zzzzz = contr_buffer_xxph[ph_off + 20];

            #pragma omp simd aligned(cd_x, g_0_xxxxx, g_0_xxxxxx, g_0_xxxxxy, g_0_xxxxxz, g_0_xxxxy, g_0_xxxxyy, g_0_xxxxyz, g_0_xxxxz, g_0_xxxxzz, g_0_xxxyy, g_0_xxxyyy, g_0_xxxyyz, g_0_xxxyz, g_0_xxxyzz, g_0_xxxzz, g_0_xxxzzz, g_0_xxyyy, g_0_xxyyyy, g_0_xxyyyz, g_0_xxyyz, g_0_xxyyzz, g_0_xxyzz, g_0_xxyzzz, g_0_xxzzz, g_0_xxzzzz, g_0_xyyyy, g_0_xyyyyy, g_0_xyyyyz, g_0_xyyyz, g_0_xyyyzz, g_0_xyyzz, g_0_xyyzzz, g_0_xyzzz, g_0_xyzzzz, g_0_xzzzz, g_0_xzzzzz, g_0_yyyyy, g_0_yyyyz, g_0_yyyzz, g_0_yyzzz, g_0_yzzzz, g_0_zzzzz, g_x_xxxxx, g_x_xxxxy, g_x_xxxxz, g_x_xxxyy, g_x_xxxyz, g_x_xxxzz, g_x_xxyyy, g_x_xxyyz, g_x_xxyzz, g_x_xxzzz, g_x_xyyyy, g_x_xyyyz, g_x_xyyzz, g_x_xyzzz, g_x_xzzzz, g_x_yyyyy, g_x_yyyyz, g_x_yyyzz, g_x_yyzzz, g_x_yzzzz, g_x_zzzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_x_xxxxx[k] = -g_0_xxxxx[k] * cd_x[k] + g_0_xxxxxx[k];

                g_x_xxxxy[k] = -g_0_xxxxy[k] * cd_x[k] + g_0_xxxxxy[k];

                g_x_xxxxz[k] = -g_0_xxxxz[k] * cd_x[k] + g_0_xxxxxz[k];

                g_x_xxxyy[k] = -g_0_xxxyy[k] * cd_x[k] + g_0_xxxxyy[k];

                g_x_xxxyz[k] = -g_0_xxxyz[k] * cd_x[k] + g_0_xxxxyz[k];

                g_x_xxxzz[k] = -g_0_xxxzz[k] * cd_x[k] + g_0_xxxxzz[k];

                g_x_xxyyy[k] = -g_0_xxyyy[k] * cd_x[k] + g_0_xxxyyy[k];

                g_x_xxyyz[k] = -g_0_xxyyz[k] * cd_x[k] + g_0_xxxyyz[k];

                g_x_xxyzz[k] = -g_0_xxyzz[k] * cd_x[k] + g_0_xxxyzz[k];

                g_x_xxzzz[k] = -g_0_xxzzz[k] * cd_x[k] + g_0_xxxzzz[k];

                g_x_xyyyy[k] = -g_0_xyyyy[k] * cd_x[k] + g_0_xxyyyy[k];

                g_x_xyyyz[k] = -g_0_xyyyz[k] * cd_x[k] + g_0_xxyyyz[k];

                g_x_xyyzz[k] = -g_0_xyyzz[k] * cd_x[k] + g_0_xxyyzz[k];

                g_x_xyzzz[k] = -g_0_xyzzz[k] * cd_x[k] + g_0_xxyzzz[k];

                g_x_xzzzz[k] = -g_0_xzzzz[k] * cd_x[k] + g_0_xxzzzz[k];

                g_x_yyyyy[k] = -g_0_yyyyy[k] * cd_x[k] + g_0_xyyyyy[k];

                g_x_yyyyz[k] = -g_0_yyyyz[k] * cd_x[k] + g_0_xyyyyz[k];

                g_x_yyyzz[k] = -g_0_yyyzz[k] * cd_x[k] + g_0_xyyyzz[k];

                g_x_yyzzz[k] = -g_0_yyzzz[k] * cd_x[k] + g_0_xyyzzz[k];

                g_x_yzzzz[k] = -g_0_yzzzz[k] * cd_x[k] + g_0_xyzzzz[k];

                g_x_zzzzz[k] = -g_0_zzzzz[k] * cd_x[k] + g_0_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : contr_buffer_xxph

            auto g_y_xxxxx = contr_buffer_xxph[ph_off + 21];

            auto g_y_xxxxy = contr_buffer_xxph[ph_off + 22];

            auto g_y_xxxxz = contr_buffer_xxph[ph_off + 23];

            auto g_y_xxxyy = contr_buffer_xxph[ph_off + 24];

            auto g_y_xxxyz = contr_buffer_xxph[ph_off + 25];

            auto g_y_xxxzz = contr_buffer_xxph[ph_off + 26];

            auto g_y_xxyyy = contr_buffer_xxph[ph_off + 27];

            auto g_y_xxyyz = contr_buffer_xxph[ph_off + 28];

            auto g_y_xxyzz = contr_buffer_xxph[ph_off + 29];

            auto g_y_xxzzz = contr_buffer_xxph[ph_off + 30];

            auto g_y_xyyyy = contr_buffer_xxph[ph_off + 31];

            auto g_y_xyyyz = contr_buffer_xxph[ph_off + 32];

            auto g_y_xyyzz = contr_buffer_xxph[ph_off + 33];

            auto g_y_xyzzz = contr_buffer_xxph[ph_off + 34];

            auto g_y_xzzzz = contr_buffer_xxph[ph_off + 35];

            auto g_y_yyyyy = contr_buffer_xxph[ph_off + 36];

            auto g_y_yyyyz = contr_buffer_xxph[ph_off + 37];

            auto g_y_yyyzz = contr_buffer_xxph[ph_off + 38];

            auto g_y_yyzzz = contr_buffer_xxph[ph_off + 39];

            auto g_y_yzzzz = contr_buffer_xxph[ph_off + 40];

            auto g_y_zzzzz = contr_buffer_xxph[ph_off + 41];

            #pragma omp simd aligned(cd_y, g_0_xxxxx, g_0_xxxxxy, g_0_xxxxy, g_0_xxxxyy, g_0_xxxxyz, g_0_xxxxz, g_0_xxxyy, g_0_xxxyyy, g_0_xxxyyz, g_0_xxxyz, g_0_xxxyzz, g_0_xxxzz, g_0_xxyyy, g_0_xxyyyy, g_0_xxyyyz, g_0_xxyyz, g_0_xxyyzz, g_0_xxyzz, g_0_xxyzzz, g_0_xxzzz, g_0_xyyyy, g_0_xyyyyy, g_0_xyyyyz, g_0_xyyyz, g_0_xyyyzz, g_0_xyyzz, g_0_xyyzzz, g_0_xyzzz, g_0_xyzzzz, g_0_xzzzz, g_0_yyyyy, g_0_yyyyyy, g_0_yyyyyz, g_0_yyyyz, g_0_yyyyzz, g_0_yyyzz, g_0_yyyzzz, g_0_yyzzz, g_0_yyzzzz, g_0_yzzzz, g_0_yzzzzz, g_0_zzzzz, g_y_xxxxx, g_y_xxxxy, g_y_xxxxz, g_y_xxxyy, g_y_xxxyz, g_y_xxxzz, g_y_xxyyy, g_y_xxyyz, g_y_xxyzz, g_y_xxzzz, g_y_xyyyy, g_y_xyyyz, g_y_xyyzz, g_y_xyzzz, g_y_xzzzz, g_y_yyyyy, g_y_yyyyz, g_y_yyyzz, g_y_yyzzz, g_y_yzzzz, g_y_zzzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_y_xxxxx[k] = -g_0_xxxxx[k] * cd_y[k] + g_0_xxxxxy[k];

                g_y_xxxxy[k] = -g_0_xxxxy[k] * cd_y[k] + g_0_xxxxyy[k];

                g_y_xxxxz[k] = -g_0_xxxxz[k] * cd_y[k] + g_0_xxxxyz[k];

                g_y_xxxyy[k] = -g_0_xxxyy[k] * cd_y[k] + g_0_xxxyyy[k];

                g_y_xxxyz[k] = -g_0_xxxyz[k] * cd_y[k] + g_0_xxxyyz[k];

                g_y_xxxzz[k] = -g_0_xxxzz[k] * cd_y[k] + g_0_xxxyzz[k];

                g_y_xxyyy[k] = -g_0_xxyyy[k] * cd_y[k] + g_0_xxyyyy[k];

                g_y_xxyyz[k] = -g_0_xxyyz[k] * cd_y[k] + g_0_xxyyyz[k];

                g_y_xxyzz[k] = -g_0_xxyzz[k] * cd_y[k] + g_0_xxyyzz[k];

                g_y_xxzzz[k] = -g_0_xxzzz[k] * cd_y[k] + g_0_xxyzzz[k];

                g_y_xyyyy[k] = -g_0_xyyyy[k] * cd_y[k] + g_0_xyyyyy[k];

                g_y_xyyyz[k] = -g_0_xyyyz[k] * cd_y[k] + g_0_xyyyyz[k];

                g_y_xyyzz[k] = -g_0_xyyzz[k] * cd_y[k] + g_0_xyyyzz[k];

                g_y_xyzzz[k] = -g_0_xyzzz[k] * cd_y[k] + g_0_xyyzzz[k];

                g_y_xzzzz[k] = -g_0_xzzzz[k] * cd_y[k] + g_0_xyzzzz[k];

                g_y_yyyyy[k] = -g_0_yyyyy[k] * cd_y[k] + g_0_yyyyyy[k];

                g_y_yyyyz[k] = -g_0_yyyyz[k] * cd_y[k] + g_0_yyyyyz[k];

                g_y_yyyzz[k] = -g_0_yyyzz[k] * cd_y[k] + g_0_yyyyzz[k];

                g_y_yyzzz[k] = -g_0_yyzzz[k] * cd_y[k] + g_0_yyyzzz[k];

                g_y_yzzzz[k] = -g_0_yzzzz[k] * cd_y[k] + g_0_yyzzzz[k];

                g_y_zzzzz[k] = -g_0_zzzzz[k] * cd_y[k] + g_0_yzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : contr_buffer_xxph

            auto g_z_xxxxx = contr_buffer_xxph[ph_off + 42];

            auto g_z_xxxxy = contr_buffer_xxph[ph_off + 43];

            auto g_z_xxxxz = contr_buffer_xxph[ph_off + 44];

            auto g_z_xxxyy = contr_buffer_xxph[ph_off + 45];

            auto g_z_xxxyz = contr_buffer_xxph[ph_off + 46];

            auto g_z_xxxzz = contr_buffer_xxph[ph_off + 47];

            auto g_z_xxyyy = contr_buffer_xxph[ph_off + 48];

            auto g_z_xxyyz = contr_buffer_xxph[ph_off + 49];

            auto g_z_xxyzz = contr_buffer_xxph[ph_off + 50];

            auto g_z_xxzzz = contr_buffer_xxph[ph_off + 51];

            auto g_z_xyyyy = contr_buffer_xxph[ph_off + 52];

            auto g_z_xyyyz = contr_buffer_xxph[ph_off + 53];

            auto g_z_xyyzz = contr_buffer_xxph[ph_off + 54];

            auto g_z_xyzzz = contr_buffer_xxph[ph_off + 55];

            auto g_z_xzzzz = contr_buffer_xxph[ph_off + 56];

            auto g_z_yyyyy = contr_buffer_xxph[ph_off + 57];

            auto g_z_yyyyz = contr_buffer_xxph[ph_off + 58];

            auto g_z_yyyzz = contr_buffer_xxph[ph_off + 59];

            auto g_z_yyzzz = contr_buffer_xxph[ph_off + 60];

            auto g_z_yzzzz = contr_buffer_xxph[ph_off + 61];

            auto g_z_zzzzz = contr_buffer_xxph[ph_off + 62];

            #pragma omp simd aligned(cd_z, g_0_xxxxx, g_0_xxxxxz, g_0_xxxxy, g_0_xxxxyz, g_0_xxxxz, g_0_xxxxzz, g_0_xxxyy, g_0_xxxyyz, g_0_xxxyz, g_0_xxxyzz, g_0_xxxzz, g_0_xxxzzz, g_0_xxyyy, g_0_xxyyyz, g_0_xxyyz, g_0_xxyyzz, g_0_xxyzz, g_0_xxyzzz, g_0_xxzzz, g_0_xxzzzz, g_0_xyyyy, g_0_xyyyyz, g_0_xyyyz, g_0_xyyyzz, g_0_xyyzz, g_0_xyyzzz, g_0_xyzzz, g_0_xyzzzz, g_0_xzzzz, g_0_xzzzzz, g_0_yyyyy, g_0_yyyyyz, g_0_yyyyz, g_0_yyyyzz, g_0_yyyzz, g_0_yyyzzz, g_0_yyzzz, g_0_yyzzzz, g_0_yzzzz, g_0_yzzzzz, g_0_zzzzz, g_0_zzzzzz, g_z_xxxxx, g_z_xxxxy, g_z_xxxxz, g_z_xxxyy, g_z_xxxyz, g_z_xxxzz, g_z_xxyyy, g_z_xxyyz, g_z_xxyzz, g_z_xxzzz, g_z_xyyyy, g_z_xyyyz, g_z_xyyzz, g_z_xyzzz, g_z_xzzzz, g_z_yyyyy, g_z_yyyyz, g_z_yyyzz, g_z_yyzzz, g_z_yzzzz, g_z_zzzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_z_xxxxx[k] = -g_0_xxxxx[k] * cd_z[k] + g_0_xxxxxz[k];

                g_z_xxxxy[k] = -g_0_xxxxy[k] * cd_z[k] + g_0_xxxxyz[k];

                g_z_xxxxz[k] = -g_0_xxxxz[k] * cd_z[k] + g_0_xxxxzz[k];

                g_z_xxxyy[k] = -g_0_xxxyy[k] * cd_z[k] + g_0_xxxyyz[k];

                g_z_xxxyz[k] = -g_0_xxxyz[k] * cd_z[k] + g_0_xxxyzz[k];

                g_z_xxxzz[k] = -g_0_xxxzz[k] * cd_z[k] + g_0_xxxzzz[k];

                g_z_xxyyy[k] = -g_0_xxyyy[k] * cd_z[k] + g_0_xxyyyz[k];

                g_z_xxyyz[k] = -g_0_xxyyz[k] * cd_z[k] + g_0_xxyyzz[k];

                g_z_xxyzz[k] = -g_0_xxyzz[k] * cd_z[k] + g_0_xxyzzz[k];

                g_z_xxzzz[k] = -g_0_xxzzz[k] * cd_z[k] + g_0_xxzzzz[k];

                g_z_xyyyy[k] = -g_0_xyyyy[k] * cd_z[k] + g_0_xyyyyz[k];

                g_z_xyyyz[k] = -g_0_xyyyz[k] * cd_z[k] + g_0_xyyyzz[k];

                g_z_xyyzz[k] = -g_0_xyyzz[k] * cd_z[k] + g_0_xyyzzz[k];

                g_z_xyzzz[k] = -g_0_xyzzz[k] * cd_z[k] + g_0_xyzzzz[k];

                g_z_xzzzz[k] = -g_0_xzzzz[k] * cd_z[k] + g_0_xzzzzz[k];

                g_z_yyyyy[k] = -g_0_yyyyy[k] * cd_z[k] + g_0_yyyyyz[k];

                g_z_yyyyz[k] = -g_0_yyyyz[k] * cd_z[k] + g_0_yyyyzz[k];

                g_z_yyyzz[k] = -g_0_yyyzz[k] * cd_z[k] + g_0_yyyzzz[k];

                g_z_yyzzz[k] = -g_0_yyzzz[k] * cd_z[k] + g_0_yyzzzz[k];

                g_z_yzzzz[k] = -g_0_yzzzz[k] * cd_z[k] + g_0_yzzzzz[k];

                g_z_zzzzz[k] = -g_0_zzzzz[k] * cd_z[k] + g_0_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace

