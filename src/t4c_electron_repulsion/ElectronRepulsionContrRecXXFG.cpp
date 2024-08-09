#include "ElectronRepulsionContrRecXXFG.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_hrr_electron_repulsion_xxfg(CSimdArray<double>& contr_buffer_xxfg,
                                     const CSimdArray<double>& contr_buffer_xxdg,
                                     const CSimdArray<double>& contr_buffer_xxdh,
                                     const double* cd_x,
                                     const double* cd_y,
                                     const double* cd_z,
                                     const int a_angmom,
                                     const int b_angmom) -> void
{
    const auto ndims = contr_buffer_xxfg.number_of_columns();

    const auto acomps = tensor::number_of_cartesian_components(a_angmom);

    const auto bcomps = tensor::number_of_cartesian_components(b_angmom);

    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
            /// Set up components of auxilary buffer : contr_buffer_xxdg

            const auto dg_off = (i * bcomps + j) * 90;

            auto g_xx_xxxx = contr_buffer_xxdg[dg_off + 0];

            auto g_xx_xxxy = contr_buffer_xxdg[dg_off + 1];

            auto g_xx_xxxz = contr_buffer_xxdg[dg_off + 2];

            auto g_xx_xxyy = contr_buffer_xxdg[dg_off + 3];

            auto g_xx_xxyz = contr_buffer_xxdg[dg_off + 4];

            auto g_xx_xxzz = contr_buffer_xxdg[dg_off + 5];

            auto g_xx_xyyy = contr_buffer_xxdg[dg_off + 6];

            auto g_xx_xyyz = contr_buffer_xxdg[dg_off + 7];

            auto g_xx_xyzz = contr_buffer_xxdg[dg_off + 8];

            auto g_xx_xzzz = contr_buffer_xxdg[dg_off + 9];

            auto g_xx_yyyy = contr_buffer_xxdg[dg_off + 10];

            auto g_xx_yyyz = contr_buffer_xxdg[dg_off + 11];

            auto g_xx_yyzz = contr_buffer_xxdg[dg_off + 12];

            auto g_xx_yzzz = contr_buffer_xxdg[dg_off + 13];

            auto g_xx_zzzz = contr_buffer_xxdg[dg_off + 14];

            auto g_xy_xxxx = contr_buffer_xxdg[dg_off + 15];

            auto g_xy_xxxy = contr_buffer_xxdg[dg_off + 16];

            auto g_xy_xxxz = contr_buffer_xxdg[dg_off + 17];

            auto g_xy_xxyy = contr_buffer_xxdg[dg_off + 18];

            auto g_xy_xxyz = contr_buffer_xxdg[dg_off + 19];

            auto g_xy_xxzz = contr_buffer_xxdg[dg_off + 20];

            auto g_xy_xyyy = contr_buffer_xxdg[dg_off + 21];

            auto g_xy_xyyz = contr_buffer_xxdg[dg_off + 22];

            auto g_xy_xyzz = contr_buffer_xxdg[dg_off + 23];

            auto g_xy_xzzz = contr_buffer_xxdg[dg_off + 24];

            auto g_xy_yyyy = contr_buffer_xxdg[dg_off + 25];

            auto g_xy_yyyz = contr_buffer_xxdg[dg_off + 26];

            auto g_xy_yyzz = contr_buffer_xxdg[dg_off + 27];

            auto g_xy_yzzz = contr_buffer_xxdg[dg_off + 28];

            auto g_xy_zzzz = contr_buffer_xxdg[dg_off + 29];

            auto g_xz_xxxx = contr_buffer_xxdg[dg_off + 30];

            auto g_xz_xxxy = contr_buffer_xxdg[dg_off + 31];

            auto g_xz_xxxz = contr_buffer_xxdg[dg_off + 32];

            auto g_xz_xxyy = contr_buffer_xxdg[dg_off + 33];

            auto g_xz_xxyz = contr_buffer_xxdg[dg_off + 34];

            auto g_xz_xxzz = contr_buffer_xxdg[dg_off + 35];

            auto g_xz_xyyy = contr_buffer_xxdg[dg_off + 36];

            auto g_xz_xyyz = contr_buffer_xxdg[dg_off + 37];

            auto g_xz_xyzz = contr_buffer_xxdg[dg_off + 38];

            auto g_xz_xzzz = contr_buffer_xxdg[dg_off + 39];

            auto g_xz_yyyy = contr_buffer_xxdg[dg_off + 40];

            auto g_xz_yyyz = contr_buffer_xxdg[dg_off + 41];

            auto g_xz_yyzz = contr_buffer_xxdg[dg_off + 42];

            auto g_xz_yzzz = contr_buffer_xxdg[dg_off + 43];

            auto g_xz_zzzz = contr_buffer_xxdg[dg_off + 44];

            auto g_yy_xxxx = contr_buffer_xxdg[dg_off + 45];

            auto g_yy_xxxy = contr_buffer_xxdg[dg_off + 46];

            auto g_yy_xxxz = contr_buffer_xxdg[dg_off + 47];

            auto g_yy_xxyy = contr_buffer_xxdg[dg_off + 48];

            auto g_yy_xxyz = contr_buffer_xxdg[dg_off + 49];

            auto g_yy_xxzz = contr_buffer_xxdg[dg_off + 50];

            auto g_yy_xyyy = contr_buffer_xxdg[dg_off + 51];

            auto g_yy_xyyz = contr_buffer_xxdg[dg_off + 52];

            auto g_yy_xyzz = contr_buffer_xxdg[dg_off + 53];

            auto g_yy_xzzz = contr_buffer_xxdg[dg_off + 54];

            auto g_yy_yyyy = contr_buffer_xxdg[dg_off + 55];

            auto g_yy_yyyz = contr_buffer_xxdg[dg_off + 56];

            auto g_yy_yyzz = contr_buffer_xxdg[dg_off + 57];

            auto g_yy_yzzz = contr_buffer_xxdg[dg_off + 58];

            auto g_yy_zzzz = contr_buffer_xxdg[dg_off + 59];

            auto g_yz_xxxx = contr_buffer_xxdg[dg_off + 60];

            auto g_yz_xxxy = contr_buffer_xxdg[dg_off + 61];

            auto g_yz_xxxz = contr_buffer_xxdg[dg_off + 62];

            auto g_yz_xxyy = contr_buffer_xxdg[dg_off + 63];

            auto g_yz_xxyz = contr_buffer_xxdg[dg_off + 64];

            auto g_yz_xxzz = contr_buffer_xxdg[dg_off + 65];

            auto g_yz_xyyy = contr_buffer_xxdg[dg_off + 66];

            auto g_yz_xyyz = contr_buffer_xxdg[dg_off + 67];

            auto g_yz_xyzz = contr_buffer_xxdg[dg_off + 68];

            auto g_yz_xzzz = contr_buffer_xxdg[dg_off + 69];

            auto g_yz_yyyy = contr_buffer_xxdg[dg_off + 70];

            auto g_yz_yyyz = contr_buffer_xxdg[dg_off + 71];

            auto g_yz_yyzz = contr_buffer_xxdg[dg_off + 72];

            auto g_yz_yzzz = contr_buffer_xxdg[dg_off + 73];

            auto g_yz_zzzz = contr_buffer_xxdg[dg_off + 74];

            auto g_zz_xxxx = contr_buffer_xxdg[dg_off + 75];

            auto g_zz_xxxy = contr_buffer_xxdg[dg_off + 76];

            auto g_zz_xxxz = contr_buffer_xxdg[dg_off + 77];

            auto g_zz_xxyy = contr_buffer_xxdg[dg_off + 78];

            auto g_zz_xxyz = contr_buffer_xxdg[dg_off + 79];

            auto g_zz_xxzz = contr_buffer_xxdg[dg_off + 80];

            auto g_zz_xyyy = contr_buffer_xxdg[dg_off + 81];

            auto g_zz_xyyz = contr_buffer_xxdg[dg_off + 82];

            auto g_zz_xyzz = contr_buffer_xxdg[dg_off + 83];

            auto g_zz_xzzz = contr_buffer_xxdg[dg_off + 84];

            auto g_zz_yyyy = contr_buffer_xxdg[dg_off + 85];

            auto g_zz_yyyz = contr_buffer_xxdg[dg_off + 86];

            auto g_zz_yyzz = contr_buffer_xxdg[dg_off + 87];

            auto g_zz_yzzz = contr_buffer_xxdg[dg_off + 88];

            auto g_zz_zzzz = contr_buffer_xxdg[dg_off + 89];

            /// Set up components of auxilary buffer : contr_buffer_xxdh

            const auto dh_off = (i * bcomps + j) * 126;

            auto g_xx_xxxxx = contr_buffer_xxdh[dh_off + 0];

            auto g_xx_xxxxy = contr_buffer_xxdh[dh_off + 1];

            auto g_xx_xxxxz = contr_buffer_xxdh[dh_off + 2];

            auto g_xx_xxxyy = contr_buffer_xxdh[dh_off + 3];

            auto g_xx_xxxyz = contr_buffer_xxdh[dh_off + 4];

            auto g_xx_xxxzz = contr_buffer_xxdh[dh_off + 5];

            auto g_xx_xxyyy = contr_buffer_xxdh[dh_off + 6];

            auto g_xx_xxyyz = contr_buffer_xxdh[dh_off + 7];

            auto g_xx_xxyzz = contr_buffer_xxdh[dh_off + 8];

            auto g_xx_xxzzz = contr_buffer_xxdh[dh_off + 9];

            auto g_xx_xyyyy = contr_buffer_xxdh[dh_off + 10];

            auto g_xx_xyyyz = contr_buffer_xxdh[dh_off + 11];

            auto g_xx_xyyzz = contr_buffer_xxdh[dh_off + 12];

            auto g_xx_xyzzz = contr_buffer_xxdh[dh_off + 13];

            auto g_xx_xzzzz = contr_buffer_xxdh[dh_off + 14];

            auto g_xy_xxxxx = contr_buffer_xxdh[dh_off + 21];

            auto g_xy_xxxxy = contr_buffer_xxdh[dh_off + 22];

            auto g_xy_xxxxz = contr_buffer_xxdh[dh_off + 23];

            auto g_xy_xxxyy = contr_buffer_xxdh[dh_off + 24];

            auto g_xy_xxxyz = contr_buffer_xxdh[dh_off + 25];

            auto g_xy_xxxzz = contr_buffer_xxdh[dh_off + 26];

            auto g_xy_xxyyy = contr_buffer_xxdh[dh_off + 27];

            auto g_xy_xxyyz = contr_buffer_xxdh[dh_off + 28];

            auto g_xy_xxyzz = contr_buffer_xxdh[dh_off + 29];

            auto g_xy_xxzzz = contr_buffer_xxdh[dh_off + 30];

            auto g_xy_xyyyy = contr_buffer_xxdh[dh_off + 31];

            auto g_xy_xyyyz = contr_buffer_xxdh[dh_off + 32];

            auto g_xy_xyyzz = contr_buffer_xxdh[dh_off + 33];

            auto g_xy_xyzzz = contr_buffer_xxdh[dh_off + 34];

            auto g_xy_xzzzz = contr_buffer_xxdh[dh_off + 35];

            auto g_xz_xxxxx = contr_buffer_xxdh[dh_off + 42];

            auto g_xz_xxxxy = contr_buffer_xxdh[dh_off + 43];

            auto g_xz_xxxxz = contr_buffer_xxdh[dh_off + 44];

            auto g_xz_xxxyy = contr_buffer_xxdh[dh_off + 45];

            auto g_xz_xxxyz = contr_buffer_xxdh[dh_off + 46];

            auto g_xz_xxxzz = contr_buffer_xxdh[dh_off + 47];

            auto g_xz_xxyyy = contr_buffer_xxdh[dh_off + 48];

            auto g_xz_xxyyz = contr_buffer_xxdh[dh_off + 49];

            auto g_xz_xxyzz = contr_buffer_xxdh[dh_off + 50];

            auto g_xz_xxzzz = contr_buffer_xxdh[dh_off + 51];

            auto g_xz_xyyyy = contr_buffer_xxdh[dh_off + 52];

            auto g_xz_xyyyz = contr_buffer_xxdh[dh_off + 53];

            auto g_xz_xyyzz = contr_buffer_xxdh[dh_off + 54];

            auto g_xz_xyzzz = contr_buffer_xxdh[dh_off + 55];

            auto g_xz_xzzzz = contr_buffer_xxdh[dh_off + 56];

            auto g_yy_xxxxx = contr_buffer_xxdh[dh_off + 63];

            auto g_yy_xxxxy = contr_buffer_xxdh[dh_off + 64];

            auto g_yy_xxxxz = contr_buffer_xxdh[dh_off + 65];

            auto g_yy_xxxyy = contr_buffer_xxdh[dh_off + 66];

            auto g_yy_xxxyz = contr_buffer_xxdh[dh_off + 67];

            auto g_yy_xxxzz = contr_buffer_xxdh[dh_off + 68];

            auto g_yy_xxyyy = contr_buffer_xxdh[dh_off + 69];

            auto g_yy_xxyyz = contr_buffer_xxdh[dh_off + 70];

            auto g_yy_xxyzz = contr_buffer_xxdh[dh_off + 71];

            auto g_yy_xxzzz = contr_buffer_xxdh[dh_off + 72];

            auto g_yy_xyyyy = contr_buffer_xxdh[dh_off + 73];

            auto g_yy_xyyyz = contr_buffer_xxdh[dh_off + 74];

            auto g_yy_xyyzz = contr_buffer_xxdh[dh_off + 75];

            auto g_yy_xyzzz = contr_buffer_xxdh[dh_off + 76];

            auto g_yy_xzzzz = contr_buffer_xxdh[dh_off + 77];

            auto g_yy_yyyyy = contr_buffer_xxdh[dh_off + 78];

            auto g_yy_yyyyz = contr_buffer_xxdh[dh_off + 79];

            auto g_yy_yyyzz = contr_buffer_xxdh[dh_off + 80];

            auto g_yy_yyzzz = contr_buffer_xxdh[dh_off + 81];

            auto g_yy_yzzzz = contr_buffer_xxdh[dh_off + 82];

            auto g_yz_xxxxx = contr_buffer_xxdh[dh_off + 84];

            auto g_yz_xxxxy = contr_buffer_xxdh[dh_off + 85];

            auto g_yz_xxxxz = contr_buffer_xxdh[dh_off + 86];

            auto g_yz_xxxyy = contr_buffer_xxdh[dh_off + 87];

            auto g_yz_xxxyz = contr_buffer_xxdh[dh_off + 88];

            auto g_yz_xxxzz = contr_buffer_xxdh[dh_off + 89];

            auto g_yz_xxyyy = contr_buffer_xxdh[dh_off + 90];

            auto g_yz_xxyyz = contr_buffer_xxdh[dh_off + 91];

            auto g_yz_xxyzz = contr_buffer_xxdh[dh_off + 92];

            auto g_yz_xxzzz = contr_buffer_xxdh[dh_off + 93];

            auto g_yz_xyyyy = contr_buffer_xxdh[dh_off + 94];

            auto g_yz_xyyyz = contr_buffer_xxdh[dh_off + 95];

            auto g_yz_xyyzz = contr_buffer_xxdh[dh_off + 96];

            auto g_yz_xyzzz = contr_buffer_xxdh[dh_off + 97];

            auto g_yz_xzzzz = contr_buffer_xxdh[dh_off + 98];

            auto g_yz_yyyyy = contr_buffer_xxdh[dh_off + 99];

            auto g_yz_yyyyz = contr_buffer_xxdh[dh_off + 100];

            auto g_yz_yyyzz = contr_buffer_xxdh[dh_off + 101];

            auto g_yz_yyzzz = contr_buffer_xxdh[dh_off + 102];

            auto g_yz_yzzzz = contr_buffer_xxdh[dh_off + 103];

            auto g_zz_xxxxx = contr_buffer_xxdh[dh_off + 105];

            auto g_zz_xxxxy = contr_buffer_xxdh[dh_off + 106];

            auto g_zz_xxxxz = contr_buffer_xxdh[dh_off + 107];

            auto g_zz_xxxyy = contr_buffer_xxdh[dh_off + 108];

            auto g_zz_xxxyz = contr_buffer_xxdh[dh_off + 109];

            auto g_zz_xxxzz = contr_buffer_xxdh[dh_off + 110];

            auto g_zz_xxyyy = contr_buffer_xxdh[dh_off + 111];

            auto g_zz_xxyyz = contr_buffer_xxdh[dh_off + 112];

            auto g_zz_xxyzz = contr_buffer_xxdh[dh_off + 113];

            auto g_zz_xxzzz = contr_buffer_xxdh[dh_off + 114];

            auto g_zz_xyyyy = contr_buffer_xxdh[dh_off + 115];

            auto g_zz_xyyyz = contr_buffer_xxdh[dh_off + 116];

            auto g_zz_xyyzz = contr_buffer_xxdh[dh_off + 117];

            auto g_zz_xyzzz = contr_buffer_xxdh[dh_off + 118];

            auto g_zz_xzzzz = contr_buffer_xxdh[dh_off + 119];

            auto g_zz_yyyyy = contr_buffer_xxdh[dh_off + 120];

            auto g_zz_yyyyz = contr_buffer_xxdh[dh_off + 121];

            auto g_zz_yyyzz = contr_buffer_xxdh[dh_off + 122];

            auto g_zz_yyzzz = contr_buffer_xxdh[dh_off + 123];

            auto g_zz_yzzzz = contr_buffer_xxdh[dh_off + 124];

            auto g_zz_zzzzz = contr_buffer_xxdh[dh_off + 125];

            /// set up bra offset for contr_buffer_xxfg

            const auto fg_off = (i * bcomps + j) * 150;

            /// Set up 0-15 components of targeted buffer : contr_buffer_xxfg

            auto g_xxx_xxxx = contr_buffer_xxfg[fg_off + 0];

            auto g_xxx_xxxy = contr_buffer_xxfg[fg_off + 1];

            auto g_xxx_xxxz = contr_buffer_xxfg[fg_off + 2];

            auto g_xxx_xxyy = contr_buffer_xxfg[fg_off + 3];

            auto g_xxx_xxyz = contr_buffer_xxfg[fg_off + 4];

            auto g_xxx_xxzz = contr_buffer_xxfg[fg_off + 5];

            auto g_xxx_xyyy = contr_buffer_xxfg[fg_off + 6];

            auto g_xxx_xyyz = contr_buffer_xxfg[fg_off + 7];

            auto g_xxx_xyzz = contr_buffer_xxfg[fg_off + 8];

            auto g_xxx_xzzz = contr_buffer_xxfg[fg_off + 9];

            auto g_xxx_yyyy = contr_buffer_xxfg[fg_off + 10];

            auto g_xxx_yyyz = contr_buffer_xxfg[fg_off + 11];

            auto g_xxx_yyzz = contr_buffer_xxfg[fg_off + 12];

            auto g_xxx_yzzz = contr_buffer_xxfg[fg_off + 13];

            auto g_xxx_zzzz = contr_buffer_xxfg[fg_off + 14];

            #pragma omp simd aligned(cd_x, g_xx_xxxx, g_xx_xxxxx, g_xx_xxxxy, g_xx_xxxxz, g_xx_xxxy, g_xx_xxxyy, g_xx_xxxyz, g_xx_xxxz, g_xx_xxxzz, g_xx_xxyy, g_xx_xxyyy, g_xx_xxyyz, g_xx_xxyz, g_xx_xxyzz, g_xx_xxzz, g_xx_xxzzz, g_xx_xyyy, g_xx_xyyyy, g_xx_xyyyz, g_xx_xyyz, g_xx_xyyzz, g_xx_xyzz, g_xx_xyzzz, g_xx_xzzz, g_xx_xzzzz, g_xx_yyyy, g_xx_yyyz, g_xx_yyzz, g_xx_yzzz, g_xx_zzzz, g_xxx_xxxx, g_xxx_xxxy, g_xxx_xxxz, g_xxx_xxyy, g_xxx_xxyz, g_xxx_xxzz, g_xxx_xyyy, g_xxx_xyyz, g_xxx_xyzz, g_xxx_xzzz, g_xxx_yyyy, g_xxx_yyyz, g_xxx_yyzz, g_xxx_yzzz, g_xxx_zzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_xxx_xxxx[k] = -g_xx_xxxx[k] * cd_x[k] + g_xx_xxxxx[k];

                g_xxx_xxxy[k] = -g_xx_xxxy[k] * cd_x[k] + g_xx_xxxxy[k];

                g_xxx_xxxz[k] = -g_xx_xxxz[k] * cd_x[k] + g_xx_xxxxz[k];

                g_xxx_xxyy[k] = -g_xx_xxyy[k] * cd_x[k] + g_xx_xxxyy[k];

                g_xxx_xxyz[k] = -g_xx_xxyz[k] * cd_x[k] + g_xx_xxxyz[k];

                g_xxx_xxzz[k] = -g_xx_xxzz[k] * cd_x[k] + g_xx_xxxzz[k];

                g_xxx_xyyy[k] = -g_xx_xyyy[k] * cd_x[k] + g_xx_xxyyy[k];

                g_xxx_xyyz[k] = -g_xx_xyyz[k] * cd_x[k] + g_xx_xxyyz[k];

                g_xxx_xyzz[k] = -g_xx_xyzz[k] * cd_x[k] + g_xx_xxyzz[k];

                g_xxx_xzzz[k] = -g_xx_xzzz[k] * cd_x[k] + g_xx_xxzzz[k];

                g_xxx_yyyy[k] = -g_xx_yyyy[k] * cd_x[k] + g_xx_xyyyy[k];

                g_xxx_yyyz[k] = -g_xx_yyyz[k] * cd_x[k] + g_xx_xyyyz[k];

                g_xxx_yyzz[k] = -g_xx_yyzz[k] * cd_x[k] + g_xx_xyyzz[k];

                g_xxx_yzzz[k] = -g_xx_yzzz[k] * cd_x[k] + g_xx_xyzzz[k];

                g_xxx_zzzz[k] = -g_xx_zzzz[k] * cd_x[k] + g_xx_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : contr_buffer_xxfg

            auto g_xxy_xxxx = contr_buffer_xxfg[fg_off + 15];

            auto g_xxy_xxxy = contr_buffer_xxfg[fg_off + 16];

            auto g_xxy_xxxz = contr_buffer_xxfg[fg_off + 17];

            auto g_xxy_xxyy = contr_buffer_xxfg[fg_off + 18];

            auto g_xxy_xxyz = contr_buffer_xxfg[fg_off + 19];

            auto g_xxy_xxzz = contr_buffer_xxfg[fg_off + 20];

            auto g_xxy_xyyy = contr_buffer_xxfg[fg_off + 21];

            auto g_xxy_xyyz = contr_buffer_xxfg[fg_off + 22];

            auto g_xxy_xyzz = contr_buffer_xxfg[fg_off + 23];

            auto g_xxy_xzzz = contr_buffer_xxfg[fg_off + 24];

            auto g_xxy_yyyy = contr_buffer_xxfg[fg_off + 25];

            auto g_xxy_yyyz = contr_buffer_xxfg[fg_off + 26];

            auto g_xxy_yyzz = contr_buffer_xxfg[fg_off + 27];

            auto g_xxy_yzzz = contr_buffer_xxfg[fg_off + 28];

            auto g_xxy_zzzz = contr_buffer_xxfg[fg_off + 29];

            #pragma omp simd aligned(cd_x, g_xxy_xxxx, g_xxy_xxxy, g_xxy_xxxz, g_xxy_xxyy, g_xxy_xxyz, g_xxy_xxzz, g_xxy_xyyy, g_xxy_xyyz, g_xxy_xyzz, g_xxy_xzzz, g_xxy_yyyy, g_xxy_yyyz, g_xxy_yyzz, g_xxy_yzzz, g_xxy_zzzz, g_xy_xxxx, g_xy_xxxxx, g_xy_xxxxy, g_xy_xxxxz, g_xy_xxxy, g_xy_xxxyy, g_xy_xxxyz, g_xy_xxxz, g_xy_xxxzz, g_xy_xxyy, g_xy_xxyyy, g_xy_xxyyz, g_xy_xxyz, g_xy_xxyzz, g_xy_xxzz, g_xy_xxzzz, g_xy_xyyy, g_xy_xyyyy, g_xy_xyyyz, g_xy_xyyz, g_xy_xyyzz, g_xy_xyzz, g_xy_xyzzz, g_xy_xzzz, g_xy_xzzzz, g_xy_yyyy, g_xy_yyyz, g_xy_yyzz, g_xy_yzzz, g_xy_zzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_xxy_xxxx[k] = -g_xy_xxxx[k] * cd_x[k] + g_xy_xxxxx[k];

                g_xxy_xxxy[k] = -g_xy_xxxy[k] * cd_x[k] + g_xy_xxxxy[k];

                g_xxy_xxxz[k] = -g_xy_xxxz[k] * cd_x[k] + g_xy_xxxxz[k];

                g_xxy_xxyy[k] = -g_xy_xxyy[k] * cd_x[k] + g_xy_xxxyy[k];

                g_xxy_xxyz[k] = -g_xy_xxyz[k] * cd_x[k] + g_xy_xxxyz[k];

                g_xxy_xxzz[k] = -g_xy_xxzz[k] * cd_x[k] + g_xy_xxxzz[k];

                g_xxy_xyyy[k] = -g_xy_xyyy[k] * cd_x[k] + g_xy_xxyyy[k];

                g_xxy_xyyz[k] = -g_xy_xyyz[k] * cd_x[k] + g_xy_xxyyz[k];

                g_xxy_xyzz[k] = -g_xy_xyzz[k] * cd_x[k] + g_xy_xxyzz[k];

                g_xxy_xzzz[k] = -g_xy_xzzz[k] * cd_x[k] + g_xy_xxzzz[k];

                g_xxy_yyyy[k] = -g_xy_yyyy[k] * cd_x[k] + g_xy_xyyyy[k];

                g_xxy_yyyz[k] = -g_xy_yyyz[k] * cd_x[k] + g_xy_xyyyz[k];

                g_xxy_yyzz[k] = -g_xy_yyzz[k] * cd_x[k] + g_xy_xyyzz[k];

                g_xxy_yzzz[k] = -g_xy_yzzz[k] * cd_x[k] + g_xy_xyzzz[k];

                g_xxy_zzzz[k] = -g_xy_zzzz[k] * cd_x[k] + g_xy_xzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : contr_buffer_xxfg

            auto g_xxz_xxxx = contr_buffer_xxfg[fg_off + 30];

            auto g_xxz_xxxy = contr_buffer_xxfg[fg_off + 31];

            auto g_xxz_xxxz = contr_buffer_xxfg[fg_off + 32];

            auto g_xxz_xxyy = contr_buffer_xxfg[fg_off + 33];

            auto g_xxz_xxyz = contr_buffer_xxfg[fg_off + 34];

            auto g_xxz_xxzz = contr_buffer_xxfg[fg_off + 35];

            auto g_xxz_xyyy = contr_buffer_xxfg[fg_off + 36];

            auto g_xxz_xyyz = contr_buffer_xxfg[fg_off + 37];

            auto g_xxz_xyzz = contr_buffer_xxfg[fg_off + 38];

            auto g_xxz_xzzz = contr_buffer_xxfg[fg_off + 39];

            auto g_xxz_yyyy = contr_buffer_xxfg[fg_off + 40];

            auto g_xxz_yyyz = contr_buffer_xxfg[fg_off + 41];

            auto g_xxz_yyzz = contr_buffer_xxfg[fg_off + 42];

            auto g_xxz_yzzz = contr_buffer_xxfg[fg_off + 43];

            auto g_xxz_zzzz = contr_buffer_xxfg[fg_off + 44];

            #pragma omp simd aligned(cd_x, g_xxz_xxxx, g_xxz_xxxy, g_xxz_xxxz, g_xxz_xxyy, g_xxz_xxyz, g_xxz_xxzz, g_xxz_xyyy, g_xxz_xyyz, g_xxz_xyzz, g_xxz_xzzz, g_xxz_yyyy, g_xxz_yyyz, g_xxz_yyzz, g_xxz_yzzz, g_xxz_zzzz, g_xz_xxxx, g_xz_xxxxx, g_xz_xxxxy, g_xz_xxxxz, g_xz_xxxy, g_xz_xxxyy, g_xz_xxxyz, g_xz_xxxz, g_xz_xxxzz, g_xz_xxyy, g_xz_xxyyy, g_xz_xxyyz, g_xz_xxyz, g_xz_xxyzz, g_xz_xxzz, g_xz_xxzzz, g_xz_xyyy, g_xz_xyyyy, g_xz_xyyyz, g_xz_xyyz, g_xz_xyyzz, g_xz_xyzz, g_xz_xyzzz, g_xz_xzzz, g_xz_xzzzz, g_xz_yyyy, g_xz_yyyz, g_xz_yyzz, g_xz_yzzz, g_xz_zzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_xxz_xxxx[k] = -g_xz_xxxx[k] * cd_x[k] + g_xz_xxxxx[k];

                g_xxz_xxxy[k] = -g_xz_xxxy[k] * cd_x[k] + g_xz_xxxxy[k];

                g_xxz_xxxz[k] = -g_xz_xxxz[k] * cd_x[k] + g_xz_xxxxz[k];

                g_xxz_xxyy[k] = -g_xz_xxyy[k] * cd_x[k] + g_xz_xxxyy[k];

                g_xxz_xxyz[k] = -g_xz_xxyz[k] * cd_x[k] + g_xz_xxxyz[k];

                g_xxz_xxzz[k] = -g_xz_xxzz[k] * cd_x[k] + g_xz_xxxzz[k];

                g_xxz_xyyy[k] = -g_xz_xyyy[k] * cd_x[k] + g_xz_xxyyy[k];

                g_xxz_xyyz[k] = -g_xz_xyyz[k] * cd_x[k] + g_xz_xxyyz[k];

                g_xxz_xyzz[k] = -g_xz_xyzz[k] * cd_x[k] + g_xz_xxyzz[k];

                g_xxz_xzzz[k] = -g_xz_xzzz[k] * cd_x[k] + g_xz_xxzzz[k];

                g_xxz_yyyy[k] = -g_xz_yyyy[k] * cd_x[k] + g_xz_xyyyy[k];

                g_xxz_yyyz[k] = -g_xz_yyyz[k] * cd_x[k] + g_xz_xyyyz[k];

                g_xxz_yyzz[k] = -g_xz_yyzz[k] * cd_x[k] + g_xz_xyyzz[k];

                g_xxz_yzzz[k] = -g_xz_yzzz[k] * cd_x[k] + g_xz_xyzzz[k];

                g_xxz_zzzz[k] = -g_xz_zzzz[k] * cd_x[k] + g_xz_xzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : contr_buffer_xxfg

            auto g_xyy_xxxx = contr_buffer_xxfg[fg_off + 45];

            auto g_xyy_xxxy = contr_buffer_xxfg[fg_off + 46];

            auto g_xyy_xxxz = contr_buffer_xxfg[fg_off + 47];

            auto g_xyy_xxyy = contr_buffer_xxfg[fg_off + 48];

            auto g_xyy_xxyz = contr_buffer_xxfg[fg_off + 49];

            auto g_xyy_xxzz = contr_buffer_xxfg[fg_off + 50];

            auto g_xyy_xyyy = contr_buffer_xxfg[fg_off + 51];

            auto g_xyy_xyyz = contr_buffer_xxfg[fg_off + 52];

            auto g_xyy_xyzz = contr_buffer_xxfg[fg_off + 53];

            auto g_xyy_xzzz = contr_buffer_xxfg[fg_off + 54];

            auto g_xyy_yyyy = contr_buffer_xxfg[fg_off + 55];

            auto g_xyy_yyyz = contr_buffer_xxfg[fg_off + 56];

            auto g_xyy_yyzz = contr_buffer_xxfg[fg_off + 57];

            auto g_xyy_yzzz = contr_buffer_xxfg[fg_off + 58];

            auto g_xyy_zzzz = contr_buffer_xxfg[fg_off + 59];

            #pragma omp simd aligned(cd_x, g_xyy_xxxx, g_xyy_xxxy, g_xyy_xxxz, g_xyy_xxyy, g_xyy_xxyz, g_xyy_xxzz, g_xyy_xyyy, g_xyy_xyyz, g_xyy_xyzz, g_xyy_xzzz, g_xyy_yyyy, g_xyy_yyyz, g_xyy_yyzz, g_xyy_yzzz, g_xyy_zzzz, g_yy_xxxx, g_yy_xxxxx, g_yy_xxxxy, g_yy_xxxxz, g_yy_xxxy, g_yy_xxxyy, g_yy_xxxyz, g_yy_xxxz, g_yy_xxxzz, g_yy_xxyy, g_yy_xxyyy, g_yy_xxyyz, g_yy_xxyz, g_yy_xxyzz, g_yy_xxzz, g_yy_xxzzz, g_yy_xyyy, g_yy_xyyyy, g_yy_xyyyz, g_yy_xyyz, g_yy_xyyzz, g_yy_xyzz, g_yy_xyzzz, g_yy_xzzz, g_yy_xzzzz, g_yy_yyyy, g_yy_yyyz, g_yy_yyzz, g_yy_yzzz, g_yy_zzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_xyy_xxxx[k] = -g_yy_xxxx[k] * cd_x[k] + g_yy_xxxxx[k];

                g_xyy_xxxy[k] = -g_yy_xxxy[k] * cd_x[k] + g_yy_xxxxy[k];

                g_xyy_xxxz[k] = -g_yy_xxxz[k] * cd_x[k] + g_yy_xxxxz[k];

                g_xyy_xxyy[k] = -g_yy_xxyy[k] * cd_x[k] + g_yy_xxxyy[k];

                g_xyy_xxyz[k] = -g_yy_xxyz[k] * cd_x[k] + g_yy_xxxyz[k];

                g_xyy_xxzz[k] = -g_yy_xxzz[k] * cd_x[k] + g_yy_xxxzz[k];

                g_xyy_xyyy[k] = -g_yy_xyyy[k] * cd_x[k] + g_yy_xxyyy[k];

                g_xyy_xyyz[k] = -g_yy_xyyz[k] * cd_x[k] + g_yy_xxyyz[k];

                g_xyy_xyzz[k] = -g_yy_xyzz[k] * cd_x[k] + g_yy_xxyzz[k];

                g_xyy_xzzz[k] = -g_yy_xzzz[k] * cd_x[k] + g_yy_xxzzz[k];

                g_xyy_yyyy[k] = -g_yy_yyyy[k] * cd_x[k] + g_yy_xyyyy[k];

                g_xyy_yyyz[k] = -g_yy_yyyz[k] * cd_x[k] + g_yy_xyyyz[k];

                g_xyy_yyzz[k] = -g_yy_yyzz[k] * cd_x[k] + g_yy_xyyzz[k];

                g_xyy_yzzz[k] = -g_yy_yzzz[k] * cd_x[k] + g_yy_xyzzz[k];

                g_xyy_zzzz[k] = -g_yy_zzzz[k] * cd_x[k] + g_yy_xzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : contr_buffer_xxfg

            auto g_xyz_xxxx = contr_buffer_xxfg[fg_off + 60];

            auto g_xyz_xxxy = contr_buffer_xxfg[fg_off + 61];

            auto g_xyz_xxxz = contr_buffer_xxfg[fg_off + 62];

            auto g_xyz_xxyy = contr_buffer_xxfg[fg_off + 63];

            auto g_xyz_xxyz = contr_buffer_xxfg[fg_off + 64];

            auto g_xyz_xxzz = contr_buffer_xxfg[fg_off + 65];

            auto g_xyz_xyyy = contr_buffer_xxfg[fg_off + 66];

            auto g_xyz_xyyz = contr_buffer_xxfg[fg_off + 67];

            auto g_xyz_xyzz = contr_buffer_xxfg[fg_off + 68];

            auto g_xyz_xzzz = contr_buffer_xxfg[fg_off + 69];

            auto g_xyz_yyyy = contr_buffer_xxfg[fg_off + 70];

            auto g_xyz_yyyz = contr_buffer_xxfg[fg_off + 71];

            auto g_xyz_yyzz = contr_buffer_xxfg[fg_off + 72];

            auto g_xyz_yzzz = contr_buffer_xxfg[fg_off + 73];

            auto g_xyz_zzzz = contr_buffer_xxfg[fg_off + 74];

            #pragma omp simd aligned(cd_x, g_xyz_xxxx, g_xyz_xxxy, g_xyz_xxxz, g_xyz_xxyy, g_xyz_xxyz, g_xyz_xxzz, g_xyz_xyyy, g_xyz_xyyz, g_xyz_xyzz, g_xyz_xzzz, g_xyz_yyyy, g_xyz_yyyz, g_xyz_yyzz, g_xyz_yzzz, g_xyz_zzzz, g_yz_xxxx, g_yz_xxxxx, g_yz_xxxxy, g_yz_xxxxz, g_yz_xxxy, g_yz_xxxyy, g_yz_xxxyz, g_yz_xxxz, g_yz_xxxzz, g_yz_xxyy, g_yz_xxyyy, g_yz_xxyyz, g_yz_xxyz, g_yz_xxyzz, g_yz_xxzz, g_yz_xxzzz, g_yz_xyyy, g_yz_xyyyy, g_yz_xyyyz, g_yz_xyyz, g_yz_xyyzz, g_yz_xyzz, g_yz_xyzzz, g_yz_xzzz, g_yz_xzzzz, g_yz_yyyy, g_yz_yyyz, g_yz_yyzz, g_yz_yzzz, g_yz_zzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_xyz_xxxx[k] = -g_yz_xxxx[k] * cd_x[k] + g_yz_xxxxx[k];

                g_xyz_xxxy[k] = -g_yz_xxxy[k] * cd_x[k] + g_yz_xxxxy[k];

                g_xyz_xxxz[k] = -g_yz_xxxz[k] * cd_x[k] + g_yz_xxxxz[k];

                g_xyz_xxyy[k] = -g_yz_xxyy[k] * cd_x[k] + g_yz_xxxyy[k];

                g_xyz_xxyz[k] = -g_yz_xxyz[k] * cd_x[k] + g_yz_xxxyz[k];

                g_xyz_xxzz[k] = -g_yz_xxzz[k] * cd_x[k] + g_yz_xxxzz[k];

                g_xyz_xyyy[k] = -g_yz_xyyy[k] * cd_x[k] + g_yz_xxyyy[k];

                g_xyz_xyyz[k] = -g_yz_xyyz[k] * cd_x[k] + g_yz_xxyyz[k];

                g_xyz_xyzz[k] = -g_yz_xyzz[k] * cd_x[k] + g_yz_xxyzz[k];

                g_xyz_xzzz[k] = -g_yz_xzzz[k] * cd_x[k] + g_yz_xxzzz[k];

                g_xyz_yyyy[k] = -g_yz_yyyy[k] * cd_x[k] + g_yz_xyyyy[k];

                g_xyz_yyyz[k] = -g_yz_yyyz[k] * cd_x[k] + g_yz_xyyyz[k];

                g_xyz_yyzz[k] = -g_yz_yyzz[k] * cd_x[k] + g_yz_xyyzz[k];

                g_xyz_yzzz[k] = -g_yz_yzzz[k] * cd_x[k] + g_yz_xyzzz[k];

                g_xyz_zzzz[k] = -g_yz_zzzz[k] * cd_x[k] + g_yz_xzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : contr_buffer_xxfg

            auto g_xzz_xxxx = contr_buffer_xxfg[fg_off + 75];

            auto g_xzz_xxxy = contr_buffer_xxfg[fg_off + 76];

            auto g_xzz_xxxz = contr_buffer_xxfg[fg_off + 77];

            auto g_xzz_xxyy = contr_buffer_xxfg[fg_off + 78];

            auto g_xzz_xxyz = contr_buffer_xxfg[fg_off + 79];

            auto g_xzz_xxzz = contr_buffer_xxfg[fg_off + 80];

            auto g_xzz_xyyy = contr_buffer_xxfg[fg_off + 81];

            auto g_xzz_xyyz = contr_buffer_xxfg[fg_off + 82];

            auto g_xzz_xyzz = contr_buffer_xxfg[fg_off + 83];

            auto g_xzz_xzzz = contr_buffer_xxfg[fg_off + 84];

            auto g_xzz_yyyy = contr_buffer_xxfg[fg_off + 85];

            auto g_xzz_yyyz = contr_buffer_xxfg[fg_off + 86];

            auto g_xzz_yyzz = contr_buffer_xxfg[fg_off + 87];

            auto g_xzz_yzzz = contr_buffer_xxfg[fg_off + 88];

            auto g_xzz_zzzz = contr_buffer_xxfg[fg_off + 89];

            #pragma omp simd aligned(cd_x, g_xzz_xxxx, g_xzz_xxxy, g_xzz_xxxz, g_xzz_xxyy, g_xzz_xxyz, g_xzz_xxzz, g_xzz_xyyy, g_xzz_xyyz, g_xzz_xyzz, g_xzz_xzzz, g_xzz_yyyy, g_xzz_yyyz, g_xzz_yyzz, g_xzz_yzzz, g_xzz_zzzz, g_zz_xxxx, g_zz_xxxxx, g_zz_xxxxy, g_zz_xxxxz, g_zz_xxxy, g_zz_xxxyy, g_zz_xxxyz, g_zz_xxxz, g_zz_xxxzz, g_zz_xxyy, g_zz_xxyyy, g_zz_xxyyz, g_zz_xxyz, g_zz_xxyzz, g_zz_xxzz, g_zz_xxzzz, g_zz_xyyy, g_zz_xyyyy, g_zz_xyyyz, g_zz_xyyz, g_zz_xyyzz, g_zz_xyzz, g_zz_xyzzz, g_zz_xzzz, g_zz_xzzzz, g_zz_yyyy, g_zz_yyyz, g_zz_yyzz, g_zz_yzzz, g_zz_zzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_xzz_xxxx[k] = -g_zz_xxxx[k] * cd_x[k] + g_zz_xxxxx[k];

                g_xzz_xxxy[k] = -g_zz_xxxy[k] * cd_x[k] + g_zz_xxxxy[k];

                g_xzz_xxxz[k] = -g_zz_xxxz[k] * cd_x[k] + g_zz_xxxxz[k];

                g_xzz_xxyy[k] = -g_zz_xxyy[k] * cd_x[k] + g_zz_xxxyy[k];

                g_xzz_xxyz[k] = -g_zz_xxyz[k] * cd_x[k] + g_zz_xxxyz[k];

                g_xzz_xxzz[k] = -g_zz_xxzz[k] * cd_x[k] + g_zz_xxxzz[k];

                g_xzz_xyyy[k] = -g_zz_xyyy[k] * cd_x[k] + g_zz_xxyyy[k];

                g_xzz_xyyz[k] = -g_zz_xyyz[k] * cd_x[k] + g_zz_xxyyz[k];

                g_xzz_xyzz[k] = -g_zz_xyzz[k] * cd_x[k] + g_zz_xxyzz[k];

                g_xzz_xzzz[k] = -g_zz_xzzz[k] * cd_x[k] + g_zz_xxzzz[k];

                g_xzz_yyyy[k] = -g_zz_yyyy[k] * cd_x[k] + g_zz_xyyyy[k];

                g_xzz_yyyz[k] = -g_zz_yyyz[k] * cd_x[k] + g_zz_xyyyz[k];

                g_xzz_yyzz[k] = -g_zz_yyzz[k] * cd_x[k] + g_zz_xyyzz[k];

                g_xzz_yzzz[k] = -g_zz_yzzz[k] * cd_x[k] + g_zz_xyzzz[k];

                g_xzz_zzzz[k] = -g_zz_zzzz[k] * cd_x[k] + g_zz_xzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : contr_buffer_xxfg

            auto g_yyy_xxxx = contr_buffer_xxfg[fg_off + 90];

            auto g_yyy_xxxy = contr_buffer_xxfg[fg_off + 91];

            auto g_yyy_xxxz = contr_buffer_xxfg[fg_off + 92];

            auto g_yyy_xxyy = contr_buffer_xxfg[fg_off + 93];

            auto g_yyy_xxyz = contr_buffer_xxfg[fg_off + 94];

            auto g_yyy_xxzz = contr_buffer_xxfg[fg_off + 95];

            auto g_yyy_xyyy = contr_buffer_xxfg[fg_off + 96];

            auto g_yyy_xyyz = contr_buffer_xxfg[fg_off + 97];

            auto g_yyy_xyzz = contr_buffer_xxfg[fg_off + 98];

            auto g_yyy_xzzz = contr_buffer_xxfg[fg_off + 99];

            auto g_yyy_yyyy = contr_buffer_xxfg[fg_off + 100];

            auto g_yyy_yyyz = contr_buffer_xxfg[fg_off + 101];

            auto g_yyy_yyzz = contr_buffer_xxfg[fg_off + 102];

            auto g_yyy_yzzz = contr_buffer_xxfg[fg_off + 103];

            auto g_yyy_zzzz = contr_buffer_xxfg[fg_off + 104];

            #pragma omp simd aligned(cd_y, g_yy_xxxx, g_yy_xxxxy, g_yy_xxxy, g_yy_xxxyy, g_yy_xxxyz, g_yy_xxxz, g_yy_xxyy, g_yy_xxyyy, g_yy_xxyyz, g_yy_xxyz, g_yy_xxyzz, g_yy_xxzz, g_yy_xyyy, g_yy_xyyyy, g_yy_xyyyz, g_yy_xyyz, g_yy_xyyzz, g_yy_xyzz, g_yy_xyzzz, g_yy_xzzz, g_yy_yyyy, g_yy_yyyyy, g_yy_yyyyz, g_yy_yyyz, g_yy_yyyzz, g_yy_yyzz, g_yy_yyzzz, g_yy_yzzz, g_yy_yzzzz, g_yy_zzzz, g_yyy_xxxx, g_yyy_xxxy, g_yyy_xxxz, g_yyy_xxyy, g_yyy_xxyz, g_yyy_xxzz, g_yyy_xyyy, g_yyy_xyyz, g_yyy_xyzz, g_yyy_xzzz, g_yyy_yyyy, g_yyy_yyyz, g_yyy_yyzz, g_yyy_yzzz, g_yyy_zzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_yyy_xxxx[k] = -g_yy_xxxx[k] * cd_y[k] + g_yy_xxxxy[k];

                g_yyy_xxxy[k] = -g_yy_xxxy[k] * cd_y[k] + g_yy_xxxyy[k];

                g_yyy_xxxz[k] = -g_yy_xxxz[k] * cd_y[k] + g_yy_xxxyz[k];

                g_yyy_xxyy[k] = -g_yy_xxyy[k] * cd_y[k] + g_yy_xxyyy[k];

                g_yyy_xxyz[k] = -g_yy_xxyz[k] * cd_y[k] + g_yy_xxyyz[k];

                g_yyy_xxzz[k] = -g_yy_xxzz[k] * cd_y[k] + g_yy_xxyzz[k];

                g_yyy_xyyy[k] = -g_yy_xyyy[k] * cd_y[k] + g_yy_xyyyy[k];

                g_yyy_xyyz[k] = -g_yy_xyyz[k] * cd_y[k] + g_yy_xyyyz[k];

                g_yyy_xyzz[k] = -g_yy_xyzz[k] * cd_y[k] + g_yy_xyyzz[k];

                g_yyy_xzzz[k] = -g_yy_xzzz[k] * cd_y[k] + g_yy_xyzzz[k];

                g_yyy_yyyy[k] = -g_yy_yyyy[k] * cd_y[k] + g_yy_yyyyy[k];

                g_yyy_yyyz[k] = -g_yy_yyyz[k] * cd_y[k] + g_yy_yyyyz[k];

                g_yyy_yyzz[k] = -g_yy_yyzz[k] * cd_y[k] + g_yy_yyyzz[k];

                g_yyy_yzzz[k] = -g_yy_yzzz[k] * cd_y[k] + g_yy_yyzzz[k];

                g_yyy_zzzz[k] = -g_yy_zzzz[k] * cd_y[k] + g_yy_yzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : contr_buffer_xxfg

            auto g_yyz_xxxx = contr_buffer_xxfg[fg_off + 105];

            auto g_yyz_xxxy = contr_buffer_xxfg[fg_off + 106];

            auto g_yyz_xxxz = contr_buffer_xxfg[fg_off + 107];

            auto g_yyz_xxyy = contr_buffer_xxfg[fg_off + 108];

            auto g_yyz_xxyz = contr_buffer_xxfg[fg_off + 109];

            auto g_yyz_xxzz = contr_buffer_xxfg[fg_off + 110];

            auto g_yyz_xyyy = contr_buffer_xxfg[fg_off + 111];

            auto g_yyz_xyyz = contr_buffer_xxfg[fg_off + 112];

            auto g_yyz_xyzz = contr_buffer_xxfg[fg_off + 113];

            auto g_yyz_xzzz = contr_buffer_xxfg[fg_off + 114];

            auto g_yyz_yyyy = contr_buffer_xxfg[fg_off + 115];

            auto g_yyz_yyyz = contr_buffer_xxfg[fg_off + 116];

            auto g_yyz_yyzz = contr_buffer_xxfg[fg_off + 117];

            auto g_yyz_yzzz = contr_buffer_xxfg[fg_off + 118];

            auto g_yyz_zzzz = contr_buffer_xxfg[fg_off + 119];

            #pragma omp simd aligned(cd_y, g_yyz_xxxx, g_yyz_xxxy, g_yyz_xxxz, g_yyz_xxyy, g_yyz_xxyz, g_yyz_xxzz, g_yyz_xyyy, g_yyz_xyyz, g_yyz_xyzz, g_yyz_xzzz, g_yyz_yyyy, g_yyz_yyyz, g_yyz_yyzz, g_yyz_yzzz, g_yyz_zzzz, g_yz_xxxx, g_yz_xxxxy, g_yz_xxxy, g_yz_xxxyy, g_yz_xxxyz, g_yz_xxxz, g_yz_xxyy, g_yz_xxyyy, g_yz_xxyyz, g_yz_xxyz, g_yz_xxyzz, g_yz_xxzz, g_yz_xyyy, g_yz_xyyyy, g_yz_xyyyz, g_yz_xyyz, g_yz_xyyzz, g_yz_xyzz, g_yz_xyzzz, g_yz_xzzz, g_yz_yyyy, g_yz_yyyyy, g_yz_yyyyz, g_yz_yyyz, g_yz_yyyzz, g_yz_yyzz, g_yz_yyzzz, g_yz_yzzz, g_yz_yzzzz, g_yz_zzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_yyz_xxxx[k] = -g_yz_xxxx[k] * cd_y[k] + g_yz_xxxxy[k];

                g_yyz_xxxy[k] = -g_yz_xxxy[k] * cd_y[k] + g_yz_xxxyy[k];

                g_yyz_xxxz[k] = -g_yz_xxxz[k] * cd_y[k] + g_yz_xxxyz[k];

                g_yyz_xxyy[k] = -g_yz_xxyy[k] * cd_y[k] + g_yz_xxyyy[k];

                g_yyz_xxyz[k] = -g_yz_xxyz[k] * cd_y[k] + g_yz_xxyyz[k];

                g_yyz_xxzz[k] = -g_yz_xxzz[k] * cd_y[k] + g_yz_xxyzz[k];

                g_yyz_xyyy[k] = -g_yz_xyyy[k] * cd_y[k] + g_yz_xyyyy[k];

                g_yyz_xyyz[k] = -g_yz_xyyz[k] * cd_y[k] + g_yz_xyyyz[k];

                g_yyz_xyzz[k] = -g_yz_xyzz[k] * cd_y[k] + g_yz_xyyzz[k];

                g_yyz_xzzz[k] = -g_yz_xzzz[k] * cd_y[k] + g_yz_xyzzz[k];

                g_yyz_yyyy[k] = -g_yz_yyyy[k] * cd_y[k] + g_yz_yyyyy[k];

                g_yyz_yyyz[k] = -g_yz_yyyz[k] * cd_y[k] + g_yz_yyyyz[k];

                g_yyz_yyzz[k] = -g_yz_yyzz[k] * cd_y[k] + g_yz_yyyzz[k];

                g_yyz_yzzz[k] = -g_yz_yzzz[k] * cd_y[k] + g_yz_yyzzz[k];

                g_yyz_zzzz[k] = -g_yz_zzzz[k] * cd_y[k] + g_yz_yzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : contr_buffer_xxfg

            auto g_yzz_xxxx = contr_buffer_xxfg[fg_off + 120];

            auto g_yzz_xxxy = contr_buffer_xxfg[fg_off + 121];

            auto g_yzz_xxxz = contr_buffer_xxfg[fg_off + 122];

            auto g_yzz_xxyy = contr_buffer_xxfg[fg_off + 123];

            auto g_yzz_xxyz = contr_buffer_xxfg[fg_off + 124];

            auto g_yzz_xxzz = contr_buffer_xxfg[fg_off + 125];

            auto g_yzz_xyyy = contr_buffer_xxfg[fg_off + 126];

            auto g_yzz_xyyz = contr_buffer_xxfg[fg_off + 127];

            auto g_yzz_xyzz = contr_buffer_xxfg[fg_off + 128];

            auto g_yzz_xzzz = contr_buffer_xxfg[fg_off + 129];

            auto g_yzz_yyyy = contr_buffer_xxfg[fg_off + 130];

            auto g_yzz_yyyz = contr_buffer_xxfg[fg_off + 131];

            auto g_yzz_yyzz = contr_buffer_xxfg[fg_off + 132];

            auto g_yzz_yzzz = contr_buffer_xxfg[fg_off + 133];

            auto g_yzz_zzzz = contr_buffer_xxfg[fg_off + 134];

            #pragma omp simd aligned(cd_y, g_yzz_xxxx, g_yzz_xxxy, g_yzz_xxxz, g_yzz_xxyy, g_yzz_xxyz, g_yzz_xxzz, g_yzz_xyyy, g_yzz_xyyz, g_yzz_xyzz, g_yzz_xzzz, g_yzz_yyyy, g_yzz_yyyz, g_yzz_yyzz, g_yzz_yzzz, g_yzz_zzzz, g_zz_xxxx, g_zz_xxxxy, g_zz_xxxy, g_zz_xxxyy, g_zz_xxxyz, g_zz_xxxz, g_zz_xxyy, g_zz_xxyyy, g_zz_xxyyz, g_zz_xxyz, g_zz_xxyzz, g_zz_xxzz, g_zz_xyyy, g_zz_xyyyy, g_zz_xyyyz, g_zz_xyyz, g_zz_xyyzz, g_zz_xyzz, g_zz_xyzzz, g_zz_xzzz, g_zz_yyyy, g_zz_yyyyy, g_zz_yyyyz, g_zz_yyyz, g_zz_yyyzz, g_zz_yyzz, g_zz_yyzzz, g_zz_yzzz, g_zz_yzzzz, g_zz_zzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_yzz_xxxx[k] = -g_zz_xxxx[k] * cd_y[k] + g_zz_xxxxy[k];

                g_yzz_xxxy[k] = -g_zz_xxxy[k] * cd_y[k] + g_zz_xxxyy[k];

                g_yzz_xxxz[k] = -g_zz_xxxz[k] * cd_y[k] + g_zz_xxxyz[k];

                g_yzz_xxyy[k] = -g_zz_xxyy[k] * cd_y[k] + g_zz_xxyyy[k];

                g_yzz_xxyz[k] = -g_zz_xxyz[k] * cd_y[k] + g_zz_xxyyz[k];

                g_yzz_xxzz[k] = -g_zz_xxzz[k] * cd_y[k] + g_zz_xxyzz[k];

                g_yzz_xyyy[k] = -g_zz_xyyy[k] * cd_y[k] + g_zz_xyyyy[k];

                g_yzz_xyyz[k] = -g_zz_xyyz[k] * cd_y[k] + g_zz_xyyyz[k];

                g_yzz_xyzz[k] = -g_zz_xyzz[k] * cd_y[k] + g_zz_xyyzz[k];

                g_yzz_xzzz[k] = -g_zz_xzzz[k] * cd_y[k] + g_zz_xyzzz[k];

                g_yzz_yyyy[k] = -g_zz_yyyy[k] * cd_y[k] + g_zz_yyyyy[k];

                g_yzz_yyyz[k] = -g_zz_yyyz[k] * cd_y[k] + g_zz_yyyyz[k];

                g_yzz_yyzz[k] = -g_zz_yyzz[k] * cd_y[k] + g_zz_yyyzz[k];

                g_yzz_yzzz[k] = -g_zz_yzzz[k] * cd_y[k] + g_zz_yyzzz[k];

                g_yzz_zzzz[k] = -g_zz_zzzz[k] * cd_y[k] + g_zz_yzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : contr_buffer_xxfg

            auto g_zzz_xxxx = contr_buffer_xxfg[fg_off + 135];

            auto g_zzz_xxxy = contr_buffer_xxfg[fg_off + 136];

            auto g_zzz_xxxz = contr_buffer_xxfg[fg_off + 137];

            auto g_zzz_xxyy = contr_buffer_xxfg[fg_off + 138];

            auto g_zzz_xxyz = contr_buffer_xxfg[fg_off + 139];

            auto g_zzz_xxzz = contr_buffer_xxfg[fg_off + 140];

            auto g_zzz_xyyy = contr_buffer_xxfg[fg_off + 141];

            auto g_zzz_xyyz = contr_buffer_xxfg[fg_off + 142];

            auto g_zzz_xyzz = contr_buffer_xxfg[fg_off + 143];

            auto g_zzz_xzzz = contr_buffer_xxfg[fg_off + 144];

            auto g_zzz_yyyy = contr_buffer_xxfg[fg_off + 145];

            auto g_zzz_yyyz = contr_buffer_xxfg[fg_off + 146];

            auto g_zzz_yyzz = contr_buffer_xxfg[fg_off + 147];

            auto g_zzz_yzzz = contr_buffer_xxfg[fg_off + 148];

            auto g_zzz_zzzz = contr_buffer_xxfg[fg_off + 149];

            #pragma omp simd aligned(cd_z, g_zz_xxxx, g_zz_xxxxz, g_zz_xxxy, g_zz_xxxyz, g_zz_xxxz, g_zz_xxxzz, g_zz_xxyy, g_zz_xxyyz, g_zz_xxyz, g_zz_xxyzz, g_zz_xxzz, g_zz_xxzzz, g_zz_xyyy, g_zz_xyyyz, g_zz_xyyz, g_zz_xyyzz, g_zz_xyzz, g_zz_xyzzz, g_zz_xzzz, g_zz_xzzzz, g_zz_yyyy, g_zz_yyyyz, g_zz_yyyz, g_zz_yyyzz, g_zz_yyzz, g_zz_yyzzz, g_zz_yzzz, g_zz_yzzzz, g_zz_zzzz, g_zz_zzzzz, g_zzz_xxxx, g_zzz_xxxy, g_zzz_xxxz, g_zzz_xxyy, g_zzz_xxyz, g_zzz_xxzz, g_zzz_xyyy, g_zzz_xyyz, g_zzz_xyzz, g_zzz_xzzz, g_zzz_yyyy, g_zzz_yyyz, g_zzz_yyzz, g_zzz_yzzz, g_zzz_zzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_zzz_xxxx[k] = -g_zz_xxxx[k] * cd_z[k] + g_zz_xxxxz[k];

                g_zzz_xxxy[k] = -g_zz_xxxy[k] * cd_z[k] + g_zz_xxxyz[k];

                g_zzz_xxxz[k] = -g_zz_xxxz[k] * cd_z[k] + g_zz_xxxzz[k];

                g_zzz_xxyy[k] = -g_zz_xxyy[k] * cd_z[k] + g_zz_xxyyz[k];

                g_zzz_xxyz[k] = -g_zz_xxyz[k] * cd_z[k] + g_zz_xxyzz[k];

                g_zzz_xxzz[k] = -g_zz_xxzz[k] * cd_z[k] + g_zz_xxzzz[k];

                g_zzz_xyyy[k] = -g_zz_xyyy[k] * cd_z[k] + g_zz_xyyyz[k];

                g_zzz_xyyz[k] = -g_zz_xyyz[k] * cd_z[k] + g_zz_xyyzz[k];

                g_zzz_xyzz[k] = -g_zz_xyzz[k] * cd_z[k] + g_zz_xyzzz[k];

                g_zzz_xzzz[k] = -g_zz_xzzz[k] * cd_z[k] + g_zz_xzzzz[k];

                g_zzz_yyyy[k] = -g_zz_yyyy[k] * cd_z[k] + g_zz_yyyyz[k];

                g_zzz_yyyz[k] = -g_zz_yyyz[k] * cd_z[k] + g_zz_yyyzz[k];

                g_zzz_yyzz[k] = -g_zz_yyzz[k] * cd_z[k] + g_zz_yyzzz[k];

                g_zzz_yzzz[k] = -g_zz_yzzz[k] * cd_z[k] + g_zz_yzzzz[k];

                g_zzz_zzzz[k] = -g_zz_zzzz[k] * cd_z[k] + g_zz_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

