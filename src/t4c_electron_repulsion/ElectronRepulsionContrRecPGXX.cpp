#include "ElectronRepulsionContrRecPGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_hrr_electron_repulsion_pgxx(CSimdArray<double>& contr_buffer_pgxx,
                                     const CSimdArray<double>& contr_buffer_sgxx,
                                     const CSimdArray<double>& contr_buffer_shxx,
                                     const double ab_x,
                                     const double ab_y,
                                     const double ab_z,
                                     const int c_angmom,
                                     const int d_angmom) -> void
{
    const auto ndims = contr_buffer_pgxx.number_of_columns();

    const auto ccomps = tensor::number_of_spherical_components(c_angmom);

    const auto dcomps = tensor::number_of_spherical_components(d_angmom);

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : contr_buffer_sgxx

            const auto sg_off = i * dcomps + j;

            auto g_0_xxxx = contr_buffer_sgxx[sg_off + 0 * ccomps * dcomps];

            auto g_0_xxxy = contr_buffer_sgxx[sg_off + 1 * ccomps * dcomps];

            auto g_0_xxxz = contr_buffer_sgxx[sg_off + 2 * ccomps * dcomps];

            auto g_0_xxyy = contr_buffer_sgxx[sg_off + 3 * ccomps * dcomps];

            auto g_0_xxyz = contr_buffer_sgxx[sg_off + 4 * ccomps * dcomps];

            auto g_0_xxzz = contr_buffer_sgxx[sg_off + 5 * ccomps * dcomps];

            auto g_0_xyyy = contr_buffer_sgxx[sg_off + 6 * ccomps * dcomps];

            auto g_0_xyyz = contr_buffer_sgxx[sg_off + 7 * ccomps * dcomps];

            auto g_0_xyzz = contr_buffer_sgxx[sg_off + 8 * ccomps * dcomps];

            auto g_0_xzzz = contr_buffer_sgxx[sg_off + 9 * ccomps * dcomps];

            auto g_0_yyyy = contr_buffer_sgxx[sg_off + 10 * ccomps * dcomps];

            auto g_0_yyyz = contr_buffer_sgxx[sg_off + 11 * ccomps * dcomps];

            auto g_0_yyzz = contr_buffer_sgxx[sg_off + 12 * ccomps * dcomps];

            auto g_0_yzzz = contr_buffer_sgxx[sg_off + 13 * ccomps * dcomps];

            auto g_0_zzzz = contr_buffer_sgxx[sg_off + 14 * ccomps * dcomps];

            /// Set up components of auxilary buffer : contr_buffer_shxx

            const auto sh_off = i * dcomps + j;

            auto g_0_xxxxx = contr_buffer_shxx[sh_off + 0 * ccomps * dcomps];

            auto g_0_xxxxy = contr_buffer_shxx[sh_off + 1 * ccomps * dcomps];

            auto g_0_xxxxz = contr_buffer_shxx[sh_off + 2 * ccomps * dcomps];

            auto g_0_xxxyy = contr_buffer_shxx[sh_off + 3 * ccomps * dcomps];

            auto g_0_xxxyz = contr_buffer_shxx[sh_off + 4 * ccomps * dcomps];

            auto g_0_xxxzz = contr_buffer_shxx[sh_off + 5 * ccomps * dcomps];

            auto g_0_xxyyy = contr_buffer_shxx[sh_off + 6 * ccomps * dcomps];

            auto g_0_xxyyz = contr_buffer_shxx[sh_off + 7 * ccomps * dcomps];

            auto g_0_xxyzz = contr_buffer_shxx[sh_off + 8 * ccomps * dcomps];

            auto g_0_xxzzz = contr_buffer_shxx[sh_off + 9 * ccomps * dcomps];

            auto g_0_xyyyy = contr_buffer_shxx[sh_off + 10 * ccomps * dcomps];

            auto g_0_xyyyz = contr_buffer_shxx[sh_off + 11 * ccomps * dcomps];

            auto g_0_xyyzz = contr_buffer_shxx[sh_off + 12 * ccomps * dcomps];

            auto g_0_xyzzz = contr_buffer_shxx[sh_off + 13 * ccomps * dcomps];

            auto g_0_xzzzz = contr_buffer_shxx[sh_off + 14 * ccomps * dcomps];

            auto g_0_yyyyy = contr_buffer_shxx[sh_off + 15 * ccomps * dcomps];

            auto g_0_yyyyz = contr_buffer_shxx[sh_off + 16 * ccomps * dcomps];

            auto g_0_yyyzz = contr_buffer_shxx[sh_off + 17 * ccomps * dcomps];

            auto g_0_yyzzz = contr_buffer_shxx[sh_off + 18 * ccomps * dcomps];

            auto g_0_yzzzz = contr_buffer_shxx[sh_off + 19 * ccomps * dcomps];

            auto g_0_zzzzz = contr_buffer_shxx[sh_off + 20 * ccomps * dcomps];

            /// set up bra offset for contr_buffer_pgxx

            const auto pg_off = i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : contr_buffer_pgxx

            auto g_x_xxxx = contr_buffer_pgxx[pg_off + 0 * ccomps * dcomps];

            auto g_x_xxxy = contr_buffer_pgxx[pg_off + 1 * ccomps * dcomps];

            auto g_x_xxxz = contr_buffer_pgxx[pg_off + 2 * ccomps * dcomps];

            auto g_x_xxyy = contr_buffer_pgxx[pg_off + 3 * ccomps * dcomps];

            auto g_x_xxyz = contr_buffer_pgxx[pg_off + 4 * ccomps * dcomps];

            auto g_x_xxzz = contr_buffer_pgxx[pg_off + 5 * ccomps * dcomps];

            auto g_x_xyyy = contr_buffer_pgxx[pg_off + 6 * ccomps * dcomps];

            auto g_x_xyyz = contr_buffer_pgxx[pg_off + 7 * ccomps * dcomps];

            auto g_x_xyzz = contr_buffer_pgxx[pg_off + 8 * ccomps * dcomps];

            auto g_x_xzzz = contr_buffer_pgxx[pg_off + 9 * ccomps * dcomps];

            auto g_x_yyyy = contr_buffer_pgxx[pg_off + 10 * ccomps * dcomps];

            auto g_x_yyyz = contr_buffer_pgxx[pg_off + 11 * ccomps * dcomps];

            auto g_x_yyzz = contr_buffer_pgxx[pg_off + 12 * ccomps * dcomps];

            auto g_x_yzzz = contr_buffer_pgxx[pg_off + 13 * ccomps * dcomps];

            auto g_x_zzzz = contr_buffer_pgxx[pg_off + 14 * ccomps * dcomps];

            #pragma omp simd aligned(g_0_xxxx, g_0_xxxxx, g_0_xxxxy, g_0_xxxxz, g_0_xxxy, g_0_xxxyy, g_0_xxxyz, g_0_xxxz, g_0_xxxzz, g_0_xxyy, g_0_xxyyy, g_0_xxyyz, g_0_xxyz, g_0_xxyzz, g_0_xxzz, g_0_xxzzz, g_0_xyyy, g_0_xyyyy, g_0_xyyyz, g_0_xyyz, g_0_xyyzz, g_0_xyzz, g_0_xyzzz, g_0_xzzz, g_0_xzzzz, g_0_yyyy, g_0_yyyz, g_0_yyzz, g_0_yzzz, g_0_zzzz, g_x_xxxx, g_x_xxxy, g_x_xxxz, g_x_xxyy, g_x_xxyz, g_x_xxzz, g_x_xyyy, g_x_xyyz, g_x_xyzz, g_x_xzzz, g_x_yyyy, g_x_yyyz, g_x_yyzz, g_x_yzzz, g_x_zzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_x_xxxx[k] = -g_0_xxxx[k] * ab_x + g_0_xxxxx[k];

                g_x_xxxy[k] = -g_0_xxxy[k] * ab_x + g_0_xxxxy[k];

                g_x_xxxz[k] = -g_0_xxxz[k] * ab_x + g_0_xxxxz[k];

                g_x_xxyy[k] = -g_0_xxyy[k] * ab_x + g_0_xxxyy[k];

                g_x_xxyz[k] = -g_0_xxyz[k] * ab_x + g_0_xxxyz[k];

                g_x_xxzz[k] = -g_0_xxzz[k] * ab_x + g_0_xxxzz[k];

                g_x_xyyy[k] = -g_0_xyyy[k] * ab_x + g_0_xxyyy[k];

                g_x_xyyz[k] = -g_0_xyyz[k] * ab_x + g_0_xxyyz[k];

                g_x_xyzz[k] = -g_0_xyzz[k] * ab_x + g_0_xxyzz[k];

                g_x_xzzz[k] = -g_0_xzzz[k] * ab_x + g_0_xxzzz[k];

                g_x_yyyy[k] = -g_0_yyyy[k] * ab_x + g_0_xyyyy[k];

                g_x_yyyz[k] = -g_0_yyyz[k] * ab_x + g_0_xyyyz[k];

                g_x_yyzz[k] = -g_0_yyzz[k] * ab_x + g_0_xyyzz[k];

                g_x_yzzz[k] = -g_0_yzzz[k] * ab_x + g_0_xyzzz[k];

                g_x_zzzz[k] = -g_0_zzzz[k] * ab_x + g_0_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : contr_buffer_pgxx

            auto g_y_xxxx = contr_buffer_pgxx[pg_off + 15 * ccomps * dcomps];

            auto g_y_xxxy = contr_buffer_pgxx[pg_off + 16 * ccomps * dcomps];

            auto g_y_xxxz = contr_buffer_pgxx[pg_off + 17 * ccomps * dcomps];

            auto g_y_xxyy = contr_buffer_pgxx[pg_off + 18 * ccomps * dcomps];

            auto g_y_xxyz = contr_buffer_pgxx[pg_off + 19 * ccomps * dcomps];

            auto g_y_xxzz = contr_buffer_pgxx[pg_off + 20 * ccomps * dcomps];

            auto g_y_xyyy = contr_buffer_pgxx[pg_off + 21 * ccomps * dcomps];

            auto g_y_xyyz = contr_buffer_pgxx[pg_off + 22 * ccomps * dcomps];

            auto g_y_xyzz = contr_buffer_pgxx[pg_off + 23 * ccomps * dcomps];

            auto g_y_xzzz = contr_buffer_pgxx[pg_off + 24 * ccomps * dcomps];

            auto g_y_yyyy = contr_buffer_pgxx[pg_off + 25 * ccomps * dcomps];

            auto g_y_yyyz = contr_buffer_pgxx[pg_off + 26 * ccomps * dcomps];

            auto g_y_yyzz = contr_buffer_pgxx[pg_off + 27 * ccomps * dcomps];

            auto g_y_yzzz = contr_buffer_pgxx[pg_off + 28 * ccomps * dcomps];

            auto g_y_zzzz = contr_buffer_pgxx[pg_off + 29 * ccomps * dcomps];

            #pragma omp simd aligned(g_0_xxxx, g_0_xxxxy, g_0_xxxy, g_0_xxxyy, g_0_xxxyz, g_0_xxxz, g_0_xxyy, g_0_xxyyy, g_0_xxyyz, g_0_xxyz, g_0_xxyzz, g_0_xxzz, g_0_xyyy, g_0_xyyyy, g_0_xyyyz, g_0_xyyz, g_0_xyyzz, g_0_xyzz, g_0_xyzzz, g_0_xzzz, g_0_yyyy, g_0_yyyyy, g_0_yyyyz, g_0_yyyz, g_0_yyyzz, g_0_yyzz, g_0_yyzzz, g_0_yzzz, g_0_yzzzz, g_0_zzzz, g_y_xxxx, g_y_xxxy, g_y_xxxz, g_y_xxyy, g_y_xxyz, g_y_xxzz, g_y_xyyy, g_y_xyyz, g_y_xyzz, g_y_xzzz, g_y_yyyy, g_y_yyyz, g_y_yyzz, g_y_yzzz, g_y_zzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_y_xxxx[k] = -g_0_xxxx[k] * ab_y + g_0_xxxxy[k];

                g_y_xxxy[k] = -g_0_xxxy[k] * ab_y + g_0_xxxyy[k];

                g_y_xxxz[k] = -g_0_xxxz[k] * ab_y + g_0_xxxyz[k];

                g_y_xxyy[k] = -g_0_xxyy[k] * ab_y + g_0_xxyyy[k];

                g_y_xxyz[k] = -g_0_xxyz[k] * ab_y + g_0_xxyyz[k];

                g_y_xxzz[k] = -g_0_xxzz[k] * ab_y + g_0_xxyzz[k];

                g_y_xyyy[k] = -g_0_xyyy[k] * ab_y + g_0_xyyyy[k];

                g_y_xyyz[k] = -g_0_xyyz[k] * ab_y + g_0_xyyyz[k];

                g_y_xyzz[k] = -g_0_xyzz[k] * ab_y + g_0_xyyzz[k];

                g_y_xzzz[k] = -g_0_xzzz[k] * ab_y + g_0_xyzzz[k];

                g_y_yyyy[k] = -g_0_yyyy[k] * ab_y + g_0_yyyyy[k];

                g_y_yyyz[k] = -g_0_yyyz[k] * ab_y + g_0_yyyyz[k];

                g_y_yyzz[k] = -g_0_yyzz[k] * ab_y + g_0_yyyzz[k];

                g_y_yzzz[k] = -g_0_yzzz[k] * ab_y + g_0_yyzzz[k];

                g_y_zzzz[k] = -g_0_zzzz[k] * ab_y + g_0_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : contr_buffer_pgxx

            auto g_z_xxxx = contr_buffer_pgxx[pg_off + 30 * ccomps * dcomps];

            auto g_z_xxxy = contr_buffer_pgxx[pg_off + 31 * ccomps * dcomps];

            auto g_z_xxxz = contr_buffer_pgxx[pg_off + 32 * ccomps * dcomps];

            auto g_z_xxyy = contr_buffer_pgxx[pg_off + 33 * ccomps * dcomps];

            auto g_z_xxyz = contr_buffer_pgxx[pg_off + 34 * ccomps * dcomps];

            auto g_z_xxzz = contr_buffer_pgxx[pg_off + 35 * ccomps * dcomps];

            auto g_z_xyyy = contr_buffer_pgxx[pg_off + 36 * ccomps * dcomps];

            auto g_z_xyyz = contr_buffer_pgxx[pg_off + 37 * ccomps * dcomps];

            auto g_z_xyzz = contr_buffer_pgxx[pg_off + 38 * ccomps * dcomps];

            auto g_z_xzzz = contr_buffer_pgxx[pg_off + 39 * ccomps * dcomps];

            auto g_z_yyyy = contr_buffer_pgxx[pg_off + 40 * ccomps * dcomps];

            auto g_z_yyyz = contr_buffer_pgxx[pg_off + 41 * ccomps * dcomps];

            auto g_z_yyzz = contr_buffer_pgxx[pg_off + 42 * ccomps * dcomps];

            auto g_z_yzzz = contr_buffer_pgxx[pg_off + 43 * ccomps * dcomps];

            auto g_z_zzzz = contr_buffer_pgxx[pg_off + 44 * ccomps * dcomps];

            #pragma omp simd aligned(g_0_xxxx, g_0_xxxxz, g_0_xxxy, g_0_xxxyz, g_0_xxxz, g_0_xxxzz, g_0_xxyy, g_0_xxyyz, g_0_xxyz, g_0_xxyzz, g_0_xxzz, g_0_xxzzz, g_0_xyyy, g_0_xyyyz, g_0_xyyz, g_0_xyyzz, g_0_xyzz, g_0_xyzzz, g_0_xzzz, g_0_xzzzz, g_0_yyyy, g_0_yyyyz, g_0_yyyz, g_0_yyyzz, g_0_yyzz, g_0_yyzzz, g_0_yzzz, g_0_yzzzz, g_0_zzzz, g_0_zzzzz, g_z_xxxx, g_z_xxxy, g_z_xxxz, g_z_xxyy, g_z_xxyz, g_z_xxzz, g_z_xyyy, g_z_xyyz, g_z_xyzz, g_z_xzzz, g_z_yyyy, g_z_yyyz, g_z_yyzz, g_z_yzzz, g_z_zzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_z_xxxx[k] = -g_0_xxxx[k] * ab_z + g_0_xxxxz[k];

                g_z_xxxy[k] = -g_0_xxxy[k] * ab_z + g_0_xxxyz[k];

                g_z_xxxz[k] = -g_0_xxxz[k] * ab_z + g_0_xxxzz[k];

                g_z_xxyy[k] = -g_0_xxyy[k] * ab_z + g_0_xxyyz[k];

                g_z_xxyz[k] = -g_0_xxyz[k] * ab_z + g_0_xxyzz[k];

                g_z_xxzz[k] = -g_0_xxzz[k] * ab_z + g_0_xxzzz[k];

                g_z_xyyy[k] = -g_0_xyyy[k] * ab_z + g_0_xyyyz[k];

                g_z_xyyz[k] = -g_0_xyyz[k] * ab_z + g_0_xyyzz[k];

                g_z_xyzz[k] = -g_0_xyzz[k] * ab_z + g_0_xyzzz[k];

                g_z_xzzz[k] = -g_0_xzzz[k] * ab_z + g_0_xzzzz[k];

                g_z_yyyy[k] = -g_0_yyyy[k] * ab_z + g_0_yyyyz[k];

                g_z_yyyz[k] = -g_0_yyyz[k] * ab_z + g_0_yyyzz[k];

                g_z_yyzz[k] = -g_0_yyzz[k] * ab_z + g_0_yyzzz[k];

                g_z_yzzz[k] = -g_0_yzzz[k] * ab_z + g_0_yzzzz[k];

                g_z_zzzz[k] = -g_0_zzzz[k] * ab_z + g_0_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

