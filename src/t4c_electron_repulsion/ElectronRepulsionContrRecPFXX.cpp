#include "ElectronRepulsionContrRecPFXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_hrr_electron_repulsion_pfxx(CSimdArray<double>& contr_buffer_pfxx,
                                     const CSimdArray<double>& contr_buffer_sfxx,
                                     const CSimdArray<double>& contr_buffer_sgxx,
                                     const double ab_x,
                                     const double ab_y,
                                     const double ab_z,
                                     const int c_angmom,
                                     const int d_angmom) -> void
{
    const auto ndims = contr_buffer_pfxx.number_of_columns();

    const auto ccomps = tensor::number_of_spherical_components(c_angmom);

    const auto dcomps = tensor::number_of_spherical_components(d_angmom);

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : contr_buffer_sfxx

            const auto sf_off = i * dcomps + j;

            auto g_0_xxx = contr_buffer_sfxx[sf_off + 0 * ccomps * dcomps];

            auto g_0_xxy = contr_buffer_sfxx[sf_off + 1 * ccomps * dcomps];

            auto g_0_xxz = contr_buffer_sfxx[sf_off + 2 * ccomps * dcomps];

            auto g_0_xyy = contr_buffer_sfxx[sf_off + 3 * ccomps * dcomps];

            auto g_0_xyz = contr_buffer_sfxx[sf_off + 4 * ccomps * dcomps];

            auto g_0_xzz = contr_buffer_sfxx[sf_off + 5 * ccomps * dcomps];

            auto g_0_yyy = contr_buffer_sfxx[sf_off + 6 * ccomps * dcomps];

            auto g_0_yyz = contr_buffer_sfxx[sf_off + 7 * ccomps * dcomps];

            auto g_0_yzz = contr_buffer_sfxx[sf_off + 8 * ccomps * dcomps];

            auto g_0_zzz = contr_buffer_sfxx[sf_off + 9 * ccomps * dcomps];

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

            /// set up bra offset for contr_buffer_pfxx

            const auto pf_off = i * dcomps + j;

            /// Set up 0-10 components of targeted buffer : contr_buffer_pfxx

            auto g_x_xxx = contr_buffer_pfxx[pf_off + 0 * ccomps * dcomps];

            auto g_x_xxy = contr_buffer_pfxx[pf_off + 1 * ccomps * dcomps];

            auto g_x_xxz = contr_buffer_pfxx[pf_off + 2 * ccomps * dcomps];

            auto g_x_xyy = contr_buffer_pfxx[pf_off + 3 * ccomps * dcomps];

            auto g_x_xyz = contr_buffer_pfxx[pf_off + 4 * ccomps * dcomps];

            auto g_x_xzz = contr_buffer_pfxx[pf_off + 5 * ccomps * dcomps];

            auto g_x_yyy = contr_buffer_pfxx[pf_off + 6 * ccomps * dcomps];

            auto g_x_yyz = contr_buffer_pfxx[pf_off + 7 * ccomps * dcomps];

            auto g_x_yzz = contr_buffer_pfxx[pf_off + 8 * ccomps * dcomps];

            auto g_x_zzz = contr_buffer_pfxx[pf_off + 9 * ccomps * dcomps];

            #pragma omp simd aligned(g_0_xxx, g_0_xxxx, g_0_xxxy, g_0_xxxz, g_0_xxy, g_0_xxyy, g_0_xxyz, g_0_xxz, g_0_xxzz, g_0_xyy, g_0_xyyy, g_0_xyyz, g_0_xyz, g_0_xyzz, g_0_xzz, g_0_xzzz, g_0_yyy, g_0_yyz, g_0_yzz, g_0_zzz, g_x_xxx, g_x_xxy, g_x_xxz, g_x_xyy, g_x_xyz, g_x_xzz, g_x_yyy, g_x_yyz, g_x_yzz, g_x_zzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_x_xxx[k] = -g_0_xxx[k] * ab_x + g_0_xxxx[k];

                g_x_xxy[k] = -g_0_xxy[k] * ab_x + g_0_xxxy[k];

                g_x_xxz[k] = -g_0_xxz[k] * ab_x + g_0_xxxz[k];

                g_x_xyy[k] = -g_0_xyy[k] * ab_x + g_0_xxyy[k];

                g_x_xyz[k] = -g_0_xyz[k] * ab_x + g_0_xxyz[k];

                g_x_xzz[k] = -g_0_xzz[k] * ab_x + g_0_xxzz[k];

                g_x_yyy[k] = -g_0_yyy[k] * ab_x + g_0_xyyy[k];

                g_x_yyz[k] = -g_0_yyz[k] * ab_x + g_0_xyyz[k];

                g_x_yzz[k] = -g_0_yzz[k] * ab_x + g_0_xyzz[k];

                g_x_zzz[k] = -g_0_zzz[k] * ab_x + g_0_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : contr_buffer_pfxx

            auto g_y_xxx = contr_buffer_pfxx[pf_off + 10 * ccomps * dcomps];

            auto g_y_xxy = contr_buffer_pfxx[pf_off + 11 * ccomps * dcomps];

            auto g_y_xxz = contr_buffer_pfxx[pf_off + 12 * ccomps * dcomps];

            auto g_y_xyy = contr_buffer_pfxx[pf_off + 13 * ccomps * dcomps];

            auto g_y_xyz = contr_buffer_pfxx[pf_off + 14 * ccomps * dcomps];

            auto g_y_xzz = contr_buffer_pfxx[pf_off + 15 * ccomps * dcomps];

            auto g_y_yyy = contr_buffer_pfxx[pf_off + 16 * ccomps * dcomps];

            auto g_y_yyz = contr_buffer_pfxx[pf_off + 17 * ccomps * dcomps];

            auto g_y_yzz = contr_buffer_pfxx[pf_off + 18 * ccomps * dcomps];

            auto g_y_zzz = contr_buffer_pfxx[pf_off + 19 * ccomps * dcomps];

            #pragma omp simd aligned(g_0_xxx, g_0_xxxy, g_0_xxy, g_0_xxyy, g_0_xxyz, g_0_xxz, g_0_xyy, g_0_xyyy, g_0_xyyz, g_0_xyz, g_0_xyzz, g_0_xzz, g_0_yyy, g_0_yyyy, g_0_yyyz, g_0_yyz, g_0_yyzz, g_0_yzz, g_0_yzzz, g_0_zzz, g_y_xxx, g_y_xxy, g_y_xxz, g_y_xyy, g_y_xyz, g_y_xzz, g_y_yyy, g_y_yyz, g_y_yzz, g_y_zzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_y_xxx[k] = -g_0_xxx[k] * ab_y + g_0_xxxy[k];

                g_y_xxy[k] = -g_0_xxy[k] * ab_y + g_0_xxyy[k];

                g_y_xxz[k] = -g_0_xxz[k] * ab_y + g_0_xxyz[k];

                g_y_xyy[k] = -g_0_xyy[k] * ab_y + g_0_xyyy[k];

                g_y_xyz[k] = -g_0_xyz[k] * ab_y + g_0_xyyz[k];

                g_y_xzz[k] = -g_0_xzz[k] * ab_y + g_0_xyzz[k];

                g_y_yyy[k] = -g_0_yyy[k] * ab_y + g_0_yyyy[k];

                g_y_yyz[k] = -g_0_yyz[k] * ab_y + g_0_yyyz[k];

                g_y_yzz[k] = -g_0_yzz[k] * ab_y + g_0_yyzz[k];

                g_y_zzz[k] = -g_0_zzz[k] * ab_y + g_0_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : contr_buffer_pfxx

            auto g_z_xxx = contr_buffer_pfxx[pf_off + 20 * ccomps * dcomps];

            auto g_z_xxy = contr_buffer_pfxx[pf_off + 21 * ccomps * dcomps];

            auto g_z_xxz = contr_buffer_pfxx[pf_off + 22 * ccomps * dcomps];

            auto g_z_xyy = contr_buffer_pfxx[pf_off + 23 * ccomps * dcomps];

            auto g_z_xyz = contr_buffer_pfxx[pf_off + 24 * ccomps * dcomps];

            auto g_z_xzz = contr_buffer_pfxx[pf_off + 25 * ccomps * dcomps];

            auto g_z_yyy = contr_buffer_pfxx[pf_off + 26 * ccomps * dcomps];

            auto g_z_yyz = contr_buffer_pfxx[pf_off + 27 * ccomps * dcomps];

            auto g_z_yzz = contr_buffer_pfxx[pf_off + 28 * ccomps * dcomps];

            auto g_z_zzz = contr_buffer_pfxx[pf_off + 29 * ccomps * dcomps];

            #pragma omp simd aligned(g_0_xxx, g_0_xxxz, g_0_xxy, g_0_xxyz, g_0_xxz, g_0_xxzz, g_0_xyy, g_0_xyyz, g_0_xyz, g_0_xyzz, g_0_xzz, g_0_xzzz, g_0_yyy, g_0_yyyz, g_0_yyz, g_0_yyzz, g_0_yzz, g_0_yzzz, g_0_zzz, g_0_zzzz, g_z_xxx, g_z_xxy, g_z_xxz, g_z_xyy, g_z_xyz, g_z_xzz, g_z_yyy, g_z_yyz, g_z_yzz, g_z_zzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_z_xxx[k] = -g_0_xxx[k] * ab_z + g_0_xxxz[k];

                g_z_xxy[k] = -g_0_xxy[k] * ab_z + g_0_xxyz[k];

                g_z_xxz[k] = -g_0_xxz[k] * ab_z + g_0_xxzz[k];

                g_z_xyy[k] = -g_0_xyy[k] * ab_z + g_0_xyyz[k];

                g_z_xyz[k] = -g_0_xyz[k] * ab_z + g_0_xyzz[k];

                g_z_xzz[k] = -g_0_xzz[k] * ab_z + g_0_xzzz[k];

                g_z_yyy[k] = -g_0_yyy[k] * ab_z + g_0_yyyz[k];

                g_z_yyz[k] = -g_0_yyz[k] * ab_z + g_0_yyzz[k];

                g_z_yzz[k] = -g_0_yzz[k] * ab_z + g_0_yzzz[k];

                g_z_zzz[k] = -g_0_zzz[k] * ab_z + g_0_zzzz[k];
            }
        }
    }
}

} // erirec namespace

