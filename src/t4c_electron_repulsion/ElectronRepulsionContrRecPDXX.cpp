#include "ElectronRepulsionContrRecPDXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_hrr_electron_repulsion_pdxx(CSimdArray<double>& contr_buffer_pdxx,
                                     const CSimdArray<double>& contr_buffer_sdxx,
                                     const CSimdArray<double>& contr_buffer_sfxx,
                                     const double ab_x,
                                     const double ab_y,
                                     const double ab_z,
                                     const int c_angmom,
                                     const int d_angmom) -> void
{
    const auto ndims = contr_buffer_pdxx.number_of_columns();

    const auto ccomps = tensor::number_of_spherical_components(c_angmom);

    const auto dcomps = tensor::number_of_spherical_components(d_angmom);

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : contr_buffer_sdxx

            const auto sd_off = i * dcomps + j;

            auto g_0_xx = contr_buffer_sdxx[sd_off + 0 * ccomps * dcomps];

            auto g_0_xy = contr_buffer_sdxx[sd_off + 1 * ccomps * dcomps];

            auto g_0_xz = contr_buffer_sdxx[sd_off + 2 * ccomps * dcomps];

            auto g_0_yy = contr_buffer_sdxx[sd_off + 3 * ccomps * dcomps];

            auto g_0_yz = contr_buffer_sdxx[sd_off + 4 * ccomps * dcomps];

            auto g_0_zz = contr_buffer_sdxx[sd_off + 5 * ccomps * dcomps];

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

            /// set up bra offset for contr_buffer_pdxx

            const auto pd_off = i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : contr_buffer_pdxx

            auto g_x_xx = contr_buffer_pdxx[pd_off + 0 * ccomps * dcomps];

            auto g_x_xy = contr_buffer_pdxx[pd_off + 1 * ccomps * dcomps];

            auto g_x_xz = contr_buffer_pdxx[pd_off + 2 * ccomps * dcomps];

            auto g_x_yy = contr_buffer_pdxx[pd_off + 3 * ccomps * dcomps];

            auto g_x_yz = contr_buffer_pdxx[pd_off + 4 * ccomps * dcomps];

            auto g_x_zz = contr_buffer_pdxx[pd_off + 5 * ccomps * dcomps];

            #pragma omp simd aligned(g_0_xx, g_0_xxx, g_0_xxy, g_0_xxz, g_0_xy, g_0_xyy, g_0_xyz, g_0_xz, g_0_xzz, g_0_yy, g_0_yz, g_0_zz, g_x_xx, g_x_xy, g_x_xz, g_x_yy, g_x_yz, g_x_zz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_x_xx[k] = -g_0_xx[k] * ab_x + g_0_xxx[k];

                g_x_xy[k] = -g_0_xy[k] * ab_x + g_0_xxy[k];

                g_x_xz[k] = -g_0_xz[k] * ab_x + g_0_xxz[k];

                g_x_yy[k] = -g_0_yy[k] * ab_x + g_0_xyy[k];

                g_x_yz[k] = -g_0_yz[k] * ab_x + g_0_xyz[k];

                g_x_zz[k] = -g_0_zz[k] * ab_x + g_0_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : contr_buffer_pdxx

            auto g_y_xx = contr_buffer_pdxx[pd_off + 6 * ccomps * dcomps];

            auto g_y_xy = contr_buffer_pdxx[pd_off + 7 * ccomps * dcomps];

            auto g_y_xz = contr_buffer_pdxx[pd_off + 8 * ccomps * dcomps];

            auto g_y_yy = contr_buffer_pdxx[pd_off + 9 * ccomps * dcomps];

            auto g_y_yz = contr_buffer_pdxx[pd_off + 10 * ccomps * dcomps];

            auto g_y_zz = contr_buffer_pdxx[pd_off + 11 * ccomps * dcomps];

            #pragma omp simd aligned(g_0_xx, g_0_xxy, g_0_xy, g_0_xyy, g_0_xyz, g_0_xz, g_0_yy, g_0_yyy, g_0_yyz, g_0_yz, g_0_yzz, g_0_zz, g_y_xx, g_y_xy, g_y_xz, g_y_yy, g_y_yz, g_y_zz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_y_xx[k] = -g_0_xx[k] * ab_y + g_0_xxy[k];

                g_y_xy[k] = -g_0_xy[k] * ab_y + g_0_xyy[k];

                g_y_xz[k] = -g_0_xz[k] * ab_y + g_0_xyz[k];

                g_y_yy[k] = -g_0_yy[k] * ab_y + g_0_yyy[k];

                g_y_yz[k] = -g_0_yz[k] * ab_y + g_0_yyz[k];

                g_y_zz[k] = -g_0_zz[k] * ab_y + g_0_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : contr_buffer_pdxx

            auto g_z_xx = contr_buffer_pdxx[pd_off + 12 * ccomps * dcomps];

            auto g_z_xy = contr_buffer_pdxx[pd_off + 13 * ccomps * dcomps];

            auto g_z_xz = contr_buffer_pdxx[pd_off + 14 * ccomps * dcomps];

            auto g_z_yy = contr_buffer_pdxx[pd_off + 15 * ccomps * dcomps];

            auto g_z_yz = contr_buffer_pdxx[pd_off + 16 * ccomps * dcomps];

            auto g_z_zz = contr_buffer_pdxx[pd_off + 17 * ccomps * dcomps];

            #pragma omp simd aligned(g_0_xx, g_0_xxz, g_0_xy, g_0_xyz, g_0_xz, g_0_xzz, g_0_yy, g_0_yyz, g_0_yz, g_0_yzz, g_0_zz, g_0_zzz, g_z_xx, g_z_xy, g_z_xz, g_z_yy, g_z_yz, g_z_zz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_z_xx[k] = -g_0_xx[k] * ab_z + g_0_xxz[k];

                g_z_xy[k] = -g_0_xy[k] * ab_z + g_0_xyz[k];

                g_z_xz[k] = -g_0_xz[k] * ab_z + g_0_xzz[k];

                g_z_yy[k] = -g_0_yy[k] * ab_z + g_0_yyz[k];

                g_z_yz[k] = -g_0_yz[k] * ab_z + g_0_yzz[k];

                g_z_zz[k] = -g_0_zz[k] * ab_z + g_0_zzz[k];
            }
        }
    }
}

} // erirec namespace

