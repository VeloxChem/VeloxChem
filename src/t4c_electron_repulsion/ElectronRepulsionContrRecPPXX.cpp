#include "ElectronRepulsionContrRecPPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_hrr_electron_repulsion_ppxx(CSimdArray<double>& contr_buffer_ppxx,
                                     const CSimdArray<double>& contr_buffer_spxx,
                                     const CSimdArray<double>& contr_buffer_sdxx,
                                     const double ab_x,
                                     const double ab_y,
                                     const double ab_z,
                                     const int c_angmom,
                                     const int d_angmom) -> void
{
    const auto ndims = contr_buffer_ppxx.number_of_columns();

    const auto ccomps = tensor::number_of_spherical_components(c_angmom);

    const auto dcomps = tensor::number_of_spherical_components(d_angmom);

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : contr_buffer_spxx

            const auto sp_off = i * dcomps + j;

            auto g_0_x = contr_buffer_spxx[sp_off + 0 * ccomps * dcomps];

            auto g_0_y = contr_buffer_spxx[sp_off + 1 * ccomps * dcomps];

            auto g_0_z = contr_buffer_spxx[sp_off + 2 * ccomps * dcomps];

            /// Set up components of auxilary buffer : contr_buffer_sdxx

            const auto sd_off = i * dcomps + j;

            auto g_0_xx = contr_buffer_sdxx[sd_off + 0 * ccomps * dcomps];

            auto g_0_xy = contr_buffer_sdxx[sd_off + 1 * ccomps * dcomps];

            auto g_0_xz = contr_buffer_sdxx[sd_off + 2 * ccomps * dcomps];

            auto g_0_yy = contr_buffer_sdxx[sd_off + 3 * ccomps * dcomps];

            auto g_0_yz = contr_buffer_sdxx[sd_off + 4 * ccomps * dcomps];

            auto g_0_zz = contr_buffer_sdxx[sd_off + 5 * ccomps * dcomps];

            /// set up bra offset for contr_buffer_ppxx

            const auto pp_off = i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : contr_buffer_ppxx

            auto g_x_x = contr_buffer_ppxx[pp_off + 0 * ccomps * dcomps];

            auto g_x_y = contr_buffer_ppxx[pp_off + 1 * ccomps * dcomps];

            auto g_x_z = contr_buffer_ppxx[pp_off + 2 * ccomps * dcomps];

            #pragma omp simd aligned(g_0_x, g_0_xx, g_0_xy, g_0_xz, g_0_y, g_0_z, g_x_x, g_x_y, g_x_z  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_x_x[k] = -g_0_x[k] * ab_x + g_0_xx[k];

                g_x_y[k] = -g_0_y[k] * ab_x + g_0_xy[k];

                g_x_z[k] = -g_0_z[k] * ab_x + g_0_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : contr_buffer_ppxx

            auto g_y_x = contr_buffer_ppxx[pp_off + 3 * ccomps * dcomps];

            auto g_y_y = contr_buffer_ppxx[pp_off + 4 * ccomps * dcomps];

            auto g_y_z = contr_buffer_ppxx[pp_off + 5 * ccomps * dcomps];

            #pragma omp simd aligned(g_0_x, g_0_xy, g_0_y, g_0_yy, g_0_yz, g_0_z, g_y_x, g_y_y, g_y_z  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_y_x[k] = -g_0_x[k] * ab_y + g_0_xy[k];

                g_y_y[k] = -g_0_y[k] * ab_y + g_0_yy[k];

                g_y_z[k] = -g_0_z[k] * ab_y + g_0_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : contr_buffer_ppxx

            auto g_z_x = contr_buffer_ppxx[pp_off + 6 * ccomps * dcomps];

            auto g_z_y = contr_buffer_ppxx[pp_off + 7 * ccomps * dcomps];

            auto g_z_z = contr_buffer_ppxx[pp_off + 8 * ccomps * dcomps];

            #pragma omp simd aligned(g_0_x, g_0_xz, g_0_y, g_0_yz, g_0_z, g_0_zz, g_z_x, g_z_y, g_z_z  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_z_x[k] = -g_0_x[k] * ab_z + g_0_xz[k];

                g_z_y[k] = -g_0_y[k] * ab_z + g_0_yz[k];

                g_z_z[k] = -g_0_z[k] * ab_z + g_0_zz[k];
            }
        }
    }
}

} // erirec namespace

