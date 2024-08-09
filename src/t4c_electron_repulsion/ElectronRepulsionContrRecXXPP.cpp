#include "ElectronRepulsionContrRecXXPP.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_hrr_electron_repulsion_xxpp(CSimdArray<double>& contr_buffer_xxpp,
                                     const CSimdArray<double>& contr_buffer_xxsp,
                                     const CSimdArray<double>& contr_buffer_xxsd,
                                     const double* cd_x,
                                     const double* cd_y,
                                     const double* cd_z,
                                     const int a_angmom,
                                     const int b_angmom) -> void
{
    const auto ndims = contr_buffer_xxpp.number_of_columns();

    const auto acomps = tensor::number_of_cartesian_components(a_angmom);

    const auto bcomps = tensor::number_of_cartesian_components(b_angmom);

    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
            /// Set up components of auxilary buffer : contr_buffer_xxsp

            const auto sp_off = (i * bcomps + j) * 3;

            auto g_0_x = contr_buffer_xxsp[sp_off + 0];

            auto g_0_y = contr_buffer_xxsp[sp_off + 1];

            auto g_0_z = contr_buffer_xxsp[sp_off + 2];

            /// Set up components of auxilary buffer : contr_buffer_xxsd

            const auto sd_off = (i * bcomps + j) * 6;

            auto g_0_xx = contr_buffer_xxsd[sd_off + 0];

            auto g_0_xy = contr_buffer_xxsd[sd_off + 1];

            auto g_0_xz = contr_buffer_xxsd[sd_off + 2];

            auto g_0_yy = contr_buffer_xxsd[sd_off + 3];

            auto g_0_yz = contr_buffer_xxsd[sd_off + 4];

            auto g_0_zz = contr_buffer_xxsd[sd_off + 5];

            /// set up bra offset for contr_buffer_xxpp

            const auto pp_off = (i * bcomps + j) * 9;

            /// Set up 0-3 components of targeted buffer : contr_buffer_xxpp

            auto g_x_x = contr_buffer_xxpp[pp_off + 0];

            auto g_x_y = contr_buffer_xxpp[pp_off + 1];

            auto g_x_z = contr_buffer_xxpp[pp_off + 2];

            #pragma omp simd aligned(cd_x, g_0_x, g_0_xx, g_0_xy, g_0_xz, g_0_y, g_0_z, g_x_x, g_x_y, g_x_z  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_x_x[k] = -g_0_x[k] * cd_x[k] + g_0_xx[k];

                g_x_y[k] = -g_0_y[k] * cd_x[k] + g_0_xy[k];

                g_x_z[k] = -g_0_z[k] * cd_x[k] + g_0_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : contr_buffer_xxpp

            auto g_y_x = contr_buffer_xxpp[pp_off + 3];

            auto g_y_y = contr_buffer_xxpp[pp_off + 4];

            auto g_y_z = contr_buffer_xxpp[pp_off + 5];

            #pragma omp simd aligned(cd_y, g_0_x, g_0_xy, g_0_y, g_0_yy, g_0_yz, g_0_z, g_y_x, g_y_y, g_y_z  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_y_x[k] = -g_0_x[k] * cd_y[k] + g_0_xy[k];

                g_y_y[k] = -g_0_y[k] * cd_y[k] + g_0_yy[k];

                g_y_z[k] = -g_0_z[k] * cd_y[k] + g_0_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : contr_buffer_xxpp

            auto g_z_x = contr_buffer_xxpp[pp_off + 6];

            auto g_z_y = contr_buffer_xxpp[pp_off + 7];

            auto g_z_z = contr_buffer_xxpp[pp_off + 8];

            #pragma omp simd aligned(cd_z, g_0_x, g_0_xz, g_0_y, g_0_yz, g_0_z, g_0_zz, g_z_x, g_z_y, g_z_z  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_z_x[k] = -g_0_x[k] * cd_z[k] + g_0_xz[k];

                g_z_y[k] = -g_0_y[k] * cd_z[k] + g_0_yz[k];

                g_z_z[k] = -g_0_z[k] * cd_z[k] + g_0_zz[k];
            }
        }
    }
}

} // erirec namespace

