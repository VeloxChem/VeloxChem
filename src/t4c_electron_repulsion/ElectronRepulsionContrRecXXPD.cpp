#include "ElectronRepulsionContrRecXXPD.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_hrr_electron_repulsion_xxpd(CSimdArray<double>& contr_buffer_xxpd,
                                     const CSimdArray<double>& contr_buffer_xxsd,
                                     const CSimdArray<double>& contr_buffer_xxsf,
                                     const double* cd_x,
                                     const double* cd_y,
                                     const double* cd_z,
                                     const int a_angmom,
                                     const int b_angmom) -> void
{
    const auto ndims = contr_buffer_xxpd.number_of_columns();

    const auto acomps = tensor::number_of_cartesian_components(a_angmom);

    const auto bcomps = tensor::number_of_cartesian_components(b_angmom);

    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
            /// Set up components of auxilary buffer : contr_buffer_xxsd

            const auto sd_off = (i * bcomps + j) * 6;

            auto g_0_xx = contr_buffer_xxsd[sd_off + 0];

            auto g_0_xy = contr_buffer_xxsd[sd_off + 1];

            auto g_0_xz = contr_buffer_xxsd[sd_off + 2];

            auto g_0_yy = contr_buffer_xxsd[sd_off + 3];

            auto g_0_yz = contr_buffer_xxsd[sd_off + 4];

            auto g_0_zz = contr_buffer_xxsd[sd_off + 5];

            /// Set up components of auxilary buffer : contr_buffer_xxsf

            const auto sf_off = (i * bcomps + j) * 10;

            auto g_0_xxx = contr_buffer_xxsf[sf_off + 0];

            auto g_0_xxy = contr_buffer_xxsf[sf_off + 1];

            auto g_0_xxz = contr_buffer_xxsf[sf_off + 2];

            auto g_0_xyy = contr_buffer_xxsf[sf_off + 3];

            auto g_0_xyz = contr_buffer_xxsf[sf_off + 4];

            auto g_0_xzz = contr_buffer_xxsf[sf_off + 5];

            auto g_0_yyy = contr_buffer_xxsf[sf_off + 6];

            auto g_0_yyz = contr_buffer_xxsf[sf_off + 7];

            auto g_0_yzz = contr_buffer_xxsf[sf_off + 8];

            auto g_0_zzz = contr_buffer_xxsf[sf_off + 9];

            /// set up bra offset for contr_buffer_xxpd

            const auto pd_off = (i * bcomps + j) * 18;

            /// Set up 0-6 components of targeted buffer : contr_buffer_xxpd

            auto g_x_xx = contr_buffer_xxpd[pd_off + 0];

            auto g_x_xy = contr_buffer_xxpd[pd_off + 1];

            auto g_x_xz = contr_buffer_xxpd[pd_off + 2];

            auto g_x_yy = contr_buffer_xxpd[pd_off + 3];

            auto g_x_yz = contr_buffer_xxpd[pd_off + 4];

            auto g_x_zz = contr_buffer_xxpd[pd_off + 5];

            #pragma omp simd aligned(cd_x, g_0_xx, g_0_xxx, g_0_xxy, g_0_xxz, g_0_xy, g_0_xyy, g_0_xyz, g_0_xz, g_0_xzz, g_0_yy, g_0_yz, g_0_zz, g_x_xx, g_x_xy, g_x_xz, g_x_yy, g_x_yz, g_x_zz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_x_xx[k] = -g_0_xx[k] * cd_x[k] + g_0_xxx[k];

                g_x_xy[k] = -g_0_xy[k] * cd_x[k] + g_0_xxy[k];

                g_x_xz[k] = -g_0_xz[k] * cd_x[k] + g_0_xxz[k];

                g_x_yy[k] = -g_0_yy[k] * cd_x[k] + g_0_xyy[k];

                g_x_yz[k] = -g_0_yz[k] * cd_x[k] + g_0_xyz[k];

                g_x_zz[k] = -g_0_zz[k] * cd_x[k] + g_0_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : contr_buffer_xxpd

            auto g_y_xx = contr_buffer_xxpd[pd_off + 6];

            auto g_y_xy = contr_buffer_xxpd[pd_off + 7];

            auto g_y_xz = contr_buffer_xxpd[pd_off + 8];

            auto g_y_yy = contr_buffer_xxpd[pd_off + 9];

            auto g_y_yz = contr_buffer_xxpd[pd_off + 10];

            auto g_y_zz = contr_buffer_xxpd[pd_off + 11];

            #pragma omp simd aligned(cd_y, g_0_xx, g_0_xxy, g_0_xy, g_0_xyy, g_0_xyz, g_0_xz, g_0_yy, g_0_yyy, g_0_yyz, g_0_yz, g_0_yzz, g_0_zz, g_y_xx, g_y_xy, g_y_xz, g_y_yy, g_y_yz, g_y_zz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_y_xx[k] = -g_0_xx[k] * cd_y[k] + g_0_xxy[k];

                g_y_xy[k] = -g_0_xy[k] * cd_y[k] + g_0_xyy[k];

                g_y_xz[k] = -g_0_xz[k] * cd_y[k] + g_0_xyz[k];

                g_y_yy[k] = -g_0_yy[k] * cd_y[k] + g_0_yyy[k];

                g_y_yz[k] = -g_0_yz[k] * cd_y[k] + g_0_yyz[k];

                g_y_zz[k] = -g_0_zz[k] * cd_y[k] + g_0_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : contr_buffer_xxpd

            auto g_z_xx = contr_buffer_xxpd[pd_off + 12];

            auto g_z_xy = contr_buffer_xxpd[pd_off + 13];

            auto g_z_xz = contr_buffer_xxpd[pd_off + 14];

            auto g_z_yy = contr_buffer_xxpd[pd_off + 15];

            auto g_z_yz = contr_buffer_xxpd[pd_off + 16];

            auto g_z_zz = contr_buffer_xxpd[pd_off + 17];

            #pragma omp simd aligned(cd_z, g_0_xx, g_0_xxz, g_0_xy, g_0_xyz, g_0_xz, g_0_xzz, g_0_yy, g_0_yyz, g_0_yz, g_0_yzz, g_0_zz, g_0_zzz, g_z_xx, g_z_xy, g_z_xz, g_z_yy, g_z_yz, g_z_zz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_z_xx[k] = -g_0_xx[k] * cd_z[k] + g_0_xxz[k];

                g_z_xy[k] = -g_0_xy[k] * cd_z[k] + g_0_xyz[k];

                g_z_xz[k] = -g_0_xz[k] * cd_z[k] + g_0_xzz[k];

                g_z_yy[k] = -g_0_yy[k] * cd_z[k] + g_0_yyz[k];

                g_z_yz[k] = -g_0_yz[k] * cd_z[k] + g_0_yzz[k];

                g_z_zz[k] = -g_0_zz[k] * cd_z[k] + g_0_zzz[k];
            }
        }
    }
}

} // erirec namespace

