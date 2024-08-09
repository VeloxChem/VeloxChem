#include "ElectronRepulsionContrRecXXDD.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_hrr_electron_repulsion_xxdd(CSimdArray<double>& contr_buffer_xxdd,
                                     const CSimdArray<double>& contr_buffer_xxpd,
                                     const CSimdArray<double>& contr_buffer_xxpf,
                                     const double* cd_x,
                                     const double* cd_y,
                                     const double* cd_z,
                                     const int a_angmom,
                                     const int b_angmom) -> void
{
    const auto ndims = contr_buffer_xxdd.number_of_columns();

    const auto acomps = tensor::number_of_cartesian_components(a_angmom);

    const auto bcomps = tensor::number_of_cartesian_components(b_angmom);

    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
            /// Set up components of auxilary buffer : contr_buffer_xxpd

            const auto pd_off = (i * bcomps + j) * 18;

            auto g_x_xx = contr_buffer_xxpd[pd_off + 0];

            auto g_x_xy = contr_buffer_xxpd[pd_off + 1];

            auto g_x_xz = contr_buffer_xxpd[pd_off + 2];

            auto g_x_yy = contr_buffer_xxpd[pd_off + 3];

            auto g_x_yz = contr_buffer_xxpd[pd_off + 4];

            auto g_x_zz = contr_buffer_xxpd[pd_off + 5];

            auto g_y_xx = contr_buffer_xxpd[pd_off + 6];

            auto g_y_xy = contr_buffer_xxpd[pd_off + 7];

            auto g_y_xz = contr_buffer_xxpd[pd_off + 8];

            auto g_y_yy = contr_buffer_xxpd[pd_off + 9];

            auto g_y_yz = contr_buffer_xxpd[pd_off + 10];

            auto g_y_zz = contr_buffer_xxpd[pd_off + 11];

            auto g_z_xx = contr_buffer_xxpd[pd_off + 12];

            auto g_z_xy = contr_buffer_xxpd[pd_off + 13];

            auto g_z_xz = contr_buffer_xxpd[pd_off + 14];

            auto g_z_yy = contr_buffer_xxpd[pd_off + 15];

            auto g_z_yz = contr_buffer_xxpd[pd_off + 16];

            auto g_z_zz = contr_buffer_xxpd[pd_off + 17];

            /// Set up components of auxilary buffer : contr_buffer_xxpf

            const auto pf_off = (i * bcomps + j) * 30;

            auto g_x_xxx = contr_buffer_xxpf[pf_off + 0];

            auto g_x_xxy = contr_buffer_xxpf[pf_off + 1];

            auto g_x_xxz = contr_buffer_xxpf[pf_off + 2];

            auto g_x_xyy = contr_buffer_xxpf[pf_off + 3];

            auto g_x_xyz = contr_buffer_xxpf[pf_off + 4];

            auto g_x_xzz = contr_buffer_xxpf[pf_off + 5];

            auto g_y_xxx = contr_buffer_xxpf[pf_off + 10];

            auto g_y_xxy = contr_buffer_xxpf[pf_off + 11];

            auto g_y_xxz = contr_buffer_xxpf[pf_off + 12];

            auto g_y_xyy = contr_buffer_xxpf[pf_off + 13];

            auto g_y_xyz = contr_buffer_xxpf[pf_off + 14];

            auto g_y_xzz = contr_buffer_xxpf[pf_off + 15];

            auto g_y_yyy = contr_buffer_xxpf[pf_off + 16];

            auto g_y_yyz = contr_buffer_xxpf[pf_off + 17];

            auto g_y_yzz = contr_buffer_xxpf[pf_off + 18];

            auto g_z_xxx = contr_buffer_xxpf[pf_off + 20];

            auto g_z_xxy = contr_buffer_xxpf[pf_off + 21];

            auto g_z_xxz = contr_buffer_xxpf[pf_off + 22];

            auto g_z_xyy = contr_buffer_xxpf[pf_off + 23];

            auto g_z_xyz = contr_buffer_xxpf[pf_off + 24];

            auto g_z_xzz = contr_buffer_xxpf[pf_off + 25];

            auto g_z_yyy = contr_buffer_xxpf[pf_off + 26];

            auto g_z_yyz = contr_buffer_xxpf[pf_off + 27];

            auto g_z_yzz = contr_buffer_xxpf[pf_off + 28];

            auto g_z_zzz = contr_buffer_xxpf[pf_off + 29];

            /// set up bra offset for contr_buffer_xxdd

            const auto dd_off = (i * bcomps + j) * 36;

            /// Set up 0-6 components of targeted buffer : contr_buffer_xxdd

            auto g_xx_xx = contr_buffer_xxdd[dd_off + 0];

            auto g_xx_xy = contr_buffer_xxdd[dd_off + 1];

            auto g_xx_xz = contr_buffer_xxdd[dd_off + 2];

            auto g_xx_yy = contr_buffer_xxdd[dd_off + 3];

            auto g_xx_yz = contr_buffer_xxdd[dd_off + 4];

            auto g_xx_zz = contr_buffer_xxdd[dd_off + 5];

            #pragma omp simd aligned(cd_x, g_x_xx, g_x_xxx, g_x_xxy, g_x_xxz, g_x_xy, g_x_xyy, g_x_xyz, g_x_xz, g_x_xzz, g_x_yy, g_x_yz, g_x_zz, g_xx_xx, g_xx_xy, g_xx_xz, g_xx_yy, g_xx_yz, g_xx_zz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_xx_xx[k] = -g_x_xx[k] * cd_x[k] + g_x_xxx[k];

                g_xx_xy[k] = -g_x_xy[k] * cd_x[k] + g_x_xxy[k];

                g_xx_xz[k] = -g_x_xz[k] * cd_x[k] + g_x_xxz[k];

                g_xx_yy[k] = -g_x_yy[k] * cd_x[k] + g_x_xyy[k];

                g_xx_yz[k] = -g_x_yz[k] * cd_x[k] + g_x_xyz[k];

                g_xx_zz[k] = -g_x_zz[k] * cd_x[k] + g_x_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : contr_buffer_xxdd

            auto g_xy_xx = contr_buffer_xxdd[dd_off + 6];

            auto g_xy_xy = contr_buffer_xxdd[dd_off + 7];

            auto g_xy_xz = contr_buffer_xxdd[dd_off + 8];

            auto g_xy_yy = contr_buffer_xxdd[dd_off + 9];

            auto g_xy_yz = contr_buffer_xxdd[dd_off + 10];

            auto g_xy_zz = contr_buffer_xxdd[dd_off + 11];

            #pragma omp simd aligned(cd_x, g_xy_xx, g_xy_xy, g_xy_xz, g_xy_yy, g_xy_yz, g_xy_zz, g_y_xx, g_y_xxx, g_y_xxy, g_y_xxz, g_y_xy, g_y_xyy, g_y_xyz, g_y_xz, g_y_xzz, g_y_yy, g_y_yz, g_y_zz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_xy_xx[k] = -g_y_xx[k] * cd_x[k] + g_y_xxx[k];

                g_xy_xy[k] = -g_y_xy[k] * cd_x[k] + g_y_xxy[k];

                g_xy_xz[k] = -g_y_xz[k] * cd_x[k] + g_y_xxz[k];

                g_xy_yy[k] = -g_y_yy[k] * cd_x[k] + g_y_xyy[k];

                g_xy_yz[k] = -g_y_yz[k] * cd_x[k] + g_y_xyz[k];

                g_xy_zz[k] = -g_y_zz[k] * cd_x[k] + g_y_xzz[k];
            }

            /// Set up 12-18 components of targeted buffer : contr_buffer_xxdd

            auto g_xz_xx = contr_buffer_xxdd[dd_off + 12];

            auto g_xz_xy = contr_buffer_xxdd[dd_off + 13];

            auto g_xz_xz = contr_buffer_xxdd[dd_off + 14];

            auto g_xz_yy = contr_buffer_xxdd[dd_off + 15];

            auto g_xz_yz = contr_buffer_xxdd[dd_off + 16];

            auto g_xz_zz = contr_buffer_xxdd[dd_off + 17];

            #pragma omp simd aligned(cd_x, g_xz_xx, g_xz_xy, g_xz_xz, g_xz_yy, g_xz_yz, g_xz_zz, g_z_xx, g_z_xxx, g_z_xxy, g_z_xxz, g_z_xy, g_z_xyy, g_z_xyz, g_z_xz, g_z_xzz, g_z_yy, g_z_yz, g_z_zz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_xz_xx[k] = -g_z_xx[k] * cd_x[k] + g_z_xxx[k];

                g_xz_xy[k] = -g_z_xy[k] * cd_x[k] + g_z_xxy[k];

                g_xz_xz[k] = -g_z_xz[k] * cd_x[k] + g_z_xxz[k];

                g_xz_yy[k] = -g_z_yy[k] * cd_x[k] + g_z_xyy[k];

                g_xz_yz[k] = -g_z_yz[k] * cd_x[k] + g_z_xyz[k];

                g_xz_zz[k] = -g_z_zz[k] * cd_x[k] + g_z_xzz[k];
            }

            /// Set up 18-24 components of targeted buffer : contr_buffer_xxdd

            auto g_yy_xx = contr_buffer_xxdd[dd_off + 18];

            auto g_yy_xy = contr_buffer_xxdd[dd_off + 19];

            auto g_yy_xz = contr_buffer_xxdd[dd_off + 20];

            auto g_yy_yy = contr_buffer_xxdd[dd_off + 21];

            auto g_yy_yz = contr_buffer_xxdd[dd_off + 22];

            auto g_yy_zz = contr_buffer_xxdd[dd_off + 23];

            #pragma omp simd aligned(cd_y, g_y_xx, g_y_xxy, g_y_xy, g_y_xyy, g_y_xyz, g_y_xz, g_y_yy, g_y_yyy, g_y_yyz, g_y_yz, g_y_yzz, g_y_zz, g_yy_xx, g_yy_xy, g_yy_xz, g_yy_yy, g_yy_yz, g_yy_zz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_yy_xx[k] = -g_y_xx[k] * cd_y[k] + g_y_xxy[k];

                g_yy_xy[k] = -g_y_xy[k] * cd_y[k] + g_y_xyy[k];

                g_yy_xz[k] = -g_y_xz[k] * cd_y[k] + g_y_xyz[k];

                g_yy_yy[k] = -g_y_yy[k] * cd_y[k] + g_y_yyy[k];

                g_yy_yz[k] = -g_y_yz[k] * cd_y[k] + g_y_yyz[k];

                g_yy_zz[k] = -g_y_zz[k] * cd_y[k] + g_y_yzz[k];
            }

            /// Set up 24-30 components of targeted buffer : contr_buffer_xxdd

            auto g_yz_xx = contr_buffer_xxdd[dd_off + 24];

            auto g_yz_xy = contr_buffer_xxdd[dd_off + 25];

            auto g_yz_xz = contr_buffer_xxdd[dd_off + 26];

            auto g_yz_yy = contr_buffer_xxdd[dd_off + 27];

            auto g_yz_yz = contr_buffer_xxdd[dd_off + 28];

            auto g_yz_zz = contr_buffer_xxdd[dd_off + 29];

            #pragma omp simd aligned(cd_y, g_yz_xx, g_yz_xy, g_yz_xz, g_yz_yy, g_yz_yz, g_yz_zz, g_z_xx, g_z_xxy, g_z_xy, g_z_xyy, g_z_xyz, g_z_xz, g_z_yy, g_z_yyy, g_z_yyz, g_z_yz, g_z_yzz, g_z_zz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_yz_xx[k] = -g_z_xx[k] * cd_y[k] + g_z_xxy[k];

                g_yz_xy[k] = -g_z_xy[k] * cd_y[k] + g_z_xyy[k];

                g_yz_xz[k] = -g_z_xz[k] * cd_y[k] + g_z_xyz[k];

                g_yz_yy[k] = -g_z_yy[k] * cd_y[k] + g_z_yyy[k];

                g_yz_yz[k] = -g_z_yz[k] * cd_y[k] + g_z_yyz[k];

                g_yz_zz[k] = -g_z_zz[k] * cd_y[k] + g_z_yzz[k];
            }

            /// Set up 30-36 components of targeted buffer : contr_buffer_xxdd

            auto g_zz_xx = contr_buffer_xxdd[dd_off + 30];

            auto g_zz_xy = contr_buffer_xxdd[dd_off + 31];

            auto g_zz_xz = contr_buffer_xxdd[dd_off + 32];

            auto g_zz_yy = contr_buffer_xxdd[dd_off + 33];

            auto g_zz_yz = contr_buffer_xxdd[dd_off + 34];

            auto g_zz_zz = contr_buffer_xxdd[dd_off + 35];

            #pragma omp simd aligned(cd_z, g_z_xx, g_z_xxz, g_z_xy, g_z_xyz, g_z_xz, g_z_xzz, g_z_yy, g_z_yyz, g_z_yz, g_z_yzz, g_z_zz, g_z_zzz, g_zz_xx, g_zz_xy, g_zz_xz, g_zz_yy, g_zz_yz, g_zz_zz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_zz_xx[k] = -g_z_xx[k] * cd_z[k] + g_z_xxz[k];

                g_zz_xy[k] = -g_z_xy[k] * cd_z[k] + g_z_xyz[k];

                g_zz_xz[k] = -g_z_xz[k] * cd_z[k] + g_z_xzz[k];

                g_zz_yy[k] = -g_z_yy[k] * cd_z[k] + g_z_yyz[k];

                g_zz_yz[k] = -g_z_yz[k] * cd_z[k] + g_z_yzz[k];

                g_zz_zz[k] = -g_z_zz[k] * cd_z[k] + g_z_zzz[k];
            }
        }
    }
}

} // erirec namespace

