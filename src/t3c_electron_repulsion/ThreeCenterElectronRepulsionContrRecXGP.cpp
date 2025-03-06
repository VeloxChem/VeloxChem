#include "ThreeCenterElectronRepulsionContrRecXGP.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_hrr_electron_repulsion_xgp(CSimdArray<double>& cbuffer,
                                const size_t idx_xgp,
                                const size_t idx_xfp,
                                const size_t idx_xfd,
                                const CSimdArray<double>& factors,
                                const size_t idx_cd,
                                const int a_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        /// Set up components of auxilary buffer : SFP

        const auto fp_off = idx_xfp + i * 30;

        auto g_xxx_x = cbuffer.data(fp_off + 0);

        auto g_xxx_y = cbuffer.data(fp_off + 1);

        auto g_xxx_z = cbuffer.data(fp_off + 2);

        auto g_xxy_x = cbuffer.data(fp_off + 3);

        auto g_xxy_y = cbuffer.data(fp_off + 4);

        auto g_xxy_z = cbuffer.data(fp_off + 5);

        auto g_xxz_x = cbuffer.data(fp_off + 6);

        auto g_xxz_y = cbuffer.data(fp_off + 7);

        auto g_xxz_z = cbuffer.data(fp_off + 8);

        auto g_xyy_x = cbuffer.data(fp_off + 9);

        auto g_xyy_y = cbuffer.data(fp_off + 10);

        auto g_xyy_z = cbuffer.data(fp_off + 11);

        auto g_xyz_x = cbuffer.data(fp_off + 12);

        auto g_xyz_y = cbuffer.data(fp_off + 13);

        auto g_xyz_z = cbuffer.data(fp_off + 14);

        auto g_xzz_x = cbuffer.data(fp_off + 15);

        auto g_xzz_y = cbuffer.data(fp_off + 16);

        auto g_xzz_z = cbuffer.data(fp_off + 17);

        auto g_yyy_x = cbuffer.data(fp_off + 18);

        auto g_yyy_y = cbuffer.data(fp_off + 19);

        auto g_yyy_z = cbuffer.data(fp_off + 20);

        auto g_yyz_x = cbuffer.data(fp_off + 21);

        auto g_yyz_y = cbuffer.data(fp_off + 22);

        auto g_yyz_z = cbuffer.data(fp_off + 23);

        auto g_yzz_x = cbuffer.data(fp_off + 24);

        auto g_yzz_y = cbuffer.data(fp_off + 25);

        auto g_yzz_z = cbuffer.data(fp_off + 26);

        auto g_zzz_x = cbuffer.data(fp_off + 27);

        auto g_zzz_y = cbuffer.data(fp_off + 28);

        auto g_zzz_z = cbuffer.data(fp_off + 29);

        /// Set up components of auxilary buffer : SFD

        const auto fd_off = idx_xfd + i * 60;

        auto g_xxx_xx = cbuffer.data(fd_off + 0);

        auto g_xxx_xy = cbuffer.data(fd_off + 1);

        auto g_xxx_xz = cbuffer.data(fd_off + 2);

        auto g_xxy_xx = cbuffer.data(fd_off + 6);

        auto g_xxy_xy = cbuffer.data(fd_off + 7);

        auto g_xxy_xz = cbuffer.data(fd_off + 8);

        auto g_xxz_xx = cbuffer.data(fd_off + 12);

        auto g_xxz_xy = cbuffer.data(fd_off + 13);

        auto g_xxz_xz = cbuffer.data(fd_off + 14);

        auto g_xyy_xx = cbuffer.data(fd_off + 18);

        auto g_xyy_xy = cbuffer.data(fd_off + 19);

        auto g_xyy_xz = cbuffer.data(fd_off + 20);

        auto g_xyz_xx = cbuffer.data(fd_off + 24);

        auto g_xyz_xy = cbuffer.data(fd_off + 25);

        auto g_xyz_xz = cbuffer.data(fd_off + 26);

        auto g_xzz_xx = cbuffer.data(fd_off + 30);

        auto g_xzz_xy = cbuffer.data(fd_off + 31);

        auto g_xzz_xz = cbuffer.data(fd_off + 32);

        auto g_yyy_xx = cbuffer.data(fd_off + 36);

        auto g_yyy_xy = cbuffer.data(fd_off + 37);

        auto g_yyy_xz = cbuffer.data(fd_off + 38);

        auto g_yyy_yy = cbuffer.data(fd_off + 39);

        auto g_yyy_yz = cbuffer.data(fd_off + 40);

        auto g_yyz_xx = cbuffer.data(fd_off + 42);

        auto g_yyz_xy = cbuffer.data(fd_off + 43);

        auto g_yyz_xz = cbuffer.data(fd_off + 44);

        auto g_yyz_yy = cbuffer.data(fd_off + 45);

        auto g_yyz_yz = cbuffer.data(fd_off + 46);

        auto g_yzz_xx = cbuffer.data(fd_off + 48);

        auto g_yzz_xy = cbuffer.data(fd_off + 49);

        auto g_yzz_xz = cbuffer.data(fd_off + 50);

        auto g_yzz_yy = cbuffer.data(fd_off + 51);

        auto g_yzz_yz = cbuffer.data(fd_off + 52);

        auto g_zzz_xx = cbuffer.data(fd_off + 54);

        auto g_zzz_xy = cbuffer.data(fd_off + 55);

        auto g_zzz_xz = cbuffer.data(fd_off + 56);

        auto g_zzz_yy = cbuffer.data(fd_off + 57);

        auto g_zzz_yz = cbuffer.data(fd_off + 58);

        auto g_zzz_zz = cbuffer.data(fd_off + 59);

        /// set up bra offset for contr_buffer_xgp

        const auto gp_off = idx_xgp + i * 45;

        /// Set up 0-3 components of targeted buffer : cbuffer.data(

        auto g_xxxx_x = cbuffer.data(gp_off + 0);

        auto g_xxxx_y = cbuffer.data(gp_off + 1);

        auto g_xxxx_z = cbuffer.data(gp_off + 2);

        #pragma omp simd aligned(cd_x, g_xxx_x, g_xxx_xx, g_xxx_xy, g_xxx_xz, g_xxx_y, g_xxx_z, g_xxxx_x, g_xxxx_y, g_xxxx_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxxx_x[k] = -g_xxx_x[k] * cd_x[k] + g_xxx_xx[k];

            g_xxxx_y[k] = -g_xxx_y[k] * cd_x[k] + g_xxx_xy[k];

            g_xxxx_z[k] = -g_xxx_z[k] * cd_x[k] + g_xxx_xz[k];
        }

        /// Set up 3-6 components of targeted buffer : cbuffer.data(

        auto g_xxxy_x = cbuffer.data(gp_off + 3);

        auto g_xxxy_y = cbuffer.data(gp_off + 4);

        auto g_xxxy_z = cbuffer.data(gp_off + 5);

        #pragma omp simd aligned(cd_x, g_xxxy_x, g_xxxy_y, g_xxxy_z, g_xxy_x, g_xxy_xx, g_xxy_xy, g_xxy_xz, g_xxy_y, g_xxy_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxxy_x[k] = -g_xxy_x[k] * cd_x[k] + g_xxy_xx[k];

            g_xxxy_y[k] = -g_xxy_y[k] * cd_x[k] + g_xxy_xy[k];

            g_xxxy_z[k] = -g_xxy_z[k] * cd_x[k] + g_xxy_xz[k];
        }

        /// Set up 6-9 components of targeted buffer : cbuffer.data(

        auto g_xxxz_x = cbuffer.data(gp_off + 6);

        auto g_xxxz_y = cbuffer.data(gp_off + 7);

        auto g_xxxz_z = cbuffer.data(gp_off + 8);

        #pragma omp simd aligned(cd_x, g_xxxz_x, g_xxxz_y, g_xxxz_z, g_xxz_x, g_xxz_xx, g_xxz_xy, g_xxz_xz, g_xxz_y, g_xxz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxxz_x[k] = -g_xxz_x[k] * cd_x[k] + g_xxz_xx[k];

            g_xxxz_y[k] = -g_xxz_y[k] * cd_x[k] + g_xxz_xy[k];

            g_xxxz_z[k] = -g_xxz_z[k] * cd_x[k] + g_xxz_xz[k];
        }

        /// Set up 9-12 components of targeted buffer : cbuffer.data(

        auto g_xxyy_x = cbuffer.data(gp_off + 9);

        auto g_xxyy_y = cbuffer.data(gp_off + 10);

        auto g_xxyy_z = cbuffer.data(gp_off + 11);

        #pragma omp simd aligned(cd_x, g_xxyy_x, g_xxyy_y, g_xxyy_z, g_xyy_x, g_xyy_xx, g_xyy_xy, g_xyy_xz, g_xyy_y, g_xyy_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxyy_x[k] = -g_xyy_x[k] * cd_x[k] + g_xyy_xx[k];

            g_xxyy_y[k] = -g_xyy_y[k] * cd_x[k] + g_xyy_xy[k];

            g_xxyy_z[k] = -g_xyy_z[k] * cd_x[k] + g_xyy_xz[k];
        }

        /// Set up 12-15 components of targeted buffer : cbuffer.data(

        auto g_xxyz_x = cbuffer.data(gp_off + 12);

        auto g_xxyz_y = cbuffer.data(gp_off + 13);

        auto g_xxyz_z = cbuffer.data(gp_off + 14);

        #pragma omp simd aligned(cd_x, g_xxyz_x, g_xxyz_y, g_xxyz_z, g_xyz_x, g_xyz_xx, g_xyz_xy, g_xyz_xz, g_xyz_y, g_xyz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxyz_x[k] = -g_xyz_x[k] * cd_x[k] + g_xyz_xx[k];

            g_xxyz_y[k] = -g_xyz_y[k] * cd_x[k] + g_xyz_xy[k];

            g_xxyz_z[k] = -g_xyz_z[k] * cd_x[k] + g_xyz_xz[k];
        }

        /// Set up 15-18 components of targeted buffer : cbuffer.data(

        auto g_xxzz_x = cbuffer.data(gp_off + 15);

        auto g_xxzz_y = cbuffer.data(gp_off + 16);

        auto g_xxzz_z = cbuffer.data(gp_off + 17);

        #pragma omp simd aligned(cd_x, g_xxzz_x, g_xxzz_y, g_xxzz_z, g_xzz_x, g_xzz_xx, g_xzz_xy, g_xzz_xz, g_xzz_y, g_xzz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxzz_x[k] = -g_xzz_x[k] * cd_x[k] + g_xzz_xx[k];

            g_xxzz_y[k] = -g_xzz_y[k] * cd_x[k] + g_xzz_xy[k];

            g_xxzz_z[k] = -g_xzz_z[k] * cd_x[k] + g_xzz_xz[k];
        }

        /// Set up 18-21 components of targeted buffer : cbuffer.data(

        auto g_xyyy_x = cbuffer.data(gp_off + 18);

        auto g_xyyy_y = cbuffer.data(gp_off + 19);

        auto g_xyyy_z = cbuffer.data(gp_off + 20);

        #pragma omp simd aligned(cd_x, g_xyyy_x, g_xyyy_y, g_xyyy_z, g_yyy_x, g_yyy_xx, g_yyy_xy, g_yyy_xz, g_yyy_y, g_yyy_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xyyy_x[k] = -g_yyy_x[k] * cd_x[k] + g_yyy_xx[k];

            g_xyyy_y[k] = -g_yyy_y[k] * cd_x[k] + g_yyy_xy[k];

            g_xyyy_z[k] = -g_yyy_z[k] * cd_x[k] + g_yyy_xz[k];
        }

        /// Set up 21-24 components of targeted buffer : cbuffer.data(

        auto g_xyyz_x = cbuffer.data(gp_off + 21);

        auto g_xyyz_y = cbuffer.data(gp_off + 22);

        auto g_xyyz_z = cbuffer.data(gp_off + 23);

        #pragma omp simd aligned(cd_x, g_xyyz_x, g_xyyz_y, g_xyyz_z, g_yyz_x, g_yyz_xx, g_yyz_xy, g_yyz_xz, g_yyz_y, g_yyz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xyyz_x[k] = -g_yyz_x[k] * cd_x[k] + g_yyz_xx[k];

            g_xyyz_y[k] = -g_yyz_y[k] * cd_x[k] + g_yyz_xy[k];

            g_xyyz_z[k] = -g_yyz_z[k] * cd_x[k] + g_yyz_xz[k];
        }

        /// Set up 24-27 components of targeted buffer : cbuffer.data(

        auto g_xyzz_x = cbuffer.data(gp_off + 24);

        auto g_xyzz_y = cbuffer.data(gp_off + 25);

        auto g_xyzz_z = cbuffer.data(gp_off + 26);

        #pragma omp simd aligned(cd_x, g_xyzz_x, g_xyzz_y, g_xyzz_z, g_yzz_x, g_yzz_xx, g_yzz_xy, g_yzz_xz, g_yzz_y, g_yzz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xyzz_x[k] = -g_yzz_x[k] * cd_x[k] + g_yzz_xx[k];

            g_xyzz_y[k] = -g_yzz_y[k] * cd_x[k] + g_yzz_xy[k];

            g_xyzz_z[k] = -g_yzz_z[k] * cd_x[k] + g_yzz_xz[k];
        }

        /// Set up 27-30 components of targeted buffer : cbuffer.data(

        auto g_xzzz_x = cbuffer.data(gp_off + 27);

        auto g_xzzz_y = cbuffer.data(gp_off + 28);

        auto g_xzzz_z = cbuffer.data(gp_off + 29);

        #pragma omp simd aligned(cd_x, g_xzzz_x, g_xzzz_y, g_xzzz_z, g_zzz_x, g_zzz_xx, g_zzz_xy, g_zzz_xz, g_zzz_y, g_zzz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xzzz_x[k] = -g_zzz_x[k] * cd_x[k] + g_zzz_xx[k];

            g_xzzz_y[k] = -g_zzz_y[k] * cd_x[k] + g_zzz_xy[k];

            g_xzzz_z[k] = -g_zzz_z[k] * cd_x[k] + g_zzz_xz[k];
        }

        /// Set up 30-33 components of targeted buffer : cbuffer.data(

        auto g_yyyy_x = cbuffer.data(gp_off + 30);

        auto g_yyyy_y = cbuffer.data(gp_off + 31);

        auto g_yyyy_z = cbuffer.data(gp_off + 32);

        #pragma omp simd aligned(cd_y, g_yyy_x, g_yyy_xy, g_yyy_y, g_yyy_yy, g_yyy_yz, g_yyy_z, g_yyyy_x, g_yyyy_y, g_yyyy_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yyyy_x[k] = -g_yyy_x[k] * cd_y[k] + g_yyy_xy[k];

            g_yyyy_y[k] = -g_yyy_y[k] * cd_y[k] + g_yyy_yy[k];

            g_yyyy_z[k] = -g_yyy_z[k] * cd_y[k] + g_yyy_yz[k];
        }

        /// Set up 33-36 components of targeted buffer : cbuffer.data(

        auto g_yyyz_x = cbuffer.data(gp_off + 33);

        auto g_yyyz_y = cbuffer.data(gp_off + 34);

        auto g_yyyz_z = cbuffer.data(gp_off + 35);

        #pragma omp simd aligned(cd_y, g_yyyz_x, g_yyyz_y, g_yyyz_z, g_yyz_x, g_yyz_xy, g_yyz_y, g_yyz_yy, g_yyz_yz, g_yyz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yyyz_x[k] = -g_yyz_x[k] * cd_y[k] + g_yyz_xy[k];

            g_yyyz_y[k] = -g_yyz_y[k] * cd_y[k] + g_yyz_yy[k];

            g_yyyz_z[k] = -g_yyz_z[k] * cd_y[k] + g_yyz_yz[k];
        }

        /// Set up 36-39 components of targeted buffer : cbuffer.data(

        auto g_yyzz_x = cbuffer.data(gp_off + 36);

        auto g_yyzz_y = cbuffer.data(gp_off + 37);

        auto g_yyzz_z = cbuffer.data(gp_off + 38);

        #pragma omp simd aligned(cd_y, g_yyzz_x, g_yyzz_y, g_yyzz_z, g_yzz_x, g_yzz_xy, g_yzz_y, g_yzz_yy, g_yzz_yz, g_yzz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yyzz_x[k] = -g_yzz_x[k] * cd_y[k] + g_yzz_xy[k];

            g_yyzz_y[k] = -g_yzz_y[k] * cd_y[k] + g_yzz_yy[k];

            g_yyzz_z[k] = -g_yzz_z[k] * cd_y[k] + g_yzz_yz[k];
        }

        /// Set up 39-42 components of targeted buffer : cbuffer.data(

        auto g_yzzz_x = cbuffer.data(gp_off + 39);

        auto g_yzzz_y = cbuffer.data(gp_off + 40);

        auto g_yzzz_z = cbuffer.data(gp_off + 41);

        #pragma omp simd aligned(cd_y, g_yzzz_x, g_yzzz_y, g_yzzz_z, g_zzz_x, g_zzz_xy, g_zzz_y, g_zzz_yy, g_zzz_yz, g_zzz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yzzz_x[k] = -g_zzz_x[k] * cd_y[k] + g_zzz_xy[k];

            g_yzzz_y[k] = -g_zzz_y[k] * cd_y[k] + g_zzz_yy[k];

            g_yzzz_z[k] = -g_zzz_z[k] * cd_y[k] + g_zzz_yz[k];
        }

        /// Set up 42-45 components of targeted buffer : cbuffer.data(

        auto g_zzzz_x = cbuffer.data(gp_off + 42);

        auto g_zzzz_y = cbuffer.data(gp_off + 43);

        auto g_zzzz_z = cbuffer.data(gp_off + 44);

        #pragma omp simd aligned(cd_z, g_zzz_x, g_zzz_xz, g_zzz_y, g_zzz_yz, g_zzz_z, g_zzz_zz, g_zzzz_x, g_zzzz_y, g_zzzz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_zzzz_x[k] = -g_zzz_x[k] * cd_z[k] + g_zzz_xz[k];

            g_zzzz_y[k] = -g_zzz_y[k] * cd_z[k] + g_zzz_yz[k];

            g_zzzz_z[k] = -g_zzz_z[k] * cd_z[k] + g_zzz_zz[k];
        }
    }
}

} // t3ceri namespace

