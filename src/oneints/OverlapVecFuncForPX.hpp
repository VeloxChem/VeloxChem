//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

namespace ovlvecfunc { // ovlvecfunc namespace

    // SIMD elementary functions for (P||P) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_x_x_s_0(double fx,
                               double pa_x,
                               double pb_x,
                               double s_0_0)
    {
        return s_0_0 * (0.5 * fx

                     + pa_x * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_y_s_0(double pa_x,
                               double pb_y,
                               double s_0_0)
    {
        return s_0_0 * (pa_x * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_z_s_0(double pa_x,
                               double pb_z,
                               double s_0_0)
    {
        return s_0_0 * (pa_x * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_x_s_0(double pa_y,
                               double pb_x,
                               double s_0_0)
    {
        return s_0_0 * (pa_y * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_y_s_0(double fx,
                               double pa_y,
                               double pb_y,
                               double s_0_0)
    {
        return s_0_0 * (0.5 * fx

                     + pa_y * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_z_s_0(double pa_y,
                               double pb_z,
                               double s_0_0)
    {
        return s_0_0 * (pa_y * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_x_s_0(double pa_z,
                               double pb_x,
                               double s_0_0)
    {
        return s_0_0 * (pa_z * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_y_s_0(double pa_z,
                               double pb_y,
                               double s_0_0)
    {
        return s_0_0 * (pa_z * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_z_s_0(double fx,
                               double pa_z,
                               double pb_z,
                               double s_0_0)
    {
        return s_0_0 * (0.5 * fx

                     + pa_z * pb_z);

    }

    // SIMD elementary functions for (P||D) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_x_xx_s_0(double fx,
                                double pa_x,
                                double pb_x,
                                double pb_xx,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx

                     + fx * pb_x

                     + pa_x * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_xy_s_0(double fx,
                                double pa_x,
                                double pb_xy,
                                double pb_y,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_y

                     + pa_x * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_xz_s_0(double fx,
                                double pa_x,
                                double pb_xz,
                                double pb_z,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_z

                     + pa_x * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_yy_s_0(double fx,
                                double pa_x,
                                double pb_yy,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx

                     + pa_x * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_yz_s_0(double pa_x,
                                double pb_yz,
                                double s_0_0)
    {
        return s_0_0 * (pa_x * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_zz_s_0(double fx,
                                double pa_x,
                                double pb_zz,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx

                     + pa_x * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_xx_s_0(double fx,
                                double pa_y,
                                double pb_xx,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * pa_y * fx

                     + pa_y * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_xy_s_0(double fx,
                                double pa_y,
                                double pb_x,
                                double pb_xy,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_x

                     + pa_y * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_xz_s_0(double pa_y,
                                double pb_xz,
                                double s_0_0)
    {
        return s_0_0 * (pa_y * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_yy_s_0(double fx,
                                double pa_y,
                                double pb_y,
                                double pb_yy,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * pa_y * fx

                     + fx * pb_y

                     + pa_y * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_yz_s_0(double fx,
                                double pa_y,
                                double pb_yz,
                                double pb_z,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_z

                     + pa_y * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_zz_s_0(double fx,
                                double pa_y,
                                double pb_zz,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * pa_y * fx

                     + pa_y * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_xx_s_0(double fx,
                                double pa_z,
                                double pb_xx,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * pa_z * fx

                     + pa_z * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_xy_s_0(double pa_z,
                                double pb_xy,
                                double s_0_0)
    {
        return s_0_0 * (pa_z * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_xz_s_0(double fx,
                                double pa_z,
                                double pb_x,
                                double pb_xz,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_x

                     + pa_z * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_yy_s_0(double fx,
                                double pa_z,
                                double pb_yy,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * pa_z * fx

                     + pa_z * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_yz_s_0(double fx,
                                double pa_z,
                                double pb_y,
                                double pb_yz,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_y

                     + pa_z * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_zz_s_0(double fx,
                                double pa_z,
                                double pb_z,
                                double pb_zz,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * pa_z * fx

                     + fx * pb_z

                     + pa_z * pb_zz);

    }

    // SIMD elementary functions for (D||P) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_xx_x_s_0(double fx,
                                double pa_x,
                                double pa_xx,
                                double pb_x,
                                double s_0_0)
    {
        return s_0_0 * (pa_x * fx

                     + 0.5 * fx * pb_x

                     + pa_xx * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_y_s_0(double fx,
                                double pa_xx,
                                double pb_y,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_y

                     + pa_xx * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_z_s_0(double fx,
                                double pa_xx,
                                double pb_z,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_z

                     + pa_xx * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_x_s_0(double fx,
                                double pa_xy,
                                double pa_y,
                                double pb_x,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pa_y

                     + pa_xy * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_y_s_0(double fx,
                                double pa_x,
                                double pa_xy,
                                double pb_y,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx

                     + pa_xy * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_z_s_0(double pa_xy,
                                double pb_z,
                                double s_0_0)
    {
        return s_0_0 * (pa_xy * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_x_s_0(double fx,
                                double pa_xz,
                                double pa_z,
                                double pb_x,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pa_z

                     + pa_xz * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_y_s_0(double pa_xz,
                                double pb_y,
                                double s_0_0)
    {
        return s_0_0 * (pa_xz * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_z_s_0(double fx,
                                double pa_x,
                                double pa_xz,
                                double pb_z,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx

                     + pa_xz * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_x_s_0(double fx,
                                double pa_yy,
                                double pb_x,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_x

                     + pa_yy * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_y_s_0(double fx,
                                double pa_y,
                                double pa_yy,
                                double pb_y,
                                double s_0_0)
    {
        return s_0_0 * (pa_y * fx

                     + 0.5 * fx * pb_y

                     + pa_yy * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_z_s_0(double fx,
                                double pa_yy,
                                double pb_z,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_z

                     + pa_yy * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_x_s_0(double pa_yz,
                                double pb_x,
                                double s_0_0)
    {
        return s_0_0 * (pa_yz * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_y_s_0(double fx,
                                double pa_yz,
                                double pa_z,
                                double pb_y,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pa_z

                     + pa_yz * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_z_s_0(double fx,
                                double pa_y,
                                double pa_yz,
                                double pb_z,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * pa_y * fx

                     + pa_yz * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_x_s_0(double fx,
                                double pa_zz,
                                double pb_x,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_x

                     + pa_zz * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_y_s_0(double fx,
                                double pa_zz,
                                double pb_y,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_y

                     + pa_zz * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_z_s_0(double fx,
                                double pa_z,
                                double pa_zz,
                                double pb_z,
                                double s_0_0)
    {
        return s_0_0 * (pa_z * fx

                     + 0.5 * fx * pb_z

                     + pa_zz * pb_z);

    }

    // SIMD elementary functions for (P||F) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_x_xxx_s_0(double fx,
                                 double pa_x,
                                 double pb_x,
                                 double pb_xx,
                                 double pb_xxx,
                                 double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx

                     + 1.5 * pa_x * pb_x * fx

                     + 1.5 * fx * pb_xx

                     + pa_x * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_xxy_s_0(double fx,
                                 double pa_x,
                                 double pb_xxy,
                                 double pb_xy,
                                 double pb_y,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx * pb_y

                     + fx * pb_xy

                     + pa_x * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_xxz_s_0(double fx,
                                 double pa_x,
                                 double pb_xxz,
                                 double pb_xz,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx * pb_z

                     + fx * pb_xz

                     + pa_x * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_xyy_s_0(double fx,
                                 double pa_x,
                                 double pb_x,
                                 double pb_xyy,
                                 double pb_yy,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_x * pb_x * fx

                     + 0.5 * fx * pb_yy

                     + pa_x * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_xyz_s_0(double fx,
                                 double pa_x,
                                 double pb_xyz,
                                 double pb_yz,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_yz

                     + pa_x * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_xzz_s_0(double fx,
                                 double pa_x,
                                 double pb_x,
                                 double pb_xzz,
                                 double pb_zz,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_x * pb_x * fx

                     + 0.5 * fx * pb_zz

                     + pa_x * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_yyy_s_0(double fx,
                                 double pa_x,
                                 double pb_y,
                                 double pb_yyy,
                                 double s_0_0)
    {
        return s_0_0 * (1.5 * pa_x * pb_y * fx

                     + pa_x * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_yyz_s_0(double fx,
                                 double pa_x,
                                 double pb_yyz,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx * pb_z

                     + pa_x * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_yzz_s_0(double fx,
                                 double pa_x,
                                 double pb_y,
                                 double pb_yzz,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * pb_y * fx

                     + pa_x * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_zzz_s_0(double fx,
                                 double pa_x,
                                 double pb_z,
                                 double pb_zzz,
                                 double s_0_0)
    {
        return s_0_0 * (1.5 * pa_x * pb_z * fx

                     + pa_x * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_xxx_s_0(double fx,
                                 double pa_y,
                                 double pb_x,
                                 double pb_xxx,
                                 double s_0_0)
    {
        return s_0_0 * (1.5 * pa_y * pb_x * fx

                     + pa_y * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_xxy_s_0(double fx,
                                 double pa_y,
                                 double pb_xx,
                                 double pb_xxy,
                                 double pb_y,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_y * fx * pb_y

                     + 0.5 * fx * pb_xx

                     + pa_y * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_xxz_s_0(double fx,
                                 double pa_y,
                                 double pb_xxz,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_y * fx * pb_z

                     + pa_y * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_xyy_s_0(double fx,
                                 double pa_y,
                                 double pb_x,
                                 double pb_xy,
                                 double pb_xyy,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_y * pb_x * fx

                     + fx * pb_xy

                     + pa_y * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_xyz_s_0(double fx,
                                 double pa_y,
                                 double pb_xyz,
                                 double pb_xz,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_xz

                     + pa_y * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_xzz_s_0(double fx,
                                 double pa_y,
                                 double pb_x,
                                 double pb_xzz,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_y * pb_x * fx

                     + pa_y * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_yyy_s_0(double fx,
                                 double pa_y,
                                 double pb_y,
                                 double pb_yy,
                                 double pb_yyy,
                                 double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx

                     + 1.5 * pa_y * pb_y * fx

                     + 1.5 * fx * pb_yy

                     + pa_y * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_yyz_s_0(double fx,
                                 double pa_y,
                                 double pb_yyz,
                                 double pb_yz,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_y * fx * pb_z

                     + fx * pb_yz

                     + pa_y * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_yzz_s_0(double fx,
                                 double pa_y,
                                 double pb_y,
                                 double pb_yzz,
                                 double pb_zz,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_y * pb_y * fx

                     + 0.5 * fx * pb_zz

                     + pa_y * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_zzz_s_0(double fx,
                                 double pa_y,
                                 double pb_z,
                                 double pb_zzz,
                                 double s_0_0)
    {
        return s_0_0 * (1.5 * pa_y * pb_z * fx

                     + pa_y * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_xxx_s_0(double fx,
                                 double pa_z,
                                 double pb_x,
                                 double pb_xxx,
                                 double s_0_0)
    {
        return s_0_0 * (1.5 * pa_z * pb_x * fx

                     + pa_z * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_xxy_s_0(double fx,
                                 double pa_z,
                                 double pb_xxy,
                                 double pb_y,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_z * fx * pb_y

                     + pa_z * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_xxz_s_0(double fx,
                                 double pa_z,
                                 double pb_xx,
                                 double pb_xxz,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_z * fx * pb_z

                     + 0.5 * fx * pb_xx

                     + pa_z * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_xyy_s_0(double fx,
                                 double pa_z,
                                 double pb_x,
                                 double pb_xyy,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_z * pb_x * fx

                     + pa_z * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_xyz_s_0(double fx,
                                 double pa_z,
                                 double pb_xy,
                                 double pb_xyz,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_xy

                     + pa_z * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_xzz_s_0(double fx,
                                 double pa_z,
                                 double pb_x,
                                 double pb_xz,
                                 double pb_xzz,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_z * pb_x * fx

                     + fx * pb_xz

                     + pa_z * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_yyy_s_0(double fx,
                                 double pa_z,
                                 double pb_y,
                                 double pb_yyy,
                                 double s_0_0)
    {
        return s_0_0 * (1.5 * pa_z * pb_y * fx

                     + pa_z * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_yyz_s_0(double fx,
                                 double pa_z,
                                 double pb_yy,
                                 double pb_yyz,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_z * fx * pb_z

                     + 0.5 * fx * pb_yy

                     + pa_z * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_yzz_s_0(double fx,
                                 double pa_z,
                                 double pb_y,
                                 double pb_yz,
                                 double pb_yzz,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_z * pb_y * fx

                     + fx * pb_yz

                     + pa_z * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_zzz_s_0(double fx,
                                 double pa_z,
                                 double pb_z,
                                 double pb_zz,
                                 double pb_zzz,
                                 double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx

                     + 1.5 * pa_z * pb_z * fx

                     + 1.5 * fx * pb_zz

                     + pa_z * pb_zzz);

    }

    // SIMD elementary functions for (F||P) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_x_s_0(double fx,
                                 double pa_x,
                                 double pa_xx,
                                 double pa_xxx,
                                 double pb_x,
                                 double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx

                     + 1.5 * pa_xx * fx

                     + 1.5 * pa_x * fx * pb_x

                     + pa_xxx * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_y_s_0(double fx,
                                 double pa_x,
                                 double pa_xxx,
                                 double pb_y,
                                 double s_0_0)
    {
        return s_0_0 * (1.5 * pa_x * fx * pb_y

                     + pa_xxx * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_z_s_0(double fx,
                                 double pa_x,
                                 double pa_xxx,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (1.5 * pa_x * fx * pb_z

                     + pa_xxx * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_x_s_0(double fx,
                                 double pa_xxy,
                                 double pa_xy,
                                 double pa_y,
                                 double pb_x,
                                 double s_0_0)
    {
        return s_0_0 * (pa_xy * fx

                     + 0.5 * fx * pa_y * pb_x

                     + pa_xxy * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_y_s_0(double fx,
                                 double pa_xx,
                                 double pa_xxy,
                                 double pa_y,
                                 double pb_y,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_xx * fx

                     + 0.5 * fx * pa_y * pb_y

                     + pa_xxy * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_z_s_0(double fx,
                                 double pa_xxy,
                                 double pa_y,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pa_y * pb_z

                     + pa_xxy * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_x_s_0(double fx,
                                 double pa_xxz,
                                 double pa_xz,
                                 double pa_z,
                                 double pb_x,
                                 double s_0_0)
    {
        return s_0_0 * (pa_xz * fx

                     + 0.5 * fx * pa_z * pb_x

                     + pa_xxz * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_y_s_0(double fx,
                                 double pa_xxz,
                                 double pa_z,
                                 double pb_y,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pa_z * pb_y

                     + pa_xxz * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_z_s_0(double fx,
                                 double pa_xx,
                                 double pa_xxz,
                                 double pa_z,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_xx * fx

                     + 0.5 * fx * pa_z * pb_z

                     + pa_xxz * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_x_s_0(double fx,
                                 double pa_x,
                                 double pa_xyy,
                                 double pa_yy,
                                 double pb_x,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * fx * pa_yy

                     + 0.5 * pa_x * fx * pb_x

                     + pa_xyy * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_y_s_0(double fx,
                                 double pa_x,
                                 double pa_xy,
                                 double pa_xyy,
                                 double pb_y,
                                 double s_0_0)
    {
        return s_0_0 * (pa_xy * fx

                     + 0.5 * pa_x * fx * pb_y

                     + pa_xyy * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_z_s_0(double fx,
                                 double pa_x,
                                 double pa_xyy,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx * pb_z

                     + pa_xyy * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_x_s_0(double fx,
                                 double pa_xyz,
                                 double pa_yz,
                                 double pb_x,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pa_yz

                     + pa_xyz * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_y_s_0(double fx,
                                 double pa_xyz,
                                 double pa_xz,
                                 double pb_y,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xz * fx

                     + pa_xyz * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_z_s_0(double fx,
                                 double pa_xy,
                                 double pa_xyz,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xy * fx

                     + pa_xyz * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_x_s_0(double fx,
                                 double pa_x,
                                 double pa_xzz,
                                 double pa_zz,
                                 double pb_x,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * fx * pa_zz

                     + 0.5 * pa_x * fx * pb_x

                     + pa_xzz * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_y_s_0(double fx,
                                 double pa_x,
                                 double pa_xzz,
                                 double pb_y,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx * pb_y

                     + pa_xzz * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_z_s_0(double fx,
                                 double pa_x,
                                 double pa_xz,
                                 double pa_xzz,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (pa_xz * fx

                     + 0.5 * pa_x * fx * pb_z

                     + pa_xzz * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_x_s_0(double fx,
                                 double pa_y,
                                 double pa_yyy,
                                 double pb_x,
                                 double s_0_0)
    {
        return s_0_0 * (1.5 * pa_y * fx * pb_x

                     + pa_yyy * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_y_s_0(double fx,
                                 double pa_y,
                                 double pa_yy,
                                 double pa_yyy,
                                 double pb_y,
                                 double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx

                     + 1.5 * pa_yy * fx

                     + 1.5 * pa_y * fx * pb_y

                     + pa_yyy * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_z_s_0(double fx,
                                 double pa_y,
                                 double pa_yyy,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (1.5 * pa_y * fx * pb_z

                     + pa_yyy * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_x_s_0(double fx,
                                 double pa_yyz,
                                 double pa_z,
                                 double pb_x,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pa_z * pb_x

                     + pa_yyz * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_y_s_0(double fx,
                                 double pa_yyz,
                                 double pa_yz,
                                 double pa_z,
                                 double pb_y,
                                 double s_0_0)
    {
        return s_0_0 * (pa_yz * fx

                     + 0.5 * fx * pa_z * pb_y

                     + pa_yyz * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_z_s_0(double fx,
                                 double pa_yy,
                                 double pa_yyz,
                                 double pa_z,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_yy * fx

                     + 0.5 * fx * pa_z * pb_z

                     + pa_yyz * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_x_s_0(double fx,
                                 double pa_y,
                                 double pa_yzz,
                                 double pb_x,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_y * fx * pb_x

                     + pa_yzz * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_y_s_0(double fx,
                                 double pa_y,
                                 double pa_yzz,
                                 double pa_zz,
                                 double pb_y,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * fx * pa_zz

                     + 0.5 * pa_y * fx * pb_y

                     + pa_yzz * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_z_s_0(double fx,
                                 double pa_y,
                                 double pa_yz,
                                 double pa_yzz,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (pa_yz * fx

                     + 0.5 * pa_y * fx * pb_z

                     + pa_yzz * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_x_s_0(double fx,
                                 double pa_z,
                                 double pa_zzz,
                                 double pb_x,
                                 double s_0_0)
    {
        return s_0_0 * (1.5 * pa_z * fx * pb_x

                     + pa_zzz * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_y_s_0(double fx,
                                 double pa_z,
                                 double pa_zzz,
                                 double pb_y,
                                 double s_0_0)
    {
        return s_0_0 * (1.5 * pa_z * fx * pb_y

                     + pa_zzz * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_z_s_0(double fx,
                                 double pa_z,
                                 double pa_zz,
                                 double pa_zzz,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx

                     + 1.5 * pa_zz * fx

                     + 1.5 * pa_z * fx * pb_z

                     + pa_zzz * pb_z);

    }

    // SIMD elementary functions for (P||G) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_x_xxxx_s_0(double fx,
                                  double pa_x,
                                  double pb_x,
                                  double pb_xx,
                                  double pb_xxx,
                                  double pb_xxxx,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx

                     + 3.0 * fx * fx * pb_x

                     + 3.0 * pa_x * pb_xx * fx

                     + 2.0 * fx * pb_xxx

                     + pa_x * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_xxxy_s_0(double fx,
                                  double pa_x,
                                  double pb_xxxy,
                                  double pb_xxy,
                                  double pb_xy,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_y

                     + 1.5 * pa_x * pb_xy * fx

                     + 1.5 * fx * pb_xxy

                     + pa_x * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_xxxz_s_0(double fx,
                                  double pa_x,
                                  double pb_xxxz,
                                  double pb_xxz,
                                  double pb_xz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_z

                     + 1.5 * pa_x * pb_xz * fx

                     + 1.5 * fx * pb_xxz

                     + pa_x * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_xxyy_s_0(double fx,
                                  double pa_x,
                                  double pb_x,
                                  double pb_xx,
                                  double pb_xxyy,
                                  double pb_xyy,
                                  double pb_yy,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_x * fx * fx

                     + 0.5 * fx * fx * pb_x

                     + 0.5 * pa_x * pb_xx * fx

                     + 0.5 * pa_x * fx * pb_yy

                     + fx * pb_xyy

                     + pa_x * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_xxyz_s_0(double fx,
                                  double pa_x,
                                  double pb_xxyz,
                                  double pb_xyz,
                                  double pb_yz,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx * pb_yz

                     + fx * pb_xyz

                     + pa_x * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_xxzz_s_0(double fx,
                                  double pa_x,
                                  double pb_x,
                                  double pb_xx,
                                  double pb_xxzz,
                                  double pb_xzz,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_x * fx * fx

                     + 0.5 * fx * fx * pb_x

                     + 0.5 * pa_x * pb_xx * fx

                     + 0.5 * pa_x * fx * pb_zz

                     + fx * pb_xzz

                     + pa_x * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_xyyy_s_0(double fx,
                                  double pa_x,
                                  double pb_xy,
                                  double pb_xyyy,
                                  double pb_y,
                                  double pb_yyy,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_y

                     + 1.5 * pa_x * pb_xy * fx

                     + 0.5 * fx * pb_yyy

                     + pa_x * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_xyyz_s_0(double fx,
                                  double pa_x,
                                  double pb_xyyz,
                                  double pb_xz,
                                  double pb_yyz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_z

                     + 0.5 * pa_x * pb_xz * fx

                     + 0.5 * fx * pb_yyz

                     + pa_x * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_xyzz_s_0(double fx,
                                  double pa_x,
                                  double pb_xy,
                                  double pb_xyzz,
                                  double pb_y,
                                  double pb_yzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_y

                     + 0.5 * pa_x * pb_xy * fx

                     + 0.5 * fx * pb_yzz

                     + pa_x * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_xzzz_s_0(double fx,
                                  double pa_x,
                                  double pb_xz,
                                  double pb_xzzz,
                                  double pb_z,
                                  double pb_zzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_z

                     + 1.5 * pa_x * pb_xz * fx

                     + 0.5 * fx * pb_zzz

                     + pa_x * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_yyyy_s_0(double fx,
                                  double pa_x,
                                  double pb_yy,
                                  double pb_yyyy,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx

                     + 3.0 * pa_x * pb_yy * fx

                     + pa_x * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_yyyz_s_0(double fx,
                                  double pa_x,
                                  double pb_yyyz,
                                  double pb_yz,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_x * pb_yz * fx

                     + pa_x * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_yyzz_s_0(double fx,
                                  double pa_x,
                                  double pb_yy,
                                  double pb_yyzz,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_x * fx * fx

                     + 0.5 * pa_x * pb_yy * fx

                     + 0.5 * pa_x * fx * pb_zz

                     + pa_x * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_yzzz_s_0(double fx,
                                  double pa_x,
                                  double pb_yz,
                                  double pb_yzzz,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_x * pb_yz * fx

                     + pa_x * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_zzzz_s_0(double fx,
                                  double pa_x,
                                  double pb_zz,
                                  double pb_zzzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx

                     + 3.0 * pa_x * pb_zz * fx

                     + pa_x * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_xxxx_s_0(double fx,
                                  double pa_y,
                                  double pb_xx,
                                  double pb_xxxx,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx

                     + 3.0 * pa_y * pb_xx * fx

                     + pa_y * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_xxxy_s_0(double fx,
                                  double pa_y,
                                  double pb_x,
                                  double pb_xxx,
                                  double pb_xxxy,
                                  double pb_xy,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_x

                     + 1.5 * pa_y * pb_xy * fx

                     + 0.5 * fx * pb_xxx

                     + pa_y * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_xxxz_s_0(double fx,
                                  double pa_y,
                                  double pb_xxxz,
                                  double pb_xz,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_y * pb_xz * fx

                     + pa_y * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_xxyy_s_0(double fx,
                                  double pa_y,
                                  double pb_xx,
                                  double pb_xxy,
                                  double pb_xxyy,
                                  double pb_y,
                                  double pb_yy,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_y * fx * fx

                     + 0.5 * fx * fx * pb_y

                     + 0.5 * pa_y * pb_xx * fx

                     + 0.5 * pa_y * fx * pb_yy

                     + fx * pb_xxy

                     + pa_y * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_xxyz_s_0(double fx,
                                  double pa_y,
                                  double pb_xxyz,
                                  double pb_xxz,
                                  double pb_yz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_z

                     + 0.5 * pa_y * fx * pb_yz

                     + 0.5 * fx * pb_xxz

                     + pa_y * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_xxzz_s_0(double fx,
                                  double pa_y,
                                  double pb_xx,
                                  double pb_xxzz,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_y * fx * fx

                     + 0.5 * pa_y * pb_xx * fx

                     + 0.5 * pa_y * fx * pb_zz

                     + pa_y * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_xyyy_s_0(double fx,
                                  double pa_y,
                                  double pb_x,
                                  double pb_xy,
                                  double pb_xyy,
                                  double pb_xyyy,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_x

                     + 1.5 * pa_y * pb_xy * fx

                     + 1.5 * fx * pb_xyy

                     + pa_y * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_xyyz_s_0(double fx,
                                  double pa_y,
                                  double pb_xyyz,
                                  double pb_xyz,
                                  double pb_xz,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_y * pb_xz * fx

                     + fx * pb_xyz

                     + pa_y * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_xyzz_s_0(double fx,
                                  double pa_y,
                                  double pb_x,
                                  double pb_xy,
                                  double pb_xyzz,
                                  double pb_xzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_x

                     + 0.5 * pa_y * pb_xy * fx

                     + 0.5 * fx * pb_xzz

                     + pa_y * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_xzzz_s_0(double fx,
                                  double pa_y,
                                  double pb_xz,
                                  double pb_xzzz,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_y * pb_xz * fx

                     + pa_y * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_yyyy_s_0(double fx,
                                  double pa_y,
                                  double pb_y,
                                  double pb_yy,
                                  double pb_yyy,
                                  double pb_yyyy,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx

                     + 3.0 * fx * fx * pb_y

                     + 3.0 * pa_y * pb_yy * fx

                     + 2.0 * fx * pb_yyy

                     + pa_y * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_yyyz_s_0(double fx,
                                  double pa_y,
                                  double pb_yyyz,
                                  double pb_yyz,
                                  double pb_yz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_z

                     + 1.5 * pa_y * pb_yz * fx

                     + 1.5 * fx * pb_yyz

                     + pa_y * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_yyzz_s_0(double fx,
                                  double pa_y,
                                  double pb_y,
                                  double pb_yy,
                                  double pb_yyzz,
                                  double pb_yzz,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_y * fx * fx

                     + 0.5 * fx * fx * pb_y

                     + 0.5 * pa_y * pb_yy * fx

                     + 0.5 * pa_y * fx * pb_zz

                     + fx * pb_yzz

                     + pa_y * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_yzzz_s_0(double fx,
                                  double pa_y,
                                  double pb_yz,
                                  double pb_yzzz,
                                  double pb_z,
                                  double pb_zzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_z

                     + 1.5 * pa_y * pb_yz * fx

                     + 0.5 * fx * pb_zzz

                     + pa_y * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_zzzz_s_0(double fx,
                                  double pa_y,
                                  double pb_zz,
                                  double pb_zzzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx

                     + 3.0 * pa_y * pb_zz * fx

                     + pa_y * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_xxxx_s_0(double fx,
                                  double pa_z,
                                  double pb_xx,
                                  double pb_xxxx,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_z * fx * fx

                     + 3.0 * pa_z * pb_xx * fx

                     + pa_z * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_xxxy_s_0(double fx,
                                  double pa_z,
                                  double pb_xxxy,
                                  double pb_xy,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_z * pb_xy * fx

                     + pa_z * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_xxxz_s_0(double fx,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xxx,
                                  double pb_xxxz,
                                  double pb_xz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_x

                     + 1.5 * pa_z * pb_xz * fx

                     + 0.5 * fx * pb_xxx

                     + pa_z * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_xxyy_s_0(double fx,
                                  double pa_z,
                                  double pb_xx,
                                  double pb_xxyy,
                                  double pb_yy,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_z * fx * fx

                     + 0.5 * pa_z * pb_xx * fx

                     + 0.5 * pa_z * fx * pb_yy

                     + pa_z * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_xxyz_s_0(double fx,
                                  double pa_z,
                                  double pb_xxy,
                                  double pb_xxyz,
                                  double pb_y,
                                  double pb_yz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_y

                     + 0.5 * pa_z * fx * pb_yz

                     + 0.5 * fx * pb_xxy

                     + pa_z * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_xxzz_s_0(double fx,
                                  double pa_z,
                                  double pb_xx,
                                  double pb_xxz,
                                  double pb_xxzz,
                                  double pb_z,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_z * fx * fx

                     + 0.5 * fx * fx * pb_z

                     + 0.5 * pa_z * pb_xx * fx

                     + 0.5 * pa_z * fx * pb_zz

                     + fx * pb_xxz

                     + pa_z * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_xyyy_s_0(double fx,
                                  double pa_z,
                                  double pb_xy,
                                  double pb_xyyy,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_z * pb_xy * fx

                     + pa_z * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_xyyz_s_0(double fx,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xyy,
                                  double pb_xyyz,
                                  double pb_xz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_x

                     + 0.5 * pa_z * pb_xz * fx

                     + 0.5 * fx * pb_xyy

                     + pa_z * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_xyzz_s_0(double fx,
                                  double pa_z,
                                  double pb_xy,
                                  double pb_xyz,
                                  double pb_xyzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_z * pb_xy * fx

                     + fx * pb_xyz

                     + pa_z * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_xzzz_s_0(double fx,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xz,
                                  double pb_xzz,
                                  double pb_xzzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_x

                     + 1.5 * pa_z * pb_xz * fx

                     + 1.5 * fx * pb_xzz

                     + pa_z * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_yyyy_s_0(double fx,
                                  double pa_z,
                                  double pb_yy,
                                  double pb_yyyy,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_z * fx * fx

                     + 3.0 * pa_z * pb_yy * fx

                     + pa_z * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_yyyz_s_0(double fx,
                                  double pa_z,
                                  double pb_y,
                                  double pb_yyy,
                                  double pb_yyyz,
                                  double pb_yz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_y

                     + 1.5 * pa_z * pb_yz * fx

                     + 0.5 * fx * pb_yyy

                     + pa_z * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_yyzz_s_0(double fx,
                                  double pa_z,
                                  double pb_yy,
                                  double pb_yyz,
                                  double pb_yyzz,
                                  double pb_z,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_z * fx * fx

                     + 0.5 * fx * fx * pb_z

                     + 0.5 * pa_z * pb_yy * fx

                     + 0.5 * pa_z * fx * pb_zz

                     + fx * pb_yyz

                     + pa_z * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_yzzz_s_0(double fx,
                                  double pa_z,
                                  double pb_y,
                                  double pb_yz,
                                  double pb_yzz,
                                  double pb_yzzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_y

                     + 1.5 * pa_z * pb_yz * fx

                     + 1.5 * fx * pb_yzz

                     + pa_z * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_zzzz_s_0(double fx,
                                  double pa_z,
                                  double pb_z,
                                  double pb_zz,
                                  double pb_zzz,
                                  double pb_zzzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_z * fx * fx

                     + 3.0 * fx * fx * pb_z

                     + 3.0 * pa_z * pb_zz * fx

                     + 2.0 * fx * pb_zzz

                     + pa_z * pb_zzzz);

    }

    // SIMD elementary functions for (G||P) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_x_s_0(double fx,
                                  double pa_x,
                                  double pa_xx,
                                  double pa_xxx,
                                  double pa_xxxx,
                                  double pb_x,
                                  double s_0_0)
    {
        return s_0_0 * (3.0 * pa_x * fx * fx

                     + 2.0 * pa_xxx * fx

                     + 0.75 * fx * fx * pb_x

                     + 3.0 * pa_xx * fx * pb_x

                     + pa_xxxx * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_y_s_0(double fx,
                                  double pa_xx,
                                  double pa_xxxx,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_y

                     + 3.0 * pa_xx * fx * pb_y

                     + pa_xxxx * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_z_s_0(double fx,
                                  double pa_xx,
                                  double pa_xxxx,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_z

                     + 3.0 * pa_xx * fx * pb_z

                     + pa_xxxx * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_x_s_0(double fx,
                                  double pa_xxxy,
                                  double pa_xxy,
                                  double pa_xy,
                                  double pa_y,
                                  double pb_x,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_y

                     + 1.5 * pa_xxy * fx

                     + 1.5 * pa_xy * fx * pb_x

                     + pa_xxxy * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_y_s_0(double fx,
                                  double pa_x,
                                  double pa_xxx,
                                  double pa_xxxy,
                                  double pa_xy,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx

                     + 0.5 * pa_xxx * fx

                     + 1.5 * pa_xy * fx * pb_y

                     + pa_xxxy * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_z_s_0(double fx,
                                  double pa_xxxy,
                                  double pa_xy,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xy * fx * pb_z

                     + pa_xxxy * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_x_s_0(double fx,
                                  double pa_xxxz,
                                  double pa_xxz,
                                  double pa_xz,
                                  double pa_z,
                                  double pb_x,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_z

                     + 1.5 * pa_xxz * fx

                     + 1.5 * pa_xz * fx * pb_x

                     + pa_xxxz * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_y_s_0(double fx,
                                  double pa_xxxz,
                                  double pa_xz,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xz * fx * pb_y

                     + pa_xxxz * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_z_s_0(double fx,
                                  double pa_x,
                                  double pa_xxx,
                                  double pa_xxxz,
                                  double pa_xz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx

                     + 0.5 * pa_xxx * fx

                     + 1.5 * pa_xz * fx * pb_z

                     + pa_xxxz * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_x_s_0(double fx,
                                  double pa_x,
                                  double pa_xx,
                                  double pa_xxyy,
                                  double pa_xyy,
                                  double pa_yy,
                                  double pb_x,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx * fx

                     + pa_xyy * fx

                     + 0.25 * fx * fx * pb_x

                     + 0.5 * pa_xx * fx * pb_x

                     + 0.5 * fx * pa_yy * pb_x

                     + pa_xxyy * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_y_s_0(double fx,
                                  double pa_xx,
                                  double pa_xxy,
                                  double pa_xxyy,
                                  double pa_y,
                                  double pa_yy,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * fx * fx * pa_y

                     + pa_xxy * fx

                     + 0.25 * fx * fx * pb_y

                     + 0.5 * pa_xx * fx * pb_y

                     + 0.5 * fx * pa_yy * pb_y

                     + pa_xxyy * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_z_s_0(double fx,
                                  double pa_xx,
                                  double pa_xxyy,
                                  double pa_yy,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_z

                     + 0.5 * pa_xx * fx * pb_z

                     + 0.5 * fx * pa_yy * pb_z

                     + pa_xxyy * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_x_s_0(double fx,
                                  double pa_xxyz,
                                  double pa_xyz,
                                  double pa_yz,
                                  double pb_x,
                                  double s_0_0)
    {
        return s_0_0 * (pa_xyz * fx

                     + 0.5 * fx * pa_yz * pb_x

                     + pa_xxyz * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_y_s_0(double fx,
                                  double pa_xxyz,
                                  double pa_xxz,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_z

                     + 0.5 * pa_xxz * fx

                     + 0.5 * fx * pa_yz * pb_y

                     + pa_xxyz * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_z_s_0(double fx,
                                  double pa_xxy,
                                  double pa_xxyz,
                                  double pa_y,
                                  double pa_yz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_y

                     + 0.5 * pa_xxy * fx

                     + 0.5 * fx * pa_yz * pb_z

                     + pa_xxyz * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_x_s_0(double fx,
                                  double pa_x,
                                  double pa_xx,
                                  double pa_xxzz,
                                  double pa_xzz,
                                  double pa_zz,
                                  double pb_x,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx * fx

                     + pa_xzz * fx

                     + 0.25 * fx * fx * pb_x

                     + 0.5 * pa_xx * fx * pb_x

                     + 0.5 * fx * pa_zz * pb_x

                     + pa_xxzz * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_y_s_0(double fx,
                                  double pa_xx,
                                  double pa_xxzz,
                                  double pa_zz,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_y

                     + 0.5 * pa_xx * fx * pb_y

                     + 0.5 * fx * pa_zz * pb_y

                     + pa_xxzz * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_z_s_0(double fx,
                                  double pa_xx,
                                  double pa_xxz,
                                  double pa_xxzz,
                                  double pa_z,
                                  double pa_zz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * fx * fx * pa_z

                     + pa_xxz * fx

                     + 0.25 * fx * fx * pb_z

                     + 0.5 * pa_xx * fx * pb_z

                     + 0.5 * fx * pa_zz * pb_z

                     + pa_xxzz * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_x_s_0(double fx,
                                  double pa_xy,
                                  double pa_xyyy,
                                  double pa_y,
                                  double pa_yyy,
                                  double pb_x,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_y

                     + 0.5 * fx * pa_yyy

                     + 1.5 * pa_xy * fx * pb_x

                     + pa_xyyy * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_y_s_0(double fx,
                                  double pa_x,
                                  double pa_xy,
                                  double pa_xyy,
                                  double pa_xyyy,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx

                     + 1.5 * pa_xyy * fx

                     + 1.5 * pa_xy * fx * pb_y

                     + pa_xyyy * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_z_s_0(double fx,
                                  double pa_xy,
                                  double pa_xyyy,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xy * fx * pb_z

                     + pa_xyyy * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_x_s_0(double fx,
                                  double pa_xyyz,
                                  double pa_xz,
                                  double pa_yyz,
                                  double pa_z,
                                  double pb_x,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_z

                     + 0.5 * fx * pa_yyz

                     + 0.5 * pa_xz * fx * pb_x

                     + pa_xyyz * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_y_s_0(double fx,
                                  double pa_xyyz,
                                  double pa_xyz,
                                  double pa_xz,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (pa_xyz * fx

                     + 0.5 * pa_xz * fx * pb_y

                     + pa_xyyz * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_z_s_0(double fx,
                                  double pa_x,
                                  double pa_xyy,
                                  double pa_xyyz,
                                  double pa_xz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_x * fx * fx

                     + 0.5 * pa_xyy * fx

                     + 0.5 * pa_xz * fx * pb_z

                     + pa_xyyz * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_x_s_0(double fx,
                                  double pa_xy,
                                  double pa_xyzz,
                                  double pa_y,
                                  double pa_yzz,
                                  double pb_x,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_y

                     + 0.5 * fx * pa_yzz

                     + 0.5 * pa_xy * fx * pb_x

                     + pa_xyzz * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_y_s_0(double fx,
                                  double pa_x,
                                  double pa_xy,
                                  double pa_xyzz,
                                  double pa_xzz,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_x * fx * fx

                     + 0.5 * pa_xzz * fx

                     + 0.5 * pa_xy * fx * pb_y

                     + pa_xyzz * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_z_s_0(double fx,
                                  double pa_xy,
                                  double pa_xyz,
                                  double pa_xyzz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (pa_xyz * fx

                     + 0.5 * pa_xy * fx * pb_z

                     + pa_xyzz * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_x_s_0(double fx,
                                  double pa_xz,
                                  double pa_xzzz,
                                  double pa_z,
                                  double pa_zzz,
                                  double pb_x,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_z

                     + 0.5 * fx * pa_zzz

                     + 1.5 * pa_xz * fx * pb_x

                     + pa_xzzz * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_y_s_0(double fx,
                                  double pa_xz,
                                  double pa_xzzz,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xz * fx * pb_y

                     + pa_xzzz * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_z_s_0(double fx,
                                  double pa_x,
                                  double pa_xz,
                                  double pa_xzz,
                                  double pa_xzzz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx

                     + 1.5 * pa_xzz * fx

                     + 1.5 * pa_xz * fx * pb_z

                     + pa_xzzz * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_x_s_0(double fx,
                                  double pa_yy,
                                  double pa_yyyy,
                                  double pb_x,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_x

                     + 3.0 * pa_yy * fx * pb_x

                     + pa_yyyy * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_y_s_0(double fx,
                                  double pa_y,
                                  double pa_yy,
                                  double pa_yyy,
                                  double pa_yyyy,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (3.0 * pa_y * fx * fx

                     + 2.0 * pa_yyy * fx

                     + 0.75 * fx * fx * pb_y

                     + 3.0 * pa_yy * fx * pb_y

                     + pa_yyyy * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_z_s_0(double fx,
                                  double pa_yy,
                                  double pa_yyyy,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_z

                     + 3.0 * pa_yy * fx * pb_z

                     + pa_yyyy * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_x_s_0(double fx,
                                  double pa_yyyz,
                                  double pa_yz,
                                  double pb_x,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_yz * fx * pb_x

                     + pa_yyyz * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_y_s_0(double fx,
                                  double pa_yyyz,
                                  double pa_yyz,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_z

                     + 1.5 * pa_yyz * fx

                     + 1.5 * pa_yz * fx * pb_y

                     + pa_yyyz * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_z_s_0(double fx,
                                  double pa_y,
                                  double pa_yyy,
                                  double pa_yyyz,
                                  double pa_yz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx

                     + 0.5 * pa_yyy * fx

                     + 1.5 * pa_yz * fx * pb_z

                     + pa_yyyz * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_x_s_0(double fx,
                                  double pa_yy,
                                  double pa_yyzz,
                                  double pa_zz,
                                  double pb_x,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_x

                     + 0.5 * pa_yy * fx * pb_x

                     + 0.5 * fx * pa_zz * pb_x

                     + pa_yyzz * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_y_s_0(double fx,
                                  double pa_y,
                                  double pa_yy,
                                  double pa_yyzz,
                                  double pa_yzz,
                                  double pa_zz,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_y * fx * fx

                     + pa_yzz * fx

                     + 0.25 * fx * fx * pb_y

                     + 0.5 * pa_yy * fx * pb_y

                     + 0.5 * fx * pa_zz * pb_y

                     + pa_yyzz * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_z_s_0(double fx,
                                  double pa_yy,
                                  double pa_yyz,
                                  double pa_yyzz,
                                  double pa_z,
                                  double pa_zz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * fx * fx * pa_z

                     + pa_yyz * fx

                     + 0.25 * fx * fx * pb_z

                     + 0.5 * pa_yy * fx * pb_z

                     + 0.5 * fx * pa_zz * pb_z

                     + pa_yyzz * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_x_s_0(double fx,
                                  double pa_yz,
                                  double pa_yzzz,
                                  double pb_x,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_yz * fx * pb_x

                     + pa_yzzz * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_y_s_0(double fx,
                                  double pa_yz,
                                  double pa_yzzz,
                                  double pa_z,
                                  double pa_zzz,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_z

                     + 0.5 * fx * pa_zzz

                     + 1.5 * pa_yz * fx * pb_y

                     + pa_yzzz * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_z_s_0(double fx,
                                  double pa_y,
                                  double pa_yz,
                                  double pa_yzz,
                                  double pa_yzzz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx

                     + 1.5 * pa_yzz * fx

                     + 1.5 * pa_yz * fx * pb_z

                     + pa_yzzz * pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_x_s_0(double fx,
                                  double pa_zz,
                                  double pa_zzzz,
                                  double pb_x,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_x

                     + 3.0 * pa_zz * fx * pb_x

                     + pa_zzzz * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_y_s_0(double fx,
                                  double pa_zz,
                                  double pa_zzzz,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_y

                     + 3.0 * pa_zz * fx * pb_y

                     + pa_zzzz * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_z_s_0(double fx,
                                  double pa_z,
                                  double pa_zz,
                                  double pa_zzz,
                                  double pa_zzzz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (3.0 * pa_z * fx * fx

                     + 2.0 * pa_zzz * fx

                     + 0.75 * fx * fx * pb_z

                     + 3.0 * pa_zz * fx * pb_z

                     + pa_zzzz * pb_z);

    }


} // ovlrecfunc namespace

