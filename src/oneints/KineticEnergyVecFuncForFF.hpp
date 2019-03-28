//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

namespace kinvecfunc { // kinvecfunc namespace

    // SIMD elementary functions for (F|T|F) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xxx_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx

                     + 2.25 * pa_xx * fx * fx

                     + 6.75 * pa_x * fx * fx * pb_x

                     + 2.25 * fx * fx * pb_xx

                     + 1.5 * pa_xxx * pb_x * fx

                     + 4.5 * pa_xx * fx * pb_xx

                     + 1.5 * pa_x * fx * pb_xxx

                     + pa_xxx * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xxx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double r_0_0)
    {
        return r_0_0 * (11.25 * fx * fx * fx * fz

                     - 2.25 * fx * fx * fz * fgb

                     - 2.25 * fx * fx * fz * fga

                     - 4.5 * pa_xx * fx * fz * fgb

                     + 18.0 * pa_xx * fz * fx * fx

                     + 54.0 * pa_x * fx * fx * fz * pb_x

                     - 4.5 * pa_x * fx * pb_x * fz * fgb

                     - 4.5 * pa_x * fz * fga * pb_x * fx

                     - 4.5 * fx * fz * fga * pb_xx

                     - 3.0 * pa_xxx * pb_x * fz * fgb

                     + 18.0 * fx * fx * fz * pb_xx

                     + 15.0 * pa_xxx * fz * pb_x * fx

                     + 45.0 * pa_xx * fz * fx * pb_xx

                     - 3.0 * pa_x * fz * fga * pb_xxx

                     + 15.0 * pa_x * fz * fx * pb_xxx

                     + 12.0 * pa_xxx * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xxy_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pb_xxy,
                                   double pb_xy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (2.25 * pa_x * fx * fx * pb_y

                     + 1.5 * fx * fx * pb_xy

                     + 0.5 * pa_xxx * fx * pb_y

                     + 3.0 * pa_xx * fx * pb_xy

                     + 1.5 * pa_x * fx * pb_xxy

                     + pa_xxx * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xxy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pb_xxy,
                                   double pb_xy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (18.0 * pa_x * fx * fx * fz * pb_y

                     - 1.5 * pa_x * fx * fz * fgb * pb_y

                     - 1.5 * pa_x * fz * fga * fx * pb_y

                     - 3.0 * fx * fz * fga * pb_xy

                     - pa_xxx * fz * fgb * pb_y

                     + 12.0 * fx * fx * fz * pb_xy

                     + 5.0 * pa_xxx * fz * fx * pb_y

                     + 30.0 * pa_xx * fz * fx * pb_xy

                     - 3.0 * pa_x * fz * fga * pb_xxy

                     + 15.0 * pa_x * fz * fx * pb_xxy

                     + 12.0 * pa_xxx * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xxz_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pb_xxz,
                                   double pb_xz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (2.25 * pa_x * fx * fx * pb_z

                     + 1.5 * fx * fx * pb_xz

                     + 0.5 * pa_xxx * fx * pb_z

                     + 3.0 * pa_xx * fx * pb_xz

                     + 1.5 * pa_x * fx * pb_xxz

                     + pa_xxx * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xxz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pb_xxz,
                                   double pb_xz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (18.0 * pa_x * fx * fx * fz * pb_z

                     - 1.5 * pa_x * fx * fz * fgb * pb_z

                     - 1.5 * pa_x * fz * fga * fx * pb_z

                     - 3.0 * fx * fz * fga * pb_xz

                     - pa_xxx * fz * fgb * pb_z

                     + 12.0 * fx * fx * fz * pb_xz

                     + 5.0 * pa_xxx * fz * fx * pb_z

                     + 30.0 * pa_xx * fz * fx * pb_xz

                     - 3.0 * pa_x * fz * fga * pb_xxz

                     + 15.0 * pa_x * fz * fx * pb_xxz

                     + 12.0 * pa_xxx * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xyy_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pb_x,
                                   double pb_xyy,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_xx * fx * fx

                     + 0.75 * pa_x * fx * fx * pb_x

                     + 0.75 * fx * fx * pb_yy

                     + 0.5 * pa_xxx * pb_x * fx

                     + 1.5 * pa_xx * fx * pb_yy

                     + 1.5 * pa_x * fx * pb_xyy

                     + pa_xxx * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pb_x,
                                   double pb_xyy,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     - 0.75 * fx * fx * fz * fga

                     - 1.5 * pa_xx * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * pa_xx * fz * fx * fx

                     - 1.5 * pa_x * fx * pb_x * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_x * fx

                     - 1.5 * fx * fz * fga * pb_yy

                     - pa_xxx * pb_x * fz * fgb

                     + 6.0 * pa_x * fz * fx * fx * pb_x

                     + 6.0 * fx * fx * fz * pb_yy

                     + 5.0 * pa_xxx * fz * pb_x * fx

                     + 15.0 * pa_xx * fz * fx * pb_yy

                     - 3.0 * pa_x * fz * fga * pb_xyy

                     + 15.0 * pa_x * fz * fx * pb_xyy

                     + 12.0 * pa_xxx * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xyz_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pb_xyz,
                                   double pb_yz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_yz

                     + 1.5 * pa_xx * fx * pb_yz

                     + 1.5 * pa_x * fx * pb_xyz

                     + pa_xxx * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xyz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pb_xyz,
                                   double pb_yz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fz * fga * pb_yz

                     + 6.0 * fx * fx * fz * pb_yz

                     + 15.0 * pa_xx * fz * fx * pb_yz

                     - 3.0 * pa_x * fz * fga * pb_xyz

                     + 15.0 * pa_x * fz * fx * pb_xyz

                     + 12.0 * pa_xxx * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pb_x,
                                   double pb_xzz,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_xx * fx * fx

                     + 0.75 * pa_x * fx * fx * pb_x

                     + 0.75 * fx * fx * pb_zz

                     + 0.5 * pa_xxx * pb_x * fx

                     + 1.5 * pa_xx * fx * pb_zz

                     + 1.5 * pa_x * fx * pb_xzz

                     + pa_xxx * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pb_x,
                                   double pb_xzz,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     - 0.75 * fx * fx * fz * fga

                     - 1.5 * pa_xx * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * pa_xx * fz * fx * fx

                     - 1.5 * pa_x * fx * pb_x * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_x * fx

                     - 1.5 * fx * fz * fga * pb_zz

                     - pa_xxx * pb_x * fz * fgb

                     + 6.0 * pa_x * fz * fx * fx * pb_x

                     + 6.0 * fx * fx * fz * pb_zz

                     + 5.0 * pa_xxx * fz * pb_x * fx

                     + 15.0 * pa_xx * fz * fx * pb_zz

                     - 3.0 * pa_x * fz * fga * pb_xzz

                     + 15.0 * pa_x * fz * fx * pb_xzz

                     + 12.0 * pa_xxx * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_yyy_s_0(double fx,
                                   double pa_x,
                                   double pa_xxx,
                                   double pb_y,
                                   double pb_yyy,
                                   double s_0_0)
    {
        return s_0_0 * (2.25 * pa_x * fx * fx * pb_y

                     + 1.5 * pa_xxx * pb_y * fx

                     + 1.5 * pa_x * fx * pb_yyy

                     + pa_xxx * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_yyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xxx,
                                   double pb_y,
                                   double pb_yyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_x * fx * pb_y * fz * fgb

                     - 4.5 * pa_x * fz * fga * pb_y * fx

                     - 3.0 * pa_xxx * pb_y * fz * fgb

                     + 18.0 * pa_x * fz * fx * fx * pb_y

                     + 15.0 * pa_xxx * fz * pb_y * fx

                     - 3.0 * pa_x * fz * fga * pb_yyy

                     + 15.0 * pa_x * fz * fx * pb_yyy

                     + 12.0 * pa_xxx * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_yyz_s_0(double fx,
                                   double pa_x,
                                   double pa_xxx,
                                   double pb_yyz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * pb_z

                     + 0.5 * pa_xxx * fx * pb_z

                     + 1.5 * pa_x * fx * pb_yyz

                     + pa_xxx * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_yyz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xxx,
                                   double pb_yyz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fz * fgb * pb_z

                     - 1.5 * pa_x * fz * fga * fx * pb_z

                     - pa_xxx * fz * fgb * pb_z

                     + 6.0 * pa_x * fz * fx * fx * pb_z

                     + 5.0 * pa_xxx * fz * fx * pb_z

                     - 3.0 * pa_x * fz * fga * pb_yyz

                     + 15.0 * pa_x * fz * fx * pb_yyz

                     + 12.0 * pa_xxx * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_yzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xxx,
                                   double pb_y,
                                   double pb_yzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * pb_y

                     + 0.5 * pa_xxx * pb_y * fx

                     + 1.5 * pa_x * fx * pb_yzz

                     + pa_xxx * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_yzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xxx,
                                   double pb_y,
                                   double pb_yzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * pb_y * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_y * fx

                     - pa_xxx * pb_y * fz * fgb

                     + 6.0 * pa_x * fz * fx * fx * pb_y

                     + 5.0 * pa_xxx * fz * pb_y * fx

                     - 3.0 * pa_x * fz * fga * pb_yzz

                     + 15.0 * pa_x * fz * fx * pb_yzz

                     + 12.0 * pa_xxx * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_zzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xxx,
                                   double pb_z,
                                   double pb_zzz,
                                   double s_0_0)
    {
        return s_0_0 * (2.25 * pa_x * fx * fx * pb_z

                     + 1.5 * pa_xxx * pb_z * fx

                     + 1.5 * pa_x * fx * pb_zzz

                     + pa_xxx * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_zzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xxx,
                                   double pb_z,
                                   double pb_zzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_x * fx * pb_z * fz * fgb

                     - 4.5 * pa_x * fz * fga * pb_z * fx

                     - 3.0 * pa_xxx * pb_z * fz * fgb

                     + 18.0 * pa_x * fz * fx * fx * pb_z

                     + 15.0 * pa_xxx * fz * pb_z * fx

                     - 3.0 * pa_x * fz * fga * pb_zzz

                     + 15.0 * pa_x * fz * fx * pb_zzz

                     + 12.0 * pa_xxx * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xxx_s_0(double fx,
                                   double pa_xxy,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xy * fx * fx

                     + 2.25 * fx * fx * pa_y * pb_x

                     + 1.5 * pa_xxy * pb_x * fx

                     + 3.0 * pa_xy * fx * pb_xx

                     + 0.5 * fx * pa_y * pb_xxx

                     + pa_xxy * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xxx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xxy,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xy * fx * fz * fgb

                     + 12.0 * pa_xy * fx * fx * fz

                     + 18.0 * fx * fx * fz * pa_y * pb_x

                     - 1.5 * fx * pa_y * pb_x * fz * fgb

                     - 1.5 * fz * fga * pa_y * pb_x * fx

                     - 3.0 * pa_xxy * pb_x * fz * fgb

                     + 15.0 * pa_xxy * fz * pb_x * fx

                     + 30.0 * pa_xy * fx * fz * pb_xx

                     - fz * fga * pa_y * pb_xxx

                     + 5.0 * fx * fz * pa_y * pb_xxx

                     + 12.0 * pa_xxy * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xxy_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxy,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxy,
                                   double pb_xy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.25 * pa_xx * fx * fx

                     + pa_x * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_y * pb_y

                     + 0.25 * fx * fx * pb_xx

                     + 0.5 * pa_xxy * fx * pb_y

                     + 0.5 * pa_xx * fx * pb_xx

                     + 2.0 * pa_xy * fx * pb_xy

                     + 0.5 * fx * pa_y * pb_xxy

                     + pa_xxy * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xxy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxy,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxy,
                                   double pb_xy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (2.25 * fx * fx * fx * fz

                     - 0.25 * fx * fx * fz * fgb

                     - 0.25 * fz * fga * fx * fx

                     - 0.5 * pa_xx * fx * fz * fgb

                     + 2.0 * pa_xx * fz * fx * fx

                     + 8.0 * pa_x * fx * fx * fz * pb_x

                     + 6.0 * fx * fx * fz * pa_y * pb_y

                     - 0.5 * fx * pa_y * fz * fgb * pb_y

                     - 0.5 * fz * fga * pa_y * fx * pb_y

                     - 0.5 * fz * fga * fx * pb_xx

                     - pa_xxy * fz * fgb * pb_y

                     + 2.0 * fx * fx * fz * pb_xx

                     + 5.0 * pa_xxy * fz * fx * pb_y

                     + 5.0 * pa_xx * fz * fx * pb_xx

                     + 20.0 * pa_xy * fx * fz * pb_xy

                     - fz * fga * pa_y * pb_xxy

                     + 5.0 * fx * fz * pa_y * pb_xxy

                     + 12.0 * pa_xxy * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xxz_s_0(double fx,
                                   double pa_xxy,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_xxz,
                                   double pb_xz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_y * pb_z

                     + 0.5 * pa_xxy * fx * pb_z

                     + 2.0 * pa_xy * fx * pb_xz

                     + 0.5 * fx * pa_y * pb_xxz

                     + pa_xxy * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xxz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xxy,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_xxz,
                                   double pb_xz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (6.0 * fx * fx * fz * pa_y * pb_z

                     - 0.5 * fx * pa_y * fz * fgb * pb_z

                     - 0.5 * fz * fga * pa_y * fx * pb_z

                     - pa_xxy * fz * fgb * pb_z

                     + 5.0 * pa_xxy * fz * fx * pb_z

                     + 20.0 * pa_xy * fx * fz * pb_xz

                     - fz * fga * pa_y * pb_xxz

                     + 5.0 * fx * fz * pa_y * pb_xxz

                     + 12.0 * pa_xxy * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xyy_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxy,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyy,
                                   double pb_y,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xy * fx * fx

                     + pa_x * fx * fx * pb_y

                     + 0.25 * fx * fx * pa_y * pb_x

                     + 0.5 * fx * fx * pb_xy

                     + 0.5 * pa_xxy * pb_x * fx

                     + pa_xx * fx * pb_xy

                     + pa_xy * fx * pb_yy

                     + 0.5 * fx * pa_y * pb_xyy

                     + pa_xxy * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxy,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyy,
                                   double pb_y,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_xy * fx * fz * fgb

                     + 4.0 * pa_xy * fx * fx * fz

                     + 8.0 * pa_x * fx * fx * fz * pb_y

                     - 0.5 * fx * pa_y * pb_x * fz * fgb

                     - 0.5 * fz * fga * pa_y * pb_x * fx

                     - fz * fga * fx * pb_xy

                     - pa_xxy * pb_x * fz * fgb

                     + 2.0 * fx * fx * fz * pa_y * pb_x

                     + 4.0 * fx * fx * fz * pb_xy

                     + 5.0 * pa_xxy * fz * pb_x * fx

                     + 10.0 * pa_xx * fz * fx * pb_xy

                     + 10.0 * pa_xy * fx * fz * pb_yy

                     - fz * fga * pa_y * pb_xyy

                     + 5.0 * fx * fz * pa_y * pb_xyy

                     + 12.0 * pa_xxy * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xyz_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxy,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_xyz,
                                   double pb_xz,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx * fx * pb_z

                     + 0.25 * fx * fx * pb_xz

                     + 0.5 * pa_xx * fx * pb_xz

                     + pa_xy * fx * pb_yz

                     + 0.5 * fx * pa_y * pb_xyz

                     + pa_xxy * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xyz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxy,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_xyz,
                                   double pb_xz,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (4.0 * pa_x * fx * fx * fz * pb_z

                     - 0.5 * fz * fga * fx * pb_xz

                     + 2.0 * fx * fx * fz * pb_xz

                     + 5.0 * pa_xx * fz * fx * pb_xz

                     + 10.0 * pa_xy * fx * fz * pb_yz

                     - fz * fga * pa_y * pb_xyz

                     + 5.0 * fx * fz * pa_y * pb_xyz

                     + 12.0 * pa_xxy * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xzz_s_0(double fx,
                                   double pa_xxy,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xzz,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xy * fx * fx

                     + 0.25 * fx * fx * pa_y * pb_x

                     + 0.5 * pa_xxy * pb_x * fx

                     + pa_xy * fx * pb_zz

                     + 0.5 * fx * pa_y * pb_xzz

                     + pa_xxy * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xxy,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xzz,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_xy * fx * fz * fgb

                     + 4.0 * pa_xy * fx * fx * fz

                     - 0.5 * fx * pa_y * pb_x * fz * fgb

                     - 0.5 * fz * fga * pa_y * pb_x * fx

                     - pa_xxy * pb_x * fz * fgb

                     + 2.0 * fx * fx * fz * pa_y * pb_x

                     + 5.0 * pa_xxy * fz * pb_x * fx

                     + 10.0 * pa_xy * fx * fz * pb_zz

                     - fz * fga * pa_y * pb_xzz

                     + 5.0 * fx * fz * pa_y * pb_xzz

                     + 12.0 * pa_xxy * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_yyy_s_0(double fx,
                                   double pa_xx,
                                   double pa_xxy,
                                   double pa_y,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_xx * fx * fx

                     + 0.75 * fx * fx * pa_y * pb_y

                     + 0.75 * fx * fx * pb_yy

                     + 1.5 * pa_xxy * pb_y * fx

                     + 1.5 * pa_xx * fx * pb_yy

                     + 0.5 * fx * pa_y * pb_yyy

                     + pa_xxy * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_yyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pa_xxy,
                                   double pa_y,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     - 0.75 * fz * fga * fx * fx

                     - 1.5 * pa_xx * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * pa_xx * fz * fx * fx

                     - 1.5 * fx * pa_y * pb_y * fz * fgb

                     - 1.5 * fz * fga * pa_y * pb_y * fx

                     - 1.5 * fz * fga * fx * pb_yy

                     - 3.0 * pa_xxy * pb_y * fz * fgb

                     + 6.0 * fx * fx * fz * pa_y * pb_y

                     + 6.0 * fx * fx * fz * pb_yy

                     + 15.0 * pa_xxy * fz * pb_y * fx

                     + 15.0 * pa_xx * fz * fx * pb_yy

                     - fz * fga * pa_y * pb_yyy

                     + 5.0 * fx * fz * pa_y * pb_yyy

                     + 12.0 * pa_xxy * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_yyz_s_0(double fx,
                                   double pa_xx,
                                   double pa_xxy,
                                   double pa_y,
                                   double pb_yyz,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_y * pb_z

                     + 0.5 * fx * fx * pb_yz

                     + 0.5 * pa_xxy * fx * pb_z

                     + pa_xx * fx * pb_yz

                     + 0.5 * fx * pa_y * pb_yyz

                     + pa_xxy * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_yyz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pa_xxy,
                                   double pa_y,
                                   double pb_yyz,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_y * fz * fgb * pb_z

                     - 0.5 * fz * fga * pa_y * fx * pb_z

                     - fz * fga * fx * pb_yz

                     - pa_xxy * fz * fgb * pb_z

                     + 2.0 * fx * fx * fz * pa_y * pb_z

                     + 4.0 * fx * fx * fz * pb_yz

                     + 5.0 * pa_xxy * fz * fx * pb_z

                     + 10.0 * pa_xx * fz * fx * pb_yz

                     - fz * fga * pa_y * pb_yyz

                     + 5.0 * fx * fz * pa_y * pb_yyz

                     + 12.0 * pa_xxy * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_yzz_s_0(double fx,
                                   double pa_xx,
                                   double pa_xxy,
                                   double pa_y,
                                   double pb_y,
                                   double pb_yzz,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx

                     + 0.25 * pa_xx * fx * fx

                     + 0.25 * fx * fx * pa_y * pb_y

                     + 0.25 * fx * fx * pb_zz

                     + 0.5 * pa_xxy * pb_y * fx

                     + 0.5 * pa_xx * fx * pb_zz

                     + 0.5 * fx * pa_y * pb_yzz

                     + pa_xxy * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_yzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pa_xxy,
                                   double pa_y,
                                   double pb_y,
                                   double pb_yzz,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.25 * fx * fx * fz * fgb

                     - 0.25 * fz * fga * fx * fx

                     - 0.5 * pa_xx * fx * fz * fgb

                     + 0.75 * fx * fx * fx * fz

                     + 2.0 * pa_xx * fz * fx * fx

                     - 0.5 * fx * pa_y * pb_y * fz * fgb

                     - 0.5 * fz * fga * pa_y * pb_y * fx

                     - 0.5 * fz * fga * fx * pb_zz

                     - pa_xxy * pb_y * fz * fgb

                     + 2.0 * fx * fx * fz * pa_y * pb_y

                     + 2.0 * fx * fx * fz * pb_zz

                     + 5.0 * pa_xxy * fz * pb_y * fx

                     + 5.0 * pa_xx * fz * fx * pb_zz

                     - fz * fga * pa_y * pb_yzz

                     + 5.0 * fx * fz * pa_y * pb_yzz

                     + 12.0 * pa_xxy * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_zzz_s_0(double fx,
                                   double pa_xxy,
                                   double pa_y,
                                   double pb_z,
                                   double pb_zzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_y * pb_z

                     + 1.5 * pa_xxy * pb_z * fx

                     + 0.5 * fx * pa_y * pb_zzz

                     + pa_xxy * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_zzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xxy,
                                   double pa_y,
                                   double pb_z,
                                   double pb_zzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pa_y * pb_z * fz * fgb

                     - 1.5 * fz * fga * pa_y * pb_z * fx

                     - 3.0 * pa_xxy * pb_z * fz * fgb

                     + 6.0 * fx * fx * fz * pa_y * pb_z

                     + 15.0 * pa_xxy * fz * pb_z * fx

                     - fz * fga * pa_y * pb_zzz

                     + 5.0 * fx * fz * pa_y * pb_zzz

                     + 12.0 * pa_xxy * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xxx_s_0(double fx,
                                   double pa_xxz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xz * fx * fx

                     + 2.25 * fx * fx * pa_z * pb_x

                     + 1.5 * pa_xxz * pb_x * fx

                     + 3.0 * pa_xz * fx * pb_xx

                     + 0.5 * fx * pa_z * pb_xxx

                     + pa_xxz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xxx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xxz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xz * fx * fz * fgb

                     + 12.0 * pa_xz * fx * fx * fz

                     + 18.0 * fx * fx * fz * pa_z * pb_x

                     - 1.5 * fx * pa_z * pb_x * fz * fgb

                     - 1.5 * fz * fga * pa_z * pb_x * fx

                     - 3.0 * pa_xxz * pb_x * fz * fgb

                     + 15.0 * pa_xxz * fz * pb_x * fx

                     + 30.0 * pa_xz * fx * fz * pb_xx

                     - fz * fga * pa_z * pb_xxx

                     + 5.0 * fx * fz * pa_z * pb_xxx

                     + 12.0 * pa_xxz * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xxy_s_0(double fx,
                                   double pa_xxz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_xxy,
                                   double pb_xy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_z * pb_y

                     + 0.5 * pa_xxz * fx * pb_y

                     + 2.0 * pa_xz * fx * pb_xy

                     + 0.5 * fx * pa_z * pb_xxy

                     + pa_xxz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xxy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xxz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_xxy,
                                   double pb_xy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (6.0 * fx * fx * fz * pa_z * pb_y

                     - 0.5 * fx * pa_z * fz * fgb * pb_y

                     - 0.5 * fz * fga * pa_z * fx * pb_y

                     - pa_xxz * fz * fgb * pb_y

                     + 5.0 * pa_xxz * fz * fx * pb_y

                     + 20.0 * pa_xz * fx * fz * pb_xy

                     - fz * fga * pa_z * pb_xxy

                     + 5.0 * fx * fz * pa_z * pb_xxy

                     + 12.0 * pa_xxz * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xxz_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxz,
                                   double pb_xz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.25 * pa_xx * fx * fx

                     + pa_x * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_z * pb_z

                     + 0.25 * fx * fx * pb_xx

                     + 0.5 * pa_xxz * fx * pb_z

                     + 0.5 * pa_xx * fx * pb_xx

                     + 2.0 * pa_xz * fx * pb_xz

                     + 0.5 * fx * pa_z * pb_xxz

                     + pa_xxz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xxz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxz,
                                   double pb_xz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (2.25 * fx * fx * fx * fz

                     - 0.25 * fx * fx * fz * fgb

                     - 0.25 * fz * fga * fx * fx

                     - 0.5 * pa_xx * fx * fz * fgb

                     + 2.0 * pa_xx * fz * fx * fx

                     + 8.0 * pa_x * fx * fx * fz * pb_x

                     + 6.0 * fx * fx * fz * pa_z * pb_z

                     - 0.5 * fx * pa_z * fz * fgb * pb_z

                     - 0.5 * fz * fga * pa_z * fx * pb_z

                     - 0.5 * fz * fga * fx * pb_xx

                     - pa_xxz * fz * fgb * pb_z

                     + 2.0 * fx * fx * fz * pb_xx

                     + 5.0 * pa_xxz * fz * fx * pb_z

                     + 5.0 * pa_xx * fz * fx * pb_xx

                     + 20.0 * pa_xz * fx * fz * pb_xz

                     - fz * fga * pa_z * pb_xxz

                     + 5.0 * fx * fz * pa_z * pb_xxz

                     + 12.0 * pa_xxz * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xyy_s_0(double fx,
                                   double pa_xxz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xyy,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xz * fx * fx

                     + 0.25 * fx * fx * pa_z * pb_x

                     + 0.5 * pa_xxz * pb_x * fx

                     + pa_xz * fx * pb_yy

                     + 0.5 * fx * pa_z * pb_xyy

                     + pa_xxz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xxz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xyy,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_xz * fx * fz * fgb

                     + 4.0 * pa_xz * fx * fx * fz

                     - 0.5 * fx * pa_z * pb_x * fz * fgb

                     - 0.5 * fz * fga * pa_z * pb_x * fx

                     - pa_xxz * pb_x * fz * fgb

                     + 2.0 * fx * fx * fz * pa_z * pb_x

                     + 5.0 * pa_xxz * fz * pb_x * fx

                     + 10.0 * pa_xz * fx * fz * pb_yy

                     - fz * fga * pa_z * pb_xyy

                     + 5.0 * fx * fz * pa_z * pb_xyy

                     + 12.0 * pa_xxz * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xyz_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_xy,
                                   double pb_xyz,
                                   double pb_y,
                                   double pb_yz,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx * fx * pb_y

                     + 0.25 * fx * fx * pb_xy

                     + 0.5 * pa_xx * fx * pb_xy

                     + pa_xz * fx * pb_yz

                     + 0.5 * fx * pa_z * pb_xyz

                     + pa_xxz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xyz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_xy,
                                   double pb_xyz,
                                   double pb_y,
                                   double pb_yz,
                                   double r_0_0)
    {
        return r_0_0 * (4.0 * pa_x * fx * fx * fz * pb_y

                     - 0.5 * fz * fga * fx * pb_xy

                     + 2.0 * fx * fx * fz * pb_xy

                     + 5.0 * pa_xx * fz * fx * pb_xy

                     + 10.0 * pa_xz * fx * fz * pb_yz

                     - fz * fga * pa_z * pb_xyz

                     + 5.0 * fx * fz * pa_z * pb_xyz

                     + 12.0 * pa_xxz * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_xzz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xz * fx * fx

                     + pa_x * fx * fx * pb_z

                     + 0.25 * fx * fx * pa_z * pb_x

                     + 0.5 * fx * fx * pb_xz

                     + 0.5 * pa_xxz * pb_x * fx

                     + pa_xx * fx * pb_xz

                     + pa_xz * fx * pb_zz

                     + 0.5 * fx * pa_z * pb_xzz

                     + pa_xxz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_xzz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_xz * fx * fz * fgb

                     + 4.0 * pa_xz * fx * fx * fz

                     + 8.0 * pa_x * fx * fx * fz * pb_z

                     - 0.5 * fx * pa_z * pb_x * fz * fgb

                     - 0.5 * fz * fga * pa_z * pb_x * fx

                     - fz * fga * fx * pb_xz

                     - pa_xxz * pb_x * fz * fgb

                     + 2.0 * fx * fx * fz * pa_z * pb_x

                     + 4.0 * fx * fx * fz * pb_xz

                     + 5.0 * pa_xxz * fz * pb_x * fx

                     + 10.0 * pa_xx * fz * fx * pb_xz

                     + 10.0 * pa_xz * fx * fz * pb_zz

                     - fz * fga * pa_z * pb_xzz

                     + 5.0 * fx * fz * pa_z * pb_xzz

                     + 12.0 * pa_xxz * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_yyy_s_0(double fx,
                                   double pa_xxz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yyy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_z * pb_y

                     + 1.5 * pa_xxz * pb_y * fx

                     + 0.5 * fx * pa_z * pb_yyy

                     + pa_xxz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_yyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xxz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pa_z * pb_y * fz * fgb

                     - 1.5 * fz * fga * pa_z * pb_y * fx

                     - 3.0 * pa_xxz * pb_y * fz * fgb

                     + 6.0 * fx * fx * fz * pa_z * pb_y

                     + 15.0 * pa_xxz * fz * pb_y * fx

                     - fz * fga * pa_z * pb_yyy

                     + 5.0 * fx * fz * pa_z * pb_yyy

                     + 12.0 * pa_xxz * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_yyz_s_0(double fx,
                                   double pa_xx,
                                   double pa_xxz,
                                   double pa_z,
                                   double pb_yy,
                                   double pb_yyz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx

                     + 0.25 * pa_xx * fx * fx

                     + 0.25 * fx * fx * pa_z * pb_z

                     + 0.25 * fx * fx * pb_yy

                     + 0.5 * pa_xxz * fx * pb_z

                     + 0.5 * pa_xx * fx * pb_yy

                     + 0.5 * fx * pa_z * pb_yyz

                     + pa_xxz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_yyz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pa_xxz,
                                   double pa_z,
                                   double pb_yy,
                                   double pb_yyz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.25 * fx * fx * fz * fgb

                     - 0.25 * fz * fga * fx * fx

                     - 0.5 * pa_xx * fx * fz * fgb

                     + 0.75 * fx * fx * fx * fz

                     + 2.0 * pa_xx * fz * fx * fx

                     - 0.5 * fx * pa_z * fz * fgb * pb_z

                     - 0.5 * fz * fga * pa_z * fx * pb_z

                     - 0.5 * fz * fga * fx * pb_yy

                     - pa_xxz * fz * fgb * pb_z

                     + 2.0 * fx * fx * fz * pa_z * pb_z

                     + 2.0 * fx * fx * fz * pb_yy

                     + 5.0 * pa_xxz * fz * fx * pb_z

                     + 5.0 * pa_xx * fz * fx * pb_yy

                     - fz * fga * pa_z * pb_yyz

                     + 5.0 * fx * fz * pa_z * pb_yyz

                     + 12.0 * pa_xxz * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_yzz_s_0(double fx,
                                   double pa_xx,
                                   double pa_xxz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_yzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_z * pb_y

                     + 0.5 * fx * fx * pb_yz

                     + 0.5 * pa_xxz * pb_y * fx

                     + pa_xx * fx * pb_yz

                     + 0.5 * fx * pa_z * pb_yzz

                     + pa_xxz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_yzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pa_xxz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_yzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_z * pb_y * fz * fgb

                     - 0.5 * fz * fga * pa_z * pb_y * fx

                     - fz * fga * fx * pb_yz

                     - pa_xxz * pb_y * fz * fgb

                     + 2.0 * fx * fx * fz * pa_z * pb_y

                     + 4.0 * fx * fx * fz * pb_yz

                     + 5.0 * pa_xxz * fz * pb_y * fx

                     + 10.0 * pa_xx * fz * fx * pb_yz

                     - fz * fga * pa_z * pb_yzz

                     + 5.0 * fx * fz * pa_z * pb_yzz

                     + 12.0 * pa_xxz * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_zzz_s_0(double fx,
                                   double pa_xx,
                                   double pa_xxz,
                                   double pa_z,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_xx * fx * fx

                     + 0.75 * fx * fx * pa_z * pb_z

                     + 0.75 * fx * fx * pb_zz

                     + 1.5 * pa_xxz * pb_z * fx

                     + 1.5 * pa_xx * fx * pb_zz

                     + 0.5 * fx * pa_z * pb_zzz

                     + pa_xxz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_zzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pa_xxz,
                                   double pa_z,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     - 0.75 * fz * fga * fx * fx

                     - 1.5 * pa_xx * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * pa_xx * fz * fx * fx

                     - 1.5 * fx * pa_z * pb_z * fz * fgb

                     - 1.5 * fz * fga * pa_z * pb_z * fx

                     - 1.5 * fz * fga * fx * pb_zz

                     - 3.0 * pa_xxz * pb_z * fz * fgb

                     + 6.0 * fx * fx * fz * pa_z * pb_z

                     + 6.0 * fx * fx * fz * pb_zz

                     + 15.0 * pa_xxz * fz * pb_z * fx

                     + 15.0 * pa_xx * fz * fx * pb_zz

                     - fz * fga * pa_z * pb_zzz

                     + 5.0 * fx * fz * pa_z * pb_zzz

                     + 12.0 * pa_xxz * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xxx_s_0(double fx,
                                   double pa_x,
                                   double pa_xyy,
                                   double pa_yy,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * fx * fx * pa_yy

                     + 0.75 * pa_x * fx * fx * pb_x

                     + 0.75 * fx * fx * pb_xx

                     + 1.5 * pa_xyy * pb_x * fx

                     + 1.5 * fx * pa_yy * pb_xx

                     + 0.5 * pa_x * fx * pb_xxx

                     + pa_xyy * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xxx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xyy,
                                   double pa_yy,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     - 0.75 * fx * fx * fz * fga

                     - 1.5 * fx * pa_yy * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * fx * fx * pa_yy * fz

                     - 1.5 * pa_x * fx * pb_x * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_x * fx

                     - 1.5 * fx * fz * fga * pb_xx

                     - 3.0 * pa_xyy * pb_x * fz * fgb

                     + 6.0 * pa_x * fz * fx * fx * pb_x

                     + 6.0 * fx * fx * fz * pb_xx

                     + 15.0 * pa_xyy * fz * pb_x * fx

                     + 15.0 * fx * pa_yy * fz * pb_xx

                     - pa_x * fz * fga * pb_xxx

                     + 5.0 * pa_x * fz * fx * pb_xxx

                     + 12.0 * pa_xyy * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xxy_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxy,
                                   double pb_xy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xy * fx * fx

                     + fx * fx * pa_y * pb_x

                     + 0.25 * pa_x * fx * fx * pb_y

                     + 0.5 * fx * fx * pb_xy

                     + 0.5 * pa_xyy * fx * pb_y

                     + pa_xy * fx * pb_xx

                     + fx * pa_yy * pb_xy

                     + 0.5 * pa_x * fx * pb_xxy

                     + pa_xyy * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xxy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxy,
                                   double pb_xy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_xy * fx * fz * fgb

                     + 4.0 * pa_xy * fz * fx * fx

                     + 8.0 * fx * fx * pa_y * fz * pb_x

                     - 0.5 * pa_x * fx * fz * fgb * pb_y

                     - 0.5 * pa_x * fz * fga * fx * pb_y

                     - fx * fz * fga * pb_xy

                     - pa_xyy * fz * fgb * pb_y

                     + 2.0 * pa_x * fz * fx * fx * pb_y

                     + 4.0 * fx * fx * fz * pb_xy

                     + 5.0 * pa_xyy * fz * fx * pb_y

                     + 10.0 * pa_xy * fz * fx * pb_xx

                     + 10.0 * fx * pa_yy * fz * pb_xy

                     - pa_x * fz * fga * pb_xxy

                     + 5.0 * pa_x * fz * fx * pb_xxy

                     + 12.0 * pa_xyy * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xxz_s_0(double fx,
                                   double pa_x,
                                   double pa_xyy,
                                   double pa_yy,
                                   double pb_xxz,
                                   double pb_xz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_x * fx * fx * pb_z

                     + 0.5 * fx * fx * pb_xz

                     + 0.5 * pa_xyy * fx * pb_z

                     + fx * pa_yy * pb_xz

                     + 0.5 * pa_x * fx * pb_xxz

                     + pa_xyy * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xxz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xyy,
                                   double pa_yy,
                                   double pb_xxz,
                                   double pb_xz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_x * fx * fz * fgb * pb_z

                     - 0.5 * pa_x * fz * fga * fx * pb_z

                     - fx * fz * fga * pb_xz

                     - pa_xyy * fz * fgb * pb_z

                     + 2.0 * pa_x * fz * fx * fx * pb_z

                     + 4.0 * fx * fx * fz * pb_xz

                     + 5.0 * pa_xyy * fz * fx * pb_z

                     + 10.0 * fx * pa_yy * fz * pb_xz

                     - pa_x * fz * fga * pb_xxz

                     + 5.0 * pa_x * fz * fx * pb_xxz

                     + 12.0 * pa_xyy * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xyy_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyy,
                                   double pb_y,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_x * fx * fx * pb_x

                     + 0.25 * fx * fx * pa_yy

                     + fx * fx * pa_y * pb_y

                     + 0.25 * fx * fx * pb_yy

                     + 0.5 * pa_xyy * pb_x * fx

                     + 2.0 * pa_xy * fx * pb_xy

                     + 0.5 * fx * pa_yy * pb_yy

                     + 0.5 * pa_x * fx * pb_xyy

                     + pa_xyy * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyy,
                                   double pb_y,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (2.25 * fx * fx * fx * fz

                     - 0.25 * fx * fx * fz * fgb

                     - 0.25 * fx * fx * fz * fga

                     - 0.5 * fx * pa_yy * fz * fgb

                     + 6.0 * pa_x * fx * fx * fz * pb_x

                     + 2.0 * fx * fx * pa_yy * fz

                     + 8.0 * fx * fx * pa_y * fz * pb_y

                     - 0.5 * pa_x * fx * pb_x * fz * fgb

                     - 0.5 * pa_x * fz * fga * pb_x * fx

                     - 0.5 * fx * fz * fga * pb_yy

                     - pa_xyy * pb_x * fz * fgb

                     + 2.0 * fx * fx * fz * pb_yy

                     + 5.0 * pa_xyy * fz * pb_x * fx

                     + 20.0 * pa_xy * fz * fx * pb_xy

                     + 5.0 * fx * pa_yy * fz * pb_yy

                     - pa_x * fz * fga * pb_xyy

                     + 5.0 * pa_x * fz * fx * pb_xyy

                     + 12.0 * pa_xyy * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xyz_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_xyz,
                                   double pb_xz,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * fx * fx * pa_y * pb_z

                     + 0.25 * fx * fx * pb_yz

                     + pa_xy * fx * pb_xz

                     + 0.5 * fx * pa_yy * pb_yz

                     + 0.5 * pa_x * fx * pb_xyz

                     + pa_xyy * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xyz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_xyz,
                                   double pb_xz,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (4.0 * fx * fx * pa_y * fz * pb_z

                     - 0.5 * fx * fz * fga * pb_yz

                     + 2.0 * fx * fx * fz * pb_yz

                     + 10.0 * pa_xy * fz * fx * pb_xz

                     + 5.0 * fx * pa_yy * fz * pb_yz

                     - pa_x * fz * fga * pb_xyz

                     + 5.0 * pa_x * fz * fx * pb_xyz

                     + 12.0 * pa_xyy * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xyy,
                                   double pa_yy,
                                   double pb_x,
                                   double pb_xzz,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx

                     + 0.25 * fx * fx * pa_yy

                     + 0.25 * pa_x * fx * fx * pb_x

                     + 0.25 * fx * fx * pb_zz

                     + 0.5 * pa_xyy * pb_x * fx

                     + 0.5 * fx * pa_yy * pb_zz

                     + 0.5 * pa_x * fx * pb_xzz

                     + pa_xyy * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xyy,
                                   double pa_yy,
                                   double pb_x,
                                   double pb_xzz,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.25 * fx * fx * fz * fgb

                     - 0.25 * fx * fx * fz * fga

                     - 0.5 * fx * pa_yy * fz * fgb

                     + 0.75 * fx * fx * fx * fz

                     + 2.0 * fx * fx * pa_yy * fz

                     - 0.5 * pa_x * fx * pb_x * fz * fgb

                     - 0.5 * pa_x * fz * fga * pb_x * fx

                     - 0.5 * fx * fz * fga * pb_zz

                     - pa_xyy * pb_x * fz * fgb

                     + 2.0 * pa_x * fz * fx * fx * pb_x

                     + 2.0 * fx * fx * fz * pb_zz

                     + 5.0 * pa_xyy * fz * pb_x * fx

                     + 5.0 * fx * pa_yy * fz * pb_zz

                     - pa_x * fz * fga * pb_xzz

                     + 5.0 * pa_x * fz * fx * pb_xzz

                     + 12.0 * pa_xyy * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_yyy_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xy * fx * fx

                     + 2.25 * pa_x * fx * fx * pb_y

                     + 1.5 * pa_xyy * pb_y * fx

                     + 3.0 * pa_xy * fx * pb_yy

                     + 0.5 * pa_x * fx * pb_yyy

                     + pa_xyy * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_yyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xy * fx * fz * fgb

                     + 12.0 * pa_xy * fz * fx * fx

                     + 18.0 * pa_x * fx * fx * fz * pb_y

                     - 1.5 * pa_x * fx * pb_y * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_y * fx

                     - 3.0 * pa_xyy * pb_y * fz * fgb

                     + 15.0 * pa_xyy * fz * pb_y * fx

                     + 30.0 * pa_xy * fz * fx * pb_yy

                     - pa_x * fz * fga * pb_yyy

                     + 5.0 * pa_x * fz * fx * pb_yyy

                     + 12.0 * pa_xyy * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_yyz_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pb_yyz,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * pb_z

                     + 0.5 * pa_xyy * fx * pb_z

                     + 2.0 * pa_xy * fx * pb_yz

                     + 0.5 * pa_x * fx * pb_yyz

                     + pa_xyy * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_yyz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pb_yyz,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (6.0 * pa_x * fx * fx * fz * pb_z

                     - 0.5 * pa_x * fx * fz * fgb * pb_z

                     - 0.5 * pa_x * fz * fga * fx * pb_z

                     - pa_xyy * fz * fgb * pb_z

                     + 5.0 * pa_xyy * fz * fx * pb_z

                     + 20.0 * pa_xy * fz * fx * pb_yz

                     - pa_x * fz * fga * pb_yyz

                     + 5.0 * pa_x * fz * fx * pb_yyz

                     + 12.0 * pa_xyy * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_yzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pb_y,
                                   double pb_yzz,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xy * fx * fx

                     + 0.25 * pa_x * fx * fx * pb_y

                     + 0.5 * pa_xyy * pb_y * fx

                     + pa_xy * fx * pb_zz

                     + 0.5 * pa_x * fx * pb_yzz

                     + pa_xyy * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_yzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pb_y,
                                   double pb_yzz,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_xy * fx * fz * fgb

                     + 4.0 * pa_xy * fz * fx * fx

                     - 0.5 * pa_x * fx * pb_y * fz * fgb

                     - 0.5 * pa_x * fz * fga * pb_y * fx

                     - pa_xyy * pb_y * fz * fgb

                     + 2.0 * pa_x * fz * fx * fx * pb_y

                     + 5.0 * pa_xyy * fz * pb_y * fx

                     + 10.0 * pa_xy * fz * fx * pb_zz

                     - pa_x * fz * fga * pb_yzz

                     + 5.0 * pa_x * fz * fx * pb_yzz

                     + 12.0 * pa_xyy * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_zzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xyy,
                                   double pb_z,
                                   double pb_zzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * pb_z

                     + 1.5 * pa_xyy * pb_z * fx

                     + 0.5 * pa_x * fx * pb_zzz

                     + pa_xyy * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_zzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xyy,
                                   double pb_z,
                                   double pb_zzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * pb_z * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_z * fx

                     - 3.0 * pa_xyy * pb_z * fz * fgb

                     + 6.0 * pa_x * fz * fx * fx * pb_z

                     + 15.0 * pa_xyy * fz * pb_z * fx

                     - pa_x * fz * fga * pb_zzz

                     + 5.0 * pa_x * fz * fx * pb_zzz

                     + 12.0 * pa_xyy * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xxx_s_0(double fx,
                                   double pa_xyz,
                                   double pa_yz,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_yz

                     + 1.5 * pa_xyz * pb_x * fx

                     + 1.5 * fx * pa_yz * pb_xx

                     + pa_xyz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xxx_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xyz,
                                   double pa_yz,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pa_yz * fz * fgb

                     + 6.0 * fx * fx * pa_yz * fz

                     - 3.0 * pa_xyz * pb_x * fz * fgb

                     + 15.0 * pa_xyz * fz * pb_x * fx

                     + 15.0 * fx * pa_yz * fz * pb_xx

                     + 12.0 * pa_xyz * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xxy_s_0(double fx,
                                   double pa_xyz,
                                   double pa_xz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxy,
                                   double pb_xy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_xz * fx * fx

                     + 0.5 * fx * fx * pa_z * pb_x

                     + 0.5 * pa_xyz * fx * pb_y

                     + 0.5 * pa_xz * fx * pb_xx

                     + fx * pa_yz * pb_xy

                     + pa_xyz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xxy_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xyz,
                                   double pa_xz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxy,
                                   double pb_xy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_xz * fx * fz * fgb

                     + 2.0 * pa_xz * fx * fx * fz

                     + 4.0 * fx * fx * fz * pa_z * pb_x

                     - pa_xyz * fz * fgb * pb_y

                     + 5.0 * pa_xyz * fz * fx * pb_y

                     + 5.0 * pa_xz * fx * fz * pb_xx

                     + 10.0 * fx * pa_yz * fz * pb_xy

                     + 12.0 * pa_xyz * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xxz_s_0(double fx,
                                   double pa_xy,
                                   double pa_xyz,
                                   double pa_y,
                                   double pa_yz,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxz,
                                   double pb_xz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_xy * fx * fx

                     + 0.5 * fx * fx * pa_y * pb_x

                     + 0.5 * pa_xyz * fx * pb_z

                     + 0.5 * pa_xy * fx * pb_xx

                     + fx * pa_yz * pb_xz

                     + pa_xyz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xxz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xy,
                                   double pa_xyz,
                                   double pa_y,
                                   double pa_yz,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxz,
                                   double pb_xz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_xy * fx * fz * fgb

                     + 2.0 * pa_xy * fz * fx * fx

                     + 4.0 * fx * fx * pa_y * fz * pb_x

                     - pa_xyz * fz * fgb * pb_z

                     + 5.0 * pa_xyz * fz * fx * pb_z

                     + 5.0 * pa_xy * fz * fx * pb_xx

                     + 10.0 * fx * pa_yz * fz * pb_xz

                     + 12.0 * pa_xyz * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xyy_s_0(double fx,
                                   double pa_xyz,
                                   double pa_xz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyy,
                                   double pb_y,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_yz

                     + 0.5 * fx * fx * pa_z * pb_y

                     + 0.5 * pa_xyz * pb_x * fx

                     + pa_xz * fx * pb_xy

                     + 0.5 * fx * pa_yz * pb_yy

                     + pa_xyz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xyy_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xyz,
                                   double pa_xz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyy,
                                   double pb_y,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_yz * fz * fgb

                     + 2.0 * fx * fx * pa_yz * fz

                     + 4.0 * fx * fx * fz * pa_z * pb_y

                     - pa_xyz * pb_x * fz * fgb

                     + 5.0 * pa_xyz * fz * pb_x * fx

                     + 10.0 * pa_xz * fx * fz * pb_xy

                     + 5.0 * fx * pa_yz * fz * pb_yy

                     + 12.0 * pa_xyz * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xyz_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyz,
                                   double pa_xz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyz,
                                   double pb_xz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx

                     + 0.25 * pa_x * fx * fx * pb_x

                     + 0.25 * fx * fx * pa_y * pb_y

                     + 0.25 * fx * fx * pa_z * pb_z

                     + 0.5 * pa_xy * fx * pb_xy

                     + 0.5 * pa_xz * fx * pb_xz

                     + 0.5 * fx * pa_yz * pb_yz

                     + pa_xyz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xyz_r_0(double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyz,
                                   double pa_xz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyz,
                                   double pb_xz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (0.75 * fx * fx * fx * fz

                     + 2.0 * pa_x * fx * fx * fz * pb_x

                     + 2.0 * fx * fx * pa_y * fz * pb_y

                     + 2.0 * fx * fx * fz * pa_z * pb_z

                     + 5.0 * pa_xy * fz * fx * pb_xy

                     + 5.0 * pa_xz * fx * fz * pb_xz

                     + 5.0 * fx * pa_yz * fz * pb_yz

                     + 12.0 * pa_xyz * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xzz_s_0(double fx,
                                   double pa_xy,
                                   double pa_xyz,
                                   double pa_y,
                                   double pa_yz,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_xzz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_yz

                     + 0.5 * fx * fx * pa_y * pb_z

                     + 0.5 * pa_xyz * pb_x * fx

                     + pa_xy * fx * pb_xz

                     + 0.5 * fx * pa_yz * pb_zz

                     + pa_xyz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xy,
                                   double pa_xyz,
                                   double pa_y,
                                   double pa_yz,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_xzz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_yz * fz * fgb

                     + 2.0 * fx * fx * pa_yz * fz

                     + 4.0 * fx * fx * pa_y * fz * pb_z

                     - pa_xyz * pb_x * fz * fgb

                     + 5.0 * pa_xyz * fz * pb_x * fx

                     + 10.0 * pa_xy * fz * fx * pb_xz

                     + 5.0 * fx * pa_yz * fz * pb_zz

                     + 12.0 * pa_xyz * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_yyy_s_0(double fx,
                                   double pa_xyz,
                                   double pa_xz,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xz * fx * fx

                     + 1.5 * pa_xyz * pb_y * fx

                     + 1.5 * pa_xz * fx * pb_yy

                     + pa_xyz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_yyy_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xyz,
                                   double pa_xz,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xz * fx * fz * fgb

                     + 6.0 * pa_xz * fx * fx * fz

                     - 3.0 * pa_xyz * pb_y * fz * fgb

                     + 15.0 * pa_xyz * fz * pb_y * fx

                     + 15.0 * pa_xz * fx * fz * pb_yy

                     + 12.0 * pa_xyz * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_yyz_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyz,
                                   double pa_xz,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyz,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_xy * fx * fx

                     + 0.5 * pa_x * fx * fx * pb_y

                     + 0.5 * pa_xyz * fx * pb_z

                     + 0.5 * pa_xy * fx * pb_yy

                     + pa_xz * fx * pb_yz

                     + pa_xyz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_yyz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyz,
                                   double pa_xz,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyz,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_xy * fx * fz * fgb

                     + 2.0 * pa_xy * fz * fx * fx

                     + 4.0 * pa_x * fx * fx * fz * pb_y

                     - pa_xyz * fz * fgb * pb_z

                     + 5.0 * pa_xyz * fz * fx * pb_z

                     + 5.0 * pa_xy * fz * fx * pb_yy

                     + 10.0 * pa_xz * fx * fz * pb_yz

                     + 12.0 * pa_xyz * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_yzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyz,
                                   double pa_xz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_yzz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_xz * fx * fx

                     + 0.5 * pa_x * fx * fx * pb_z

                     + 0.5 * pa_xyz * pb_y * fx

                     + pa_xy * fx * pb_yz

                     + 0.5 * pa_xz * fx * pb_zz

                     + pa_xyz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_yzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyz,
                                   double pa_xz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_yzz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_xz * fx * fz * fgb

                     + 2.0 * pa_xz * fx * fx * fz

                     + 4.0 * pa_x * fx * fx * fz * pb_z

                     - pa_xyz * pb_y * fz * fgb

                     + 5.0 * pa_xyz * fz * pb_y * fx

                     + 10.0 * pa_xy * fz * fx * pb_yz

                     + 5.0 * pa_xz * fx * fz * pb_zz

                     + 12.0 * pa_xyz * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_zzz_s_0(double fx,
                                   double pa_xy,
                                   double pa_xyz,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xy * fx * fx

                     + 1.5 * pa_xyz * pb_z * fx

                     + 1.5 * pa_xy * fx * pb_zz

                     + pa_xyz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_zzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xy,
                                   double pa_xyz,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xy * fx * fz * fgb

                     + 6.0 * pa_xy * fz * fx * fx

                     - 3.0 * pa_xyz * pb_z * fz * fgb

                     + 15.0 * pa_xyz * fz * pb_z * fx

                     + 15.0 * pa_xy * fz * fx * pb_zz

                     + 12.0 * pa_xyz * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xxx_s_0(double fx,
                                   double pa_x,
                                   double pa_xzz,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * fx * fx * pa_zz

                     + 0.75 * pa_x * fx * fx * pb_x

                     + 0.75 * fx * fx * pb_xx

                     + 1.5 * pa_xzz * pb_x * fx

                     + 1.5 * fx * pa_zz * pb_xx

                     + 0.5 * pa_x * fx * pb_xxx

                     + pa_xzz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xxx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xzz,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     - 0.75 * fx * fx * fz * fga

                     - 1.5 * fx * pa_zz * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * fx * fx * pa_zz * fz

                     - 1.5 * pa_x * fx * pb_x * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_x * fx

                     - 1.5 * fx * fz * fga * pb_xx

                     - 3.0 * pa_xzz * pb_x * fz * fgb

                     + 6.0 * pa_x * fz * fx * fx * pb_x

                     + 6.0 * fx * fx * fz * pb_xx

                     + 15.0 * pa_xzz * fz * pb_x * fx

                     + 15.0 * fx * pa_zz * fz * pb_xx

                     - pa_x * fz * fga * pb_xxx

                     + 5.0 * pa_x * fz * fx * pb_xxx

                     + 12.0 * pa_xzz * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xxy_s_0(double fx,
                                   double pa_x,
                                   double pa_xzz,
                                   double pa_zz,
                                   double pb_xxy,
                                   double pb_xy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_x * fx * fx * pb_y

                     + 0.5 * fx * fx * pb_xy

                     + 0.5 * pa_xzz * fx * pb_y

                     + fx * pa_zz * pb_xy

                     + 0.5 * pa_x * fx * pb_xxy

                     + pa_xzz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xxy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xzz,
                                   double pa_zz,
                                   double pb_xxy,
                                   double pb_xy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_x * fx * fz * fgb * pb_y

                     - 0.5 * pa_x * fz * fga * fx * pb_y

                     - fx * fz * fga * pb_xy

                     - pa_xzz * fz * fgb * pb_y

                     + 2.0 * pa_x * fz * fx * fx * pb_y

                     + 4.0 * fx * fx * fz * pb_xy

                     + 5.0 * pa_xzz * fz * fx * pb_y

                     + 10.0 * fx * pa_zz * fz * pb_xy

                     - pa_x * fz * fga * pb_xxy

                     + 5.0 * pa_x * fz * fx * pb_xxy

                     + 12.0 * pa_xzz * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xxz_s_0(double fx,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxz,
                                   double pb_xz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xz * fx * fx

                     + fx * fx * pa_z * pb_x

                     + 0.25 * pa_x * fx * fx * pb_z

                     + 0.5 * fx * fx * pb_xz

                     + 0.5 * pa_xzz * fx * pb_z

                     + pa_xz * fx * pb_xx

                     + fx * pa_zz * pb_xz

                     + 0.5 * pa_x * fx * pb_xxz

                     + pa_xzz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xxz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxz,
                                   double pb_xz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_xz * fx * fz * fgb

                     + 4.0 * pa_xz * fz * fx * fx

                     + 8.0 * fx * fx * pa_z * fz * pb_x

                     - 0.5 * pa_x * fx * fz * fgb * pb_z

                     - 0.5 * pa_x * fz * fga * fx * pb_z

                     - fx * fz * fga * pb_xz

                     - pa_xzz * fz * fgb * pb_z

                     + 2.0 * pa_x * fz * fx * fx * pb_z

                     + 4.0 * fx * fx * fz * pb_xz

                     + 5.0 * pa_xzz * fz * fx * pb_z

                     + 10.0 * pa_xz * fz * fx * pb_xx

                     + 10.0 * fx * pa_zz * fz * pb_xz

                     - pa_x * fz * fga * pb_xxz

                     + 5.0 * pa_x * fz * fx * pb_xxz

                     + 12.0 * pa_xzz * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xyy_s_0(double fx,
                                   double pa_x,
                                   double pa_xzz,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xyy,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx

                     + 0.25 * fx * fx * pa_zz

                     + 0.25 * pa_x * fx * fx * pb_x

                     + 0.25 * fx * fx * pb_yy

                     + 0.5 * pa_xzz * pb_x * fx

                     + 0.5 * fx * pa_zz * pb_yy

                     + 0.5 * pa_x * fx * pb_xyy

                     + pa_xzz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xzz,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xyy,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.25 * fx * fx * fz * fgb

                     - 0.25 * fx * fx * fz * fga

                     - 0.5 * fx * pa_zz * fz * fgb

                     + 0.75 * fx * fx * fx * fz

                     + 2.0 * fx * fx * pa_zz * fz

                     - 0.5 * pa_x * fx * pb_x * fz * fgb

                     - 0.5 * pa_x * fz * fga * pb_x * fx

                     - 0.5 * fx * fz * fga * pb_yy

                     - pa_xzz * pb_x * fz * fgb

                     + 2.0 * pa_x * fz * fx * fx * pb_x

                     + 2.0 * fx * fx * fz * pb_yy

                     + 5.0 * pa_xzz * fz * pb_x * fx

                     + 5.0 * fx * pa_zz * fz * pb_yy

                     - pa_x * fz * fga * pb_xyy

                     + 5.0 * pa_x * fz * fx * pb_xyy

                     + 12.0 * pa_xzz * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xyz_s_0(double fx,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_xy,
                                   double pb_xyz,
                                   double pb_y,
                                   double pb_yz,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * fx * fx * pa_z * pb_y

                     + 0.25 * fx * fx * pb_yz

                     + pa_xz * fx * pb_xy

                     + 0.5 * fx * pa_zz * pb_yz

                     + 0.5 * pa_x * fx * pb_xyz

                     + pa_xzz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xyz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_xy,
                                   double pb_xyz,
                                   double pb_y,
                                   double pb_yz,
                                   double r_0_0)
    {
        return r_0_0 * (4.0 * fx * fx * pa_z * fz * pb_y

                     - 0.5 * fx * fz * fga * pb_yz

                     + 2.0 * fx * fx * fz * pb_yz

                     + 10.0 * pa_xz * fz * fx * pb_xy

                     + 5.0 * fx * pa_zz * fz * pb_yz

                     - pa_x * fz * fga * pb_xyz

                     + 5.0 * pa_x * fz * fx * pb_xyz

                     + 12.0 * pa_xzz * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_xzz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_x * fx * fx * pb_x

                     + 0.25 * fx * fx * pa_zz

                     + fx * fx * pa_z * pb_z

                     + 0.25 * fx * fx * pb_zz

                     + 0.5 * pa_xzz * pb_x * fx

                     + 2.0 * pa_xz * fx * pb_xz

                     + 0.5 * fx * pa_zz * pb_zz

                     + 0.5 * pa_x * fx * pb_xzz

                     + pa_xzz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_xzz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (2.25 * fx * fx * fx * fz

                     - 0.25 * fx * fx * fz * fgb

                     - 0.25 * fx * fx * fz * fga

                     - 0.5 * fx * pa_zz * fz * fgb

                     + 6.0 * pa_x * fx * fx * fz * pb_x

                     + 2.0 * fx * fx * pa_zz * fz

                     + 8.0 * fx * fx * pa_z * fz * pb_z

                     - 0.5 * pa_x * fx * pb_x * fz * fgb

                     - 0.5 * pa_x * fz * fga * pb_x * fx

                     - 0.5 * fx * fz * fga * pb_zz

                     - pa_xzz * pb_x * fz * fgb

                     + 2.0 * fx * fx * fz * pb_zz

                     + 5.0 * pa_xzz * fz * pb_x * fx

                     + 20.0 * pa_xz * fz * fx * pb_xz

                     + 5.0 * fx * pa_zz * fz * pb_zz

                     - pa_x * fz * fga * pb_xzz

                     + 5.0 * pa_x * fz * fx * pb_xzz

                     + 12.0 * pa_xzz * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_yyy_s_0(double fx,
                                   double pa_x,
                                   double pa_xzz,
                                   double pb_y,
                                   double pb_yyy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * pb_y

                     + 1.5 * pa_xzz * pb_y * fx

                     + 0.5 * pa_x * fx * pb_yyy

                     + pa_xzz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_yyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xzz,
                                   double pb_y,
                                   double pb_yyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * pb_y * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_y * fx

                     - 3.0 * pa_xzz * pb_y * fz * fgb

                     + 6.0 * pa_x * fz * fx * fx * pb_y

                     + 15.0 * pa_xzz * fz * pb_y * fx

                     - pa_x * fz * fga * pb_yyy

                     + 5.0 * pa_x * fz * fx * pb_yyy

                     + 12.0 * pa_xzz * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_yyz_s_0(double fx,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pb_yy,
                                   double pb_yyz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xz * fx * fx

                     + 0.25 * pa_x * fx * fx * pb_z

                     + 0.5 * pa_xzz * fx * pb_z

                     + pa_xz * fx * pb_yy

                     + 0.5 * pa_x * fx * pb_yyz

                     + pa_xzz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_yyz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pb_yy,
                                   double pb_yyz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_xz * fx * fz * fgb

                     + 4.0 * pa_xz * fz * fx * fx

                     - 0.5 * pa_x * fx * fz * fgb * pb_z

                     - 0.5 * pa_x * fz * fga * fx * pb_z

                     - pa_xzz * fz * fgb * pb_z

                     + 2.0 * pa_x * fz * fx * fx * pb_z

                     + 5.0 * pa_xzz * fz * fx * pb_z

                     + 10.0 * pa_xz * fz * fx * pb_yy

                     - pa_x * fz * fga * pb_yyz

                     + 5.0 * pa_x * fz * fx * pb_yyz

                     + 12.0 * pa_xzz * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_yzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_yzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * pb_y

                     + 0.5 * pa_xzz * pb_y * fx

                     + 2.0 * pa_xz * fx * pb_yz

                     + 0.5 * pa_x * fx * pb_yzz

                     + pa_xzz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_yzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_yzz,
                                   double r_0_0)
    {
        return r_0_0 * (6.0 * pa_x * fx * fx * fz * pb_y

                     - 0.5 * pa_x * fx * pb_y * fz * fgb

                     - 0.5 * pa_x * fz * fga * pb_y * fx

                     - pa_xzz * pb_y * fz * fgb

                     + 5.0 * pa_xzz * fz * pb_y * fx

                     + 20.0 * pa_xz * fz * fx * pb_yz

                     - pa_x * fz * fga * pb_yzz

                     + 5.0 * pa_x * fz * fx * pb_yzz

                     + 12.0 * pa_xzz * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_zzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xz * fx * fx

                     + 2.25 * pa_x * fx * fx * pb_z

                     + 1.5 * pa_xzz * pb_z * fx

                     + 3.0 * pa_xz * fx * pb_zz

                     + 0.5 * pa_x * fx * pb_zzz

                     + pa_xzz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_zzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xz * fx * fz * fgb

                     + 12.0 * pa_xz * fz * fx * fx

                     + 18.0 * pa_x * fx * fx * fz * pb_z

                     - 1.5 * pa_x * fx * pb_z * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_z * fx

                     - 3.0 * pa_xzz * pb_z * fz * fgb

                     + 15.0 * pa_xzz * fz * pb_z * fx

                     + 30.0 * pa_xz * fz * fx * pb_zz

                     - pa_x * fz * fga * pb_zzz

                     + 5.0 * pa_x * fz * fx * pb_zzz

                     + 12.0 * pa_xzz * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xxx_s_0(double fx,
                                   double pa_y,
                                   double pa_yyy,
                                   double pb_x,
                                   double pb_xxx,
                                   double s_0_0)
    {
        return s_0_0 * (2.25 * pa_y * fx * fx * pb_x

                     + 1.5 * pa_yyy * pb_x * fx

                     + 1.5 * pa_y * fx * pb_xxx

                     + pa_yyy * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xxx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yyy,
                                   double pb_x,
                                   double pb_xxx,
                                   double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_y * fx * pb_x * fz * fgb

                     - 4.5 * pa_y * fz * fga * pb_x * fx

                     - 3.0 * pa_yyy * pb_x * fz * fgb

                     + 18.0 * pa_y * fz * fx * fx * pb_x

                     + 15.0 * pa_yyy * fz * pb_x * fx

                     - 3.0 * pa_y * fz * fga * pb_xxx

                     + 15.0 * pa_y * fz * fx * pb_xxx

                     + 12.0 * pa_yyy * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xxy_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pb_xx,
                                   double pb_xxy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_yy * fx * fx

                     + 0.75 * pa_y * fx * fx * pb_y

                     + 0.75 * fx * fx * pb_xx

                     + 0.5 * pa_yyy * fx * pb_y

                     + 1.5 * pa_yy * fx * pb_xx

                     + 1.5 * pa_y * fx * pb_xxy

                     + pa_yyy * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xxy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pb_xx,
                                   double pb_xxy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     - 0.75 * fx * fx * fz * fga

                     - 1.5 * pa_yy * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * pa_yy * fz * fx * fx

                     - 1.5 * pa_y * fx * fz * fgb * pb_y

                     - 1.5 * pa_y * fz * fga * fx * pb_y

                     - 1.5 * fx * fz * fga * pb_xx

                     - pa_yyy * fz * fgb * pb_y

                     + 6.0 * pa_y * fz * fx * fx * pb_y

                     + 6.0 * fx * fx * fz * pb_xx

                     + 5.0 * pa_yyy * fz * fx * pb_y

                     + 15.0 * pa_yy * fz * fx * pb_xx

                     - 3.0 * pa_y * fz * fga * pb_xxy

                     + 15.0 * pa_y * fz * fx * pb_xxy

                     + 12.0 * pa_yyy * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xxz_s_0(double fx,
                                   double pa_y,
                                   double pa_yyy,
                                   double pb_xxz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx * pb_z

                     + 0.5 * pa_yyy * fx * pb_z

                     + 1.5 * pa_y * fx * pb_xxz

                     + pa_yyy * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xxz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yyy,
                                   double pb_xxz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_y * fx * fz * fgb * pb_z

                     - 1.5 * pa_y * fz * fga * fx * pb_z

                     - pa_yyy * fz * fgb * pb_z

                     + 6.0 * pa_y * fz * fx * fx * pb_z

                     + 5.0 * pa_yyy * fz * fx * pb_z

                     - 3.0 * pa_y * fz * fga * pb_xxz

                     + 15.0 * pa_y * fz * fx * pb_xxz

                     + 12.0 * pa_yyy * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xyy_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyy,
                                   double s_0_0)
    {
        return s_0_0 * (2.25 * pa_y * fx * fx * pb_x

                     + 1.5 * fx * fx * pb_xy

                     + 0.5 * pa_yyy * pb_x * fx

                     + 3.0 * pa_yy * fx * pb_xy

                     + 1.5 * pa_y * fx * pb_xyy

                     + pa_yyy * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyy,
                                   double r_0_0)
    {
        return r_0_0 * (18.0 * pa_y * fx * fx * fz * pb_x

                     - 1.5 * pa_y * fx * pb_x * fz * fgb

                     - 1.5 * pa_y * fz * fga * pb_x * fx

                     - 3.0 * fx * fz * fga * pb_xy

                     - pa_yyy * pb_x * fz * fgb

                     + 12.0 * fx * fx * fz * pb_xy

                     + 5.0 * pa_yyy * fz * pb_x * fx

                     + 30.0 * pa_yy * fz * fx * pb_xy

                     - 3.0 * pa_y * fz * fga * pb_xyy

                     + 15.0 * pa_y * fz * fx * pb_xyy

                     + 12.0 * pa_yyy * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xyz_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pb_xyz,
                                   double pb_xz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_xz

                     + 1.5 * pa_yy * fx * pb_xz

                     + 1.5 * pa_y * fx * pb_xyz

                     + pa_yyy * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xyz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pb_xyz,
                                   double pb_xz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fz * fga * pb_xz

                     + 6.0 * fx * fx * fz * pb_xz

                     + 15.0 * pa_yy * fz * fx * pb_xz

                     - 3.0 * pa_y * fz * fga * pb_xyz

                     + 15.0 * pa_y * fz * fx * pb_xyz

                     + 12.0 * pa_yyy * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xzz_s_0(double fx,
                                   double pa_y,
                                   double pa_yyy,
                                   double pb_x,
                                   double pb_xzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx * pb_x

                     + 0.5 * pa_yyy * pb_x * fx

                     + 1.5 * pa_y * fx * pb_xzz

                     + pa_yyy * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yyy,
                                   double pb_x,
                                   double pb_xzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_y * fx * pb_x * fz * fgb

                     - 1.5 * pa_y * fz * fga * pb_x * fx

                     - pa_yyy * pb_x * fz * fgb

                     + 6.0 * pa_y * fz * fx * fx * pb_x

                     + 5.0 * pa_yyy * fz * pb_x * fx

                     - 3.0 * pa_y * fz * fga * pb_xzz

                     + 15.0 * pa_y * fz * fx * pb_xzz

                     + 12.0 * pa_yyy * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_yyy_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx

                     + 2.25 * pa_yy * fx * fx

                     + 6.75 * pa_y * fx * fx * pb_y

                     + 2.25 * fx * fx * pb_yy

                     + 1.5 * pa_yyy * pb_y * fx

                     + 4.5 * pa_yy * fx * pb_yy

                     + 1.5 * pa_y * fx * pb_yyy

                     + pa_yyy * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_yyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double r_0_0)
    {
        return r_0_0 * (11.25 * fx * fx * fx * fz

                     - 2.25 * fx * fx * fz * fgb

                     - 2.25 * fx * fx * fz * fga

                     - 4.5 * pa_yy * fx * fz * fgb

                     + 18.0 * pa_yy * fz * fx * fx

                     + 54.0 * pa_y * fx * fx * fz * pb_y

                     - 4.5 * pa_y * fx * pb_y * fz * fgb

                     - 4.5 * pa_y * fz * fga * pb_y * fx

                     - 4.5 * fx * fz * fga * pb_yy

                     - 3.0 * pa_yyy * pb_y * fz * fgb

                     + 18.0 * fx * fx * fz * pb_yy

                     + 15.0 * pa_yyy * fz * pb_y * fx

                     + 45.0 * pa_yy * fz * fx * pb_yy

                     - 3.0 * pa_y * fz * fga * pb_yyy

                     + 15.0 * pa_y * fz * fx * pb_yyy

                     + 12.0 * pa_yyy * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_yyz_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pb_yyz,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (2.25 * pa_y * fx * fx * pb_z

                     + 1.5 * fx * fx * pb_yz

                     + 0.5 * pa_yyy * fx * pb_z

                     + 3.0 * pa_yy * fx * pb_yz

                     + 1.5 * pa_y * fx * pb_yyz

                     + pa_yyy * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_yyz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pb_yyz,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (18.0 * pa_y * fx * fx * fz * pb_z

                     - 1.5 * pa_y * fx * fz * fgb * pb_z

                     - 1.5 * pa_y * fz * fga * fx * pb_z

                     - 3.0 * fx * fz * fga * pb_yz

                     - pa_yyy * fz * fgb * pb_z

                     + 12.0 * fx * fx * fz * pb_yz

                     + 5.0 * pa_yyy * fz * fx * pb_z

                     + 30.0 * pa_yy * fz * fx * pb_yz

                     - 3.0 * pa_y * fz * fga * pb_yyz

                     + 15.0 * pa_y * fz * fx * pb_yyz

                     + 12.0 * pa_yyy * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_yzz_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pb_y,
                                   double pb_yzz,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_yy * fx * fx

                     + 0.75 * pa_y * fx * fx * pb_y

                     + 0.75 * fx * fx * pb_zz

                     + 0.5 * pa_yyy * pb_y * fx

                     + 1.5 * pa_yy * fx * pb_zz

                     + 1.5 * pa_y * fx * pb_yzz

                     + pa_yyy * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_yzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pb_y,
                                   double pb_yzz,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     - 0.75 * fx * fx * fz * fga

                     - 1.5 * pa_yy * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * pa_yy * fz * fx * fx

                     - 1.5 * pa_y * fx * pb_y * fz * fgb

                     - 1.5 * pa_y * fz * fga * pb_y * fx

                     - 1.5 * fx * fz * fga * pb_zz

                     - pa_yyy * pb_y * fz * fgb

                     + 6.0 * pa_y * fz * fx * fx * pb_y

                     + 6.0 * fx * fx * fz * pb_zz

                     + 5.0 * pa_yyy * fz * pb_y * fx

                     + 15.0 * pa_yy * fz * fx * pb_zz

                     - 3.0 * pa_y * fz * fga * pb_yzz

                     + 15.0 * pa_y * fz * fx * pb_yzz

                     + 12.0 * pa_yyy * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_zzz_s_0(double fx,
                                   double pa_y,
                                   double pa_yyy,
                                   double pb_z,
                                   double pb_zzz,
                                   double s_0_0)
    {
        return s_0_0 * (2.25 * pa_y * fx * fx * pb_z

                     + 1.5 * pa_yyy * pb_z * fx

                     + 1.5 * pa_y * fx * pb_zzz

                     + pa_yyy * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_zzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yyy,
                                   double pb_z,
                                   double pb_zzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_y * fx * pb_z * fz * fgb

                     - 4.5 * pa_y * fz * fga * pb_z * fx

                     - 3.0 * pa_yyy * pb_z * fz * fgb

                     + 18.0 * pa_y * fz * fx * fx * pb_z

                     + 15.0 * pa_yyy * fz * pb_z * fx

                     - 3.0 * pa_y * fz * fga * pb_zzz

                     + 15.0 * pa_y * fz * fx * pb_zzz

                     + 12.0 * pa_yyy * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xxx_s_0(double fx,
                                   double pa_yyz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xxx,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_z * pb_x

                     + 1.5 * pa_yyz * pb_x * fx

                     + 0.5 * fx * pa_z * pb_xxx

                     + pa_yyz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xxx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yyz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xxx,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pa_z * pb_x * fz * fgb

                     - 1.5 * fz * fga * pa_z * pb_x * fx

                     - 3.0 * pa_yyz * pb_x * fz * fgb

                     + 6.0 * fx * fx * fz * pa_z * pb_x

                     + 15.0 * pa_yyz * fz * pb_x * fx

                     - fz * fga * pa_z * pb_xxx

                     + 5.0 * fx * fz * pa_z * pb_xxx

                     + 12.0 * pa_yyz * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xxy_s_0(double fx,
                                   double pa_yyz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_xx,
                                   double pb_xxy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_yz * fx * fx

                     + 0.25 * fx * fx * pa_z * pb_y

                     + 0.5 * pa_yyz * fx * pb_y

                     + pa_yz * fx * pb_xx

                     + 0.5 * fx * pa_z * pb_xxy

                     + pa_yyz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xxy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yyz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_xx,
                                   double pb_xxy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_yz * fx * fz * fgb

                     + 4.0 * pa_yz * fx * fx * fz

                     - 0.5 * fx * pa_z * fz * fgb * pb_y

                     - 0.5 * fz * fga * pa_z * fx * pb_y

                     - pa_yyz * fz * fgb * pb_y

                     + 2.0 * fx * fx * fz * pa_z * pb_y

                     + 5.0 * pa_yyz * fz * fx * pb_y

                     + 10.0 * pa_yz * fx * fz * pb_xx

                     - fz * fga * pa_z * pb_xxy

                     + 5.0 * fx * fz * pa_z * pb_xxy

                     + 12.0 * pa_yyz * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xxz_s_0(double fx,
                                   double pa_yy,
                                   double pa_yyz,
                                   double pa_z,
                                   double pb_xx,
                                   double pb_xxz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx

                     + 0.25 * pa_yy * fx * fx

                     + 0.25 * fx * fx * pa_z * pb_z

                     + 0.25 * fx * fx * pb_xx

                     + 0.5 * pa_yyz * fx * pb_z

                     + 0.5 * pa_yy * fx * pb_xx

                     + 0.5 * fx * pa_z * pb_xxz

                     + pa_yyz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xxz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yy,
                                   double pa_yyz,
                                   double pa_z,
                                   double pb_xx,
                                   double pb_xxz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.25 * fx * fx * fz * fgb

                     - 0.25 * fz * fga * fx * fx

                     - 0.5 * pa_yy * fx * fz * fgb

                     + 0.75 * fx * fx * fx * fz

                     + 2.0 * pa_yy * fz * fx * fx

                     - 0.5 * fx * pa_z * fz * fgb * pb_z

                     - 0.5 * fz * fga * pa_z * fx * pb_z

                     - 0.5 * fz * fga * fx * pb_xx

                     - pa_yyz * fz * fgb * pb_z

                     + 2.0 * fx * fx * fz * pa_z * pb_z

                     + 2.0 * fx * fx * fz * pb_xx

                     + 5.0 * pa_yyz * fz * fx * pb_z

                     + 5.0 * pa_yy * fz * fx * pb_xx

                     - fz * fga * pa_z * pb_xxz

                     + 5.0 * fx * fz * pa_z * pb_xxz

                     + 12.0 * pa_yyz * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xyy_s_0(double fx,
                                   double pa_yyz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_z * pb_x

                     + 0.5 * pa_yyz * pb_x * fx

                     + 2.0 * pa_yz * fx * pb_xy

                     + 0.5 * fx * pa_z * pb_xyy

                     + pa_yyz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yyz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyy,
                                   double r_0_0)
    {
        return r_0_0 * (6.0 * fx * fx * fz * pa_z * pb_x

                     - 0.5 * fx * pa_z * pb_x * fz * fgb

                     - 0.5 * fz * fga * pa_z * pb_x * fx

                     - pa_yyz * pb_x * fz * fgb

                     + 5.0 * pa_yyz * fz * pb_x * fx

                     + 20.0 * pa_yz * fx * fz * pb_xy

                     - fz * fga * pa_z * pb_xyy

                     + 5.0 * fx * fz * pa_z * pb_xyy

                     + 12.0 * pa_yyz * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xyz_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyz,
                                   double pb_xz,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_y * fx * fx * pb_x

                     + 0.25 * fx * fx * pb_xy

                     + 0.5 * pa_yy * fx * pb_xy

                     + pa_yz * fx * pb_xz

                     + 0.5 * fx * pa_z * pb_xyz

                     + pa_yyz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xyz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyz,
                                   double pb_xz,
                                   double r_0_0)
    {
        return r_0_0 * (4.0 * pa_y * fx * fx * fz * pb_x

                     - 0.5 * fz * fga * fx * pb_xy

                     + 2.0 * fx * fx * fz * pb_xy

                     + 5.0 * pa_yy * fz * fx * pb_xy

                     + 10.0 * pa_yz * fx * fz * pb_xz

                     - fz * fga * pa_z * pb_xyz

                     + 5.0 * fx * fz * pa_z * pb_xyz

                     + 12.0 * pa_yyz * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xzz_s_0(double fx,
                                   double pa_yy,
                                   double pa_yyz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_xzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_z * pb_x

                     + 0.5 * fx * fx * pb_xz

                     + 0.5 * pa_yyz * pb_x * fx

                     + pa_yy * fx * pb_xz

                     + 0.5 * fx * pa_z * pb_xzz

                     + pa_yyz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yy,
                                   double pa_yyz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_xzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_z * pb_x * fz * fgb

                     - 0.5 * fz * fga * pa_z * pb_x * fx

                     - fz * fga * fx * pb_xz

                     - pa_yyz * pb_x * fz * fgb

                     + 2.0 * fx * fx * fz * pa_z * pb_x

                     + 4.0 * fx * fx * fz * pb_xz

                     + 5.0 * pa_yyz * fz * pb_x * fx

                     + 10.0 * pa_yy * fz * fx * pb_xz

                     - fz * fga * pa_z * pb_xzz

                     + 5.0 * fx * fz * pa_z * pb_xzz

                     + 12.0 * pa_yyz * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_yyy_s_0(double fx,
                                   double pa_yyz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double s_0_0)
    {
        return s_0_0 * (1.5 * pa_yz * fx * fx

                     + 2.25 * fx * fx * pa_z * pb_y

                     + 1.5 * pa_yyz * pb_y * fx

                     + 3.0 * pa_yz * fx * pb_yy

                     + 0.5 * fx * pa_z * pb_yyy

                     + pa_yyz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_yyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yyz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_yz * fx * fz * fgb

                     + 12.0 * pa_yz * fx * fx * fz

                     + 18.0 * fx * fx * fz * pa_z * pb_y

                     - 1.5 * fx * pa_z * pb_y * fz * fgb

                     - 1.5 * fz * fga * pa_z * pb_y * fx

                     - 3.0 * pa_yyz * pb_y * fz * fgb

                     + 15.0 * pa_yyz * fz * pb_y * fx

                     + 30.0 * pa_yz * fx * fz * pb_yy

                     - fz * fga * pa_z * pb_yyy

                     + 5.0 * fx * fz * pa_z * pb_yyy

                     + 12.0 * pa_yyz * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_yyz_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyz,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.25 * pa_yy * fx * fx

                     + pa_y * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_z * pb_z

                     + 0.25 * fx * fx * pb_yy

                     + 0.5 * pa_yyz * fx * pb_z

                     + 0.5 * pa_yy * fx * pb_yy

                     + 2.0 * pa_yz * fx * pb_yz

                     + 0.5 * fx * pa_z * pb_yyz

                     + pa_yyz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_yyz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyz,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (2.25 * fx * fx * fx * fz

                     - 0.25 * fx * fx * fz * fgb

                     - 0.25 * fz * fga * fx * fx

                     - 0.5 * pa_yy * fx * fz * fgb

                     + 2.0 * pa_yy * fz * fx * fx

                     + 8.0 * pa_y * fx * fx * fz * pb_y

                     + 6.0 * fx * fx * fz * pa_z * pb_z

                     - 0.5 * fx * pa_z * fz * fgb * pb_z

                     - 0.5 * fz * fga * pa_z * fx * pb_z

                     - 0.5 * fz * fga * fx * pb_yy

                     - pa_yyz * fz * fgb * pb_z

                     + 2.0 * fx * fx * fz * pb_yy

                     + 5.0 * pa_yyz * fz * fx * pb_z

                     + 5.0 * pa_yy * fz * fx * pb_yy

                     + 20.0 * pa_yz * fx * fz * pb_yz

                     - fz * fga * pa_z * pb_yyz

                     + 5.0 * fx * fz * pa_z * pb_yyz

                     + 12.0 * pa_yyz * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_yzz_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_yzz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_yz * fx * fx

                     + pa_y * fx * fx * pb_z

                     + 0.25 * fx * fx * pa_z * pb_y

                     + 0.5 * fx * fx * pb_yz

                     + 0.5 * pa_yyz * pb_y * fx

                     + pa_yy * fx * pb_yz

                     + pa_yz * fx * pb_zz

                     + 0.5 * fx * pa_z * pb_yzz

                     + pa_yyz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_yzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_yzz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_yz * fx * fz * fgb

                     + 4.0 * pa_yz * fx * fx * fz

                     + 8.0 * pa_y * fx * fx * fz * pb_z

                     - 0.5 * fx * pa_z * pb_y * fz * fgb

                     - 0.5 * fz * fga * pa_z * pb_y * fx

                     - fz * fga * fx * pb_yz

                     - pa_yyz * pb_y * fz * fgb

                     + 2.0 * fx * fx * fz * pa_z * pb_y

                     + 4.0 * fx * fx * fz * pb_yz

                     + 5.0 * pa_yyz * fz * pb_y * fx

                     + 10.0 * pa_yy * fz * fx * pb_yz

                     + 10.0 * pa_yz * fx * fz * pb_zz

                     - fz * fga * pa_z * pb_yzz

                     + 5.0 * fx * fz * pa_z * pb_yzz

                     + 12.0 * pa_yyz * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_zzz_s_0(double fx,
                                   double pa_yy,
                                   double pa_yyz,
                                   double pa_z,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_yy * fx * fx

                     + 0.75 * fx * fx * pa_z * pb_z

                     + 0.75 * fx * fx * pb_zz

                     + 1.5 * pa_yyz * pb_z * fx

                     + 1.5 * pa_yy * fx * pb_zz

                     + 0.5 * fx * pa_z * pb_zzz

                     + pa_yyz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_zzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yy,
                                   double pa_yyz,
                                   double pa_z,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     - 0.75 * fz * fga * fx * fx

                     - 1.5 * pa_yy * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * pa_yy * fz * fx * fx

                     - 1.5 * fx * pa_z * pb_z * fz * fgb

                     - 1.5 * fz * fga * pa_z * pb_z * fx

                     - 1.5 * fz * fga * fx * pb_zz

                     - 3.0 * pa_yyz * pb_z * fz * fgb

                     + 6.0 * fx * fx * fz * pa_z * pb_z

                     + 6.0 * fx * fx * fz * pb_zz

                     + 15.0 * pa_yyz * fz * pb_z * fx

                     + 15.0 * pa_yy * fz * fx * pb_zz

                     - fz * fga * pa_z * pb_zzz

                     + 5.0 * fx * fz * pa_z * pb_zzz

                     + 12.0 * pa_yyz * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xxx_s_0(double fx,
                                   double pa_y,
                                   double pa_yzz,
                                   double pb_x,
                                   double pb_xxx,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx * pb_x

                     + 1.5 * pa_yzz * pb_x * fx

                     + 0.5 * pa_y * fx * pb_xxx

                     + pa_yzz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xxx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yzz,
                                   double pb_x,
                                   double pb_xxx,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_y * fx * pb_x * fz * fgb

                     - 1.5 * pa_y * fz * fga * pb_x * fx

                     - 3.0 * pa_yzz * pb_x * fz * fgb

                     + 6.0 * pa_y * fz * fx * fx * pb_x

                     + 15.0 * pa_yzz * fz * pb_x * fx

                     - pa_y * fz * fga * pb_xxx

                     + 5.0 * pa_y * fz * fx * pb_xxx

                     + 12.0 * pa_yzz * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xxy_s_0(double fx,
                                   double pa_y,
                                   double pa_yzz,
                                   double pa_zz,
                                   double pb_xx,
                                   double pb_xxy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx

                     + 0.25 * fx * fx * pa_zz

                     + 0.25 * pa_y * fx * fx * pb_y

                     + 0.25 * fx * fx * pb_xx

                     + 0.5 * pa_yzz * fx * pb_y

                     + 0.5 * fx * pa_zz * pb_xx

                     + 0.5 * pa_y * fx * pb_xxy

                     + pa_yzz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xxy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yzz,
                                   double pa_zz,
                                   double pb_xx,
                                   double pb_xxy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.25 * fx * fx * fz * fgb

                     - 0.25 * fx * fx * fz * fga

                     - 0.5 * fx * pa_zz * fz * fgb

                     + 0.75 * fx * fx * fx * fz

                     + 2.0 * fx * fx * pa_zz * fz

                     - 0.5 * pa_y * fx * fz * fgb * pb_y

                     - 0.5 * pa_y * fz * fga * fx * pb_y

                     - 0.5 * fx * fz * fga * pb_xx

                     - pa_yzz * fz * fgb * pb_y

                     + 2.0 * pa_y * fz * fx * fx * pb_y

                     + 2.0 * fx * fx * fz * pb_xx

                     + 5.0 * pa_yzz * fz * fx * pb_y

                     + 5.0 * fx * pa_zz * fz * pb_xx

                     - pa_y * fz * fga * pb_xxy

                     + 5.0 * pa_y * fz * fx * pb_xxy

                     + 12.0 * pa_yzz * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xxz_s_0(double fx,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pb_xx,
                                   double pb_xxz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_yz * fx * fx

                     + 0.25 * pa_y * fx * fx * pb_z

                     + 0.5 * pa_yzz * fx * pb_z

                     + pa_yz * fx * pb_xx

                     + 0.5 * pa_y * fx * pb_xxz

                     + pa_yzz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xxz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pb_xx,
                                   double pb_xxz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_yz * fx * fz * fgb

                     + 4.0 * pa_yz * fz * fx * fx

                     - 0.5 * pa_y * fx * fz * fgb * pb_z

                     - 0.5 * pa_y * fz * fga * fx * pb_z

                     - pa_yzz * fz * fgb * pb_z

                     + 2.0 * pa_y * fz * fx * fx * pb_z

                     + 5.0 * pa_yzz * fz * fx * pb_z

                     + 10.0 * pa_yz * fz * fx * pb_xx

                     - pa_y * fz * fga * pb_xxz

                     + 5.0 * pa_y * fz * fx * pb_xxz

                     + 12.0 * pa_yzz * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xyy_s_0(double fx,
                                   double pa_y,
                                   double pa_yzz,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyy,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_y * fx * fx * pb_x

                     + 0.5 * fx * fx * pb_xy

                     + 0.5 * pa_yzz * pb_x * fx

                     + fx * pa_zz * pb_xy

                     + 0.5 * pa_y * fx * pb_xyy

                     + pa_yzz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yzz,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_y * fx * pb_x * fz * fgb

                     - 0.5 * pa_y * fz * fga * pb_x * fx

                     - fx * fz * fga * pb_xy

                     - pa_yzz * pb_x * fz * fgb

                     + 2.0 * pa_y * fz * fx * fx * pb_x

                     + 4.0 * fx * fx * fz * pb_xy

                     + 5.0 * pa_yzz * fz * pb_x * fx

                     + 10.0 * fx * pa_zz * fz * pb_xy

                     - pa_y * fz * fga * pb_xyy

                     + 5.0 * pa_y * fz * fx * pb_xyy

                     + 12.0 * pa_yzz * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xyz_s_0(double fx,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyz,
                                   double pb_xz,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * fx * fx * pa_z * pb_x

                     + 0.25 * fx * fx * pb_xz

                     + pa_yz * fx * pb_xy

                     + 0.5 * fx * pa_zz * pb_xz

                     + 0.5 * pa_y * fx * pb_xyz

                     + pa_yzz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xyz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyz,
                                   double pb_xz,
                                   double r_0_0)
    {
        return r_0_0 * (4.0 * fx * fx * pa_z * fz * pb_x

                     - 0.5 * fx * fz * fga * pb_xz

                     + 2.0 * fx * fx * fz * pb_xz

                     + 10.0 * pa_yz * fz * fx * pb_xy

                     + 5.0 * fx * pa_zz * fz * pb_xz

                     - pa_y * fz * fga * pb_xyz

                     + 5.0 * pa_y * fz * fx * pb_xyz

                     + 12.0 * pa_yzz * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xzz_s_0(double fx,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_xzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx * pb_x

                     + 0.5 * pa_yzz * pb_x * fx

                     + 2.0 * pa_yz * fx * pb_xz

                     + 0.5 * pa_y * fx * pb_xzz

                     + pa_yzz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_xzz,
                                   double r_0_0)
    {
        return r_0_0 * (6.0 * pa_y * fx * fx * fz * pb_x

                     - 0.5 * pa_y * fx * pb_x * fz * fgb

                     - 0.5 * pa_y * fz * fga * pb_x * fx

                     - pa_yzz * pb_x * fz * fgb

                     + 5.0 * pa_yzz * fz * pb_x * fx

                     + 20.0 * pa_yz * fz * fx * pb_xz

                     - pa_y * fz * fga * pb_xzz

                     + 5.0 * pa_y * fz * fx * pb_xzz

                     + 12.0 * pa_yzz * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_yyy_s_0(double fx,
                                   double pa_y,
                                   double pa_yzz,
                                   double pa_zz,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * fx * fx * pa_zz

                     + 0.75 * pa_y * fx * fx * pb_y

                     + 0.75 * fx * fx * pb_yy

                     + 1.5 * pa_yzz * pb_y * fx

                     + 1.5 * fx * pa_zz * pb_yy

                     + 0.5 * pa_y * fx * pb_yyy

                     + pa_yzz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_yyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yzz,
                                   double pa_zz,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     - 0.75 * fx * fx * fz * fga

                     - 1.5 * fx * pa_zz * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * fx * fx * pa_zz * fz

                     - 1.5 * pa_y * fx * pb_y * fz * fgb

                     - 1.5 * pa_y * fz * fga * pb_y * fx

                     - 1.5 * fx * fz * fga * pb_yy

                     - 3.0 * pa_yzz * pb_y * fz * fgb

                     + 6.0 * pa_y * fz * fx * fx * pb_y

                     + 6.0 * fx * fx * fz * pb_yy

                     + 15.0 * pa_yzz * fz * pb_y * fx

                     + 15.0 * fx * pa_zz * fz * pb_yy

                     - pa_y * fz * fga * pb_yyy

                     + 5.0 * pa_y * fz * fx * pb_yyy

                     + 12.0 * pa_yzz * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_yyz_s_0(double fx,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyz,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_yz * fx * fx

                     + fx * fx * pa_z * pb_y

                     + 0.25 * pa_y * fx * fx * pb_z

                     + 0.5 * fx * fx * pb_yz

                     + 0.5 * pa_yzz * fx * pb_z

                     + pa_yz * fx * pb_yy

                     + fx * pa_zz * pb_yz

                     + 0.5 * pa_y * fx * pb_yyz

                     + pa_yzz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_yyz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyz,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_yz * fx * fz * fgb

                     + 4.0 * pa_yz * fz * fx * fx

                     + 8.0 * fx * fx * pa_z * fz * pb_y

                     - 0.5 * pa_y * fx * fz * fgb * pb_z

                     - 0.5 * pa_y * fz * fga * fx * pb_z

                     - fx * fz * fga * pb_yz

                     - pa_yzz * fz * fgb * pb_z

                     + 2.0 * pa_y * fz * fx * fx * pb_z

                     + 4.0 * fx * fx * fz * pb_yz

                     + 5.0 * pa_yzz * fz * fx * pb_z

                     + 10.0 * pa_yz * fz * fx * pb_yy

                     + 10.0 * fx * pa_zz * fz * pb_yz

                     - pa_y * fz * fga * pb_yyz

                     + 5.0 * pa_y * fz * fx * pb_yyz

                     + 12.0 * pa_yzz * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_yzz_s_0(double fx,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_yzz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_y * fx * fx * pb_y

                     + 0.25 * fx * fx * pa_zz

                     + fx * fx * pa_z * pb_z

                     + 0.25 * fx * fx * pb_zz

                     + 0.5 * pa_yzz * pb_y * fx

                     + 2.0 * pa_yz * fx * pb_yz

                     + 0.5 * fx * pa_zz * pb_zz

                     + 0.5 * pa_y * fx * pb_yzz

                     + pa_yzz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_yzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_yzz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (2.25 * fx * fx * fx * fz

                     - 0.25 * fx * fx * fz * fgb

                     - 0.25 * fx * fx * fz * fga

                     - 0.5 * fx * pa_zz * fz * fgb

                     + 6.0 * pa_y * fx * fx * fz * pb_y

                     + 2.0 * fx * fx * pa_zz * fz

                     + 8.0 * fx * fx * pa_z * fz * pb_z

                     - 0.5 * pa_y * fx * pb_y * fz * fgb

                     - 0.5 * pa_y * fz * fga * pb_y * fx

                     - 0.5 * fx * fz * fga * pb_zz

                     - pa_yzz * pb_y * fz * fgb

                     + 2.0 * fx * fx * fz * pb_zz

                     + 5.0 * pa_yzz * fz * pb_y * fx

                     + 20.0 * pa_yz * fz * fx * pb_yz

                     + 5.0 * fx * pa_zz * fz * pb_zz

                     - pa_y * fz * fga * pb_yzz

                     + 5.0 * pa_y * fz * fx * pb_yzz

                     + 12.0 * pa_yzz * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_zzz_s_0(double fx,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double s_0_0)
    {
        return s_0_0 * (1.5 * pa_yz * fx * fx

                     + 2.25 * pa_y * fx * fx * pb_z

                     + 1.5 * pa_yzz * pb_z * fx

                     + 3.0 * pa_yz * fx * pb_zz

                     + 0.5 * pa_y * fx * pb_zzz

                     + pa_yzz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_zzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_yz * fx * fz * fgb

                     + 12.0 * pa_yz * fz * fx * fx

                     + 18.0 * pa_y * fx * fx * fz * pb_z

                     - 1.5 * pa_y * fx * pb_z * fz * fgb

                     - 1.5 * pa_y * fz * fga * pb_z * fx

                     - 3.0 * pa_yzz * pb_z * fz * fgb

                     + 15.0 * pa_yzz * fz * pb_z * fx

                     + 30.0 * pa_yz * fz * fx * pb_zz

                     - pa_y * fz * fga * pb_zzz

                     + 5.0 * pa_y * fz * fx * pb_zzz

                     + 12.0 * pa_yzz * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xxx_s_0(double fx,
                                   double pa_z,
                                   double pa_zzz,
                                   double pb_x,
                                   double pb_xxx,
                                   double s_0_0)
    {
        return s_0_0 * (2.25 * pa_z * fx * fx * pb_x

                     + 1.5 * pa_zzz * pb_x * fx

                     + 1.5 * pa_z * fx * pb_xxx

                     + pa_zzz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xxx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zzz,
                                   double pb_x,
                                   double pb_xxx,
                                   double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_z * fx * pb_x * fz * fgb

                     - 4.5 * pa_z * fz * fga * pb_x * fx

                     - 3.0 * pa_zzz * pb_x * fz * fgb

                     + 18.0 * pa_z * fz * fx * fx * pb_x

                     + 15.0 * pa_zzz * fz * pb_x * fx

                     - 3.0 * pa_z * fz * fga * pb_xxx

                     + 15.0 * pa_z * fz * fx * pb_xxx

                     + 12.0 * pa_zzz * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xxy_s_0(double fx,
                                   double pa_z,
                                   double pa_zzz,
                                   double pb_xxy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_z * fx * fx * pb_y

                     + 0.5 * pa_zzz * fx * pb_y

                     + 1.5 * pa_z * fx * pb_xxy

                     + pa_zzz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xxy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zzz,
                                   double pb_xxy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_z * fx * fz * fgb * pb_y

                     - 1.5 * pa_z * fz * fga * fx * pb_y

                     - pa_zzz * fz * fgb * pb_y

                     + 6.0 * pa_z * fz * fx * fx * pb_y

                     + 5.0 * pa_zzz * fz * fx * pb_y

                     - 3.0 * pa_z * fz * fga * pb_xxy

                     + 15.0 * pa_z * fz * fx * pb_xxy

                     + 12.0 * pa_zzz * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xxz_s_0(double fx,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pb_xx,
                                   double pb_xxz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_zz * fx * fx

                     + 0.75 * pa_z * fx * fx * pb_z

                     + 0.75 * fx * fx * pb_xx

                     + 0.5 * pa_zzz * fx * pb_z

                     + 1.5 * pa_zz * fx * pb_xx

                     + 1.5 * pa_z * fx * pb_xxz

                     + pa_zzz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xxz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pb_xx,
                                   double pb_xxz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     - 0.75 * fx * fx * fz * fga

                     - 1.5 * pa_zz * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * pa_zz * fz * fx * fx

                     - 1.5 * pa_z * fx * fz * fgb * pb_z

                     - 1.5 * pa_z * fz * fga * fx * pb_z

                     - 1.5 * fx * fz * fga * pb_xx

                     - pa_zzz * fz * fgb * pb_z

                     + 6.0 * pa_z * fz * fx * fx * pb_z

                     + 6.0 * fx * fx * fz * pb_xx

                     + 5.0 * pa_zzz * fz * fx * pb_z

                     + 15.0 * pa_zz * fz * fx * pb_xx

                     - 3.0 * pa_z * fz * fga * pb_xxz

                     + 15.0 * pa_z * fz * fx * pb_xxz

                     + 12.0 * pa_zzz * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xyy_s_0(double fx,
                                   double pa_z,
                                   double pa_zzz,
                                   double pb_x,
                                   double pb_xyy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_z * fx * fx * pb_x

                     + 0.5 * pa_zzz * pb_x * fx

                     + 1.5 * pa_z * fx * pb_xyy

                     + pa_zzz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zzz,
                                   double pb_x,
                                   double pb_xyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_z * fx * pb_x * fz * fgb

                     - 1.5 * pa_z * fz * fga * pb_x * fx

                     - pa_zzz * pb_x * fz * fgb

                     + 6.0 * pa_z * fz * fx * fx * pb_x

                     + 5.0 * pa_zzz * fz * pb_x * fx

                     - 3.0 * pa_z * fz * fga * pb_xyy

                     + 15.0 * pa_z * fz * fx * pb_xyy

                     + 12.0 * pa_zzz * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xyz_s_0(double fx,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pb_xy,
                                   double pb_xyz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_xy

                     + 1.5 * pa_zz * fx * pb_xy

                     + 1.5 * pa_z * fx * pb_xyz

                     + pa_zzz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xyz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pb_xy,
                                   double pb_xyz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fz * fga * pb_xy

                     + 6.0 * fx * fx * fz * pb_xy

                     + 15.0 * pa_zz * fz * fx * pb_xy

                     - 3.0 * pa_z * fz * fga * pb_xyz

                     + 15.0 * pa_z * fz * fx * pb_xyz

                     + 12.0 * pa_zzz * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xzz_s_0(double fx,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_xzz,
                                   double s_0_0)
    {
        return s_0_0 * (2.25 * pa_z * fx * fx * pb_x

                     + 1.5 * fx * fx * pb_xz

                     + 0.5 * pa_zzz * pb_x * fx

                     + 3.0 * pa_zz * fx * pb_xz

                     + 1.5 * pa_z * fx * pb_xzz

                     + pa_zzz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_xzz,
                                   double r_0_0)
    {
        return r_0_0 * (18.0 * pa_z * fx * fx * fz * pb_x

                     - 1.5 * pa_z * fx * pb_x * fz * fgb

                     - 1.5 * pa_z * fz * fga * pb_x * fx

                     - 3.0 * fx * fz * fga * pb_xz

                     - pa_zzz * pb_x * fz * fgb

                     + 12.0 * fx * fx * fz * pb_xz

                     + 5.0 * pa_zzz * fz * pb_x * fx

                     + 30.0 * pa_zz * fz * fx * pb_xz

                     - 3.0 * pa_z * fz * fga * pb_xzz

                     + 15.0 * pa_z * fz * fx * pb_xzz

                     + 12.0 * pa_zzz * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_yyy_s_0(double fx,
                                   double pa_z,
                                   double pa_zzz,
                                   double pb_y,
                                   double pb_yyy,
                                   double s_0_0)
    {
        return s_0_0 * (2.25 * pa_z * fx * fx * pb_y

                     + 1.5 * pa_zzz * pb_y * fx

                     + 1.5 * pa_z * fx * pb_yyy

                     + pa_zzz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_yyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zzz,
                                   double pb_y,
                                   double pb_yyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_z * fx * pb_y * fz * fgb

                     - 4.5 * pa_z * fz * fga * pb_y * fx

                     - 3.0 * pa_zzz * pb_y * fz * fgb

                     + 18.0 * pa_z * fz * fx * fx * pb_y

                     + 15.0 * pa_zzz * fz * pb_y * fx

                     - 3.0 * pa_z * fz * fga * pb_yyy

                     + 15.0 * pa_z * fz * fx * pb_yyy

                     + 12.0 * pa_zzz * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_yyz_s_0(double fx,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pb_yy,
                                   double pb_yyz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_zz * fx * fx

                     + 0.75 * pa_z * fx * fx * pb_z

                     + 0.75 * fx * fx * pb_yy

                     + 0.5 * pa_zzz * fx * pb_z

                     + 1.5 * pa_zz * fx * pb_yy

                     + 1.5 * pa_z * fx * pb_yyz

                     + pa_zzz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_yyz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pb_yy,
                                   double pb_yyz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     - 0.75 * fx * fx * fz * fga

                     - 1.5 * pa_zz * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * pa_zz * fz * fx * fx

                     - 1.5 * pa_z * fx * fz * fgb * pb_z

                     - 1.5 * pa_z * fz * fga * fx * pb_z

                     - 1.5 * fx * fz * fga * pb_yy

                     - pa_zzz * fz * fgb * pb_z

                     + 6.0 * pa_z * fz * fx * fx * pb_z

                     + 6.0 * fx * fx * fz * pb_yy

                     + 5.0 * pa_zzz * fz * fx * pb_z

                     + 15.0 * pa_zz * fz * fx * pb_yy

                     - 3.0 * pa_z * fz * fga * pb_yyz

                     + 15.0 * pa_z * fz * fx * pb_yyz

                     + 12.0 * pa_zzz * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_yzz_s_0(double fx,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_yzz,
                                   double s_0_0)
    {
        return s_0_0 * (2.25 * pa_z * fx * fx * pb_y

                     + 1.5 * fx * fx * pb_yz

                     + 0.5 * pa_zzz * pb_y * fx

                     + 3.0 * pa_zz * fx * pb_yz

                     + 1.5 * pa_z * fx * pb_yzz

                     + pa_zzz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_yzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_yzz,
                                   double r_0_0)
    {
        return r_0_0 * (18.0 * pa_z * fx * fx * fz * pb_y

                     - 1.5 * pa_z * fx * pb_y * fz * fgb

                     - 1.5 * pa_z * fz * fga * pb_y * fx

                     - 3.0 * fx * fz * fga * pb_yz

                     - pa_zzz * pb_y * fz * fgb

                     + 12.0 * fx * fx * fz * pb_yz

                     + 5.0 * pa_zzz * fz * pb_y * fx

                     + 30.0 * pa_zz * fz * fx * pb_yz

                     - 3.0 * pa_z * fz * fga * pb_yzz

                     + 15.0 * pa_z * fz * fx * pb_yzz

                     + 12.0 * pa_zzz * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_zzz_s_0(double fx,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx

                     + 2.25 * pa_zz * fx * fx

                     + 6.75 * pa_z * fx * fx * pb_z

                     + 2.25 * fx * fx * pb_zz

                     + 1.5 * pa_zzz * pb_z * fx

                     + 4.5 * pa_zz * fx * pb_zz

                     + 1.5 * pa_z * fx * pb_zzz

                     + pa_zzz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_zzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double r_0_0)
    {
        return r_0_0 * (11.25 * fx * fx * fx * fz

                     - 2.25 * fx * fx * fz * fgb

                     - 2.25 * fx * fx * fz * fga

                     - 4.5 * pa_zz * fx * fz * fgb

                     + 18.0 * pa_zz * fz * fx * fx

                     + 54.0 * pa_z * fx * fx * fz * pb_z

                     - 4.5 * pa_z * fx * pb_z * fz * fgb

                     - 4.5 * pa_z * fz * fga * pb_z * fx

                     - 4.5 * fx * fz * fga * pb_zz

                     - 3.0 * pa_zzz * pb_z * fz * fgb

                     + 18.0 * fx * fx * fz * pb_zz

                     + 15.0 * pa_zzz * fz * pb_z * fx

                     + 45.0 * pa_zz * fz * fx * pb_zz

                     - 3.0 * pa_z * fz * fga * pb_zzz

                     + 15.0 * pa_z * fz * fx * pb_zzz

                     + 12.0 * pa_zzz * fz * pb_zzz);

    }


} // kinrecfunc namespace

