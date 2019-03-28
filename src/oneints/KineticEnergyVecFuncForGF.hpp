//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

namespace kinvecfunc { // kinvecfunc namespace

    // SIMD elementary functions for (G|T|F) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xxx_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxx,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double s_0_0)
    {
        return s_0_0 * (7.5 * pa_x * fx * fx * fx

                     + 5.625 * fx * fx * fx * pb_x

                     + 3.0 * pa_xxx * fx * fx

                     + 13.5 * pa_xx * fx * fx * pb_x

                     + 9.0 * pa_x * fx * fx * pb_xx

                     + 1.5 * pa_xxxx * pb_x * fx

                     + 6.0 * pa_xxx * fx * pb_xx

                     + 0.75 * fx * fx * pb_xxx

                     + 3.0 * pa_xx * fx * pb_xxx

                     + pa_xxxx * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxx,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double r_0_0)
    {
        return r_0_0 * (60.0 * pa_x * fx * fx * fx * fz

                     - 9.0 * pa_x * fx * fx * fz * fgb

                     - 9.0 * pa_x * fx * fx * fz * fga

                     - 13.5 * fx * fx * fz * fga * pb_x

                     - 6.0 * pa_xxx * fx * fz * fgb

                     + 45.0 * fx * fx * fx * fz * pb_x

                     + 30.0 * pa_xxx * fz * fx * fx

                     + 135.0 * pa_xx * fx * fx * fz * pb_x

                     - 2.25 * fx * fx * pb_x * fz * fgb

                     - 9.0 * pa_xx * fx * pb_x * fz * fgb

                     - 9.0 * pa_xx * fz * fga * pb_x * fx

                     - 18.0 * pa_x * fx * fz * fga * pb_xx

                     - 3.0 * pa_xxxx * pb_x * fz * fgb

                     + 90.0 * pa_x * fx * fx * fz * pb_xx

                     + 18.0 * pa_xxxx * fz * pb_x * fx

                     + 72.0 * pa_xxx * fz * fx * pb_xx

                     - 3.0 * fx * fz * fga * pb_xxx

                     - 6.0 * pa_xx * fz * fga * pb_xxx

                     + 7.5 * fx * fx * fz * pb_xxx

                     + 36.0 * pa_xx * fz * fx * pb_xxx

                     + 14.0 * pa_xxxx * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xxy_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxx,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pb_y

                     + 4.5 * pa_xx * fx * fx * pb_y

                     + 6.0 * pa_x * fx * fx * pb_xy

                     + 0.5 * pa_xxxx * fx * pb_y

                     + 4.0 * pa_xxx * fx * pb_xy

                     + 0.75 * fx * fx * pb_xxy

                     + 3.0 * pa_xx * fx * pb_xxy

                     + pa_xxxx * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxx,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * fx * fx * fz * fga * pb_y

                     + 15.0 * fx * fx * fx * fz * pb_y

                     + 45.0 * pa_xx * fx * fx * fz * pb_y

                     - 0.75 * fx * fx * fz * fgb * pb_y

                     - 3.0 * pa_xx * fx * fz * fgb * pb_y

                     - 3.0 * pa_xx * fz * fga * fx * pb_y

                     - 12.0 * pa_x * fx * fz * fga * pb_xy

                     - pa_xxxx * fz * fgb * pb_y

                     + 60.0 * pa_x * fx * fx * fz * pb_xy

                     + 6.0 * pa_xxxx * fz * fx * pb_y

                     + 48.0 * pa_xxx * fz * fx * pb_xy

                     - 3.0 * fx * fz * fga * pb_xxy

                     - 6.0 * pa_xx * fz * fga * pb_xxy

                     + 7.5 * fx * fx * fz * pb_xxy

                     + 36.0 * pa_xx * fz * fx * pb_xxy

                     + 14.0 * pa_xxxx * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xxz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxx,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pb_z

                     + 4.5 * pa_xx * fx * fx * pb_z

                     + 6.0 * pa_x * fx * fx * pb_xz

                     + 0.5 * pa_xxxx * fx * pb_z

                     + 4.0 * pa_xxx * fx * pb_xz

                     + 0.75 * fx * fx * pb_xxz

                     + 3.0 * pa_xx * fx * pb_xxz

                     + pa_xxxx * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxx,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * fx * fx * fz * fga * pb_z

                     + 15.0 * fx * fx * fx * fz * pb_z

                     + 45.0 * pa_xx * fx * fx * fz * pb_z

                     - 0.75 * fx * fx * fz * fgb * pb_z

                     - 3.0 * pa_xx * fx * fz * fgb * pb_z

                     - 3.0 * pa_xx * fz * fga * fx * pb_z

                     - 12.0 * pa_x * fx * fz * fga * pb_xz

                     - pa_xxxx * fz * fgb * pb_z

                     + 60.0 * pa_x * fx * fx * fz * pb_xz

                     + 6.0 * pa_xxxx * fz * fx * pb_z

                     + 48.0 * pa_xxx * fz * fx * pb_xz

                     - 3.0 * fx * fz * fga * pb_xxz

                     - 6.0 * pa_xx * fz * fga * pb_xxz

                     + 7.5 * fx * fx * fz * pb_xxz

                     + 36.0 * pa_xx * fz * fx * pb_xxz

                     + 14.0 * pa_xxxx * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xyy_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxx,
                                    double pb_x,
                                    double pb_xyy,
                                    double pb_yy,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_x * fx * fx * fx

                     + pa_xxx * fx * fx

                     + 0.375 * fx * fx * fx * pb_x

                     + 1.5 * pa_xx * fx * fx * pb_x

                     + 3.0 * pa_x * fx * fx * pb_yy

                     + 0.5 * pa_xxxx * pb_x * fx

                     + 2.0 * pa_xxx * fx * pb_yy

                     + 0.75 * fx * fx * pb_xyy

                     + 3.0 * pa_xx * fx * pb_xyy

                     + pa_xxxx * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxx,
                                    double pb_x,
                                    double pb_xyy,
                                    double pb_yy,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_x * fx * fx * fz * fgb

                     - 3.0 * pa_x * fx * fx * fz * fga

                     - 2.0 * pa_xxx * fx * fz * fgb

                     + 12.0 * pa_x * fx * fx * fx * fz

                     + 10.0 * pa_xxx * fz * fx * fx

                     - 0.75 * fx * fx * pb_x * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_x

                     - 3.0 * pa_xx * fx * pb_x * fz * fgb

                     - 3.0 * pa_xx * fz * fga * pb_x * fx

                     - 6.0 * pa_x * fx * fz * fga * pb_yy

                     + 3.0 * fx * fx * fx * fz * pb_x

                     - pa_xxxx * pb_x * fz * fgb

                     + 15.0 * pa_xx * fz * fx * fx * pb_x

                     + 30.0 * pa_x * fx * fx * fz * pb_yy

                     + 6.0 * pa_xxxx * fz * pb_x * fx

                     + 24.0 * pa_xxx * fz * fx * pb_yy

                     - 3.0 * fx * fz * fga * pb_xyy

                     - 6.0 * pa_xx * fz * fga * pb_xyy

                     + 7.5 * fx * fx * fz * pb_xyy

                     + 36.0 * pa_xx * fz * fx * pb_xyy

                     + 14.0 * pa_xxxx * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxx,
                                    double pb_xyz,
                                    double pb_yz,
                                    double s_0_0)
    {
        return s_0_0 * (3.0 * pa_x * fx * fx * pb_yz

                     + 2.0 * pa_xxx * fx * pb_yz

                     + 0.75 * fx * fx * pb_xyz

                     + 3.0 * pa_xx * fx * pb_xyz

                     + pa_xxxx * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xyz_r_0(double fga,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxx,
                                    double pb_xyz,
                                    double pb_yz,
                                    double r_0_0)
    {
        return r_0_0 * (- 6.0 * pa_x * fx * fz * fga * pb_yz

                     + 30.0 * pa_x * fx * fx * fz * pb_yz

                     + 24.0 * pa_xxx * fz * fx * pb_yz

                     - 3.0 * fx * fz * fga * pb_xyz

                     - 6.0 * pa_xx * fz * fga * pb_xyz

                     + 7.5 * fx * fx * fz * pb_xyz

                     + 36.0 * pa_xx * fz * fx * pb_xyz

                     + 14.0 * pa_xxxx * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxx,
                                    double pb_x,
                                    double pb_xzz,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_x * fx * fx * fx

                     + pa_xxx * fx * fx

                     + 0.375 * fx * fx * fx * pb_x

                     + 1.5 * pa_xx * fx * fx * pb_x

                     + 3.0 * pa_x * fx * fx * pb_zz

                     + 0.5 * pa_xxxx * pb_x * fx

                     + 2.0 * pa_xxx * fx * pb_zz

                     + 0.75 * fx * fx * pb_xzz

                     + 3.0 * pa_xx * fx * pb_xzz

                     + pa_xxxx * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxx,
                                    double pb_x,
                                    double pb_xzz,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_x * fx * fx * fz * fgb

                     - 3.0 * pa_x * fx * fx * fz * fga

                     - 2.0 * pa_xxx * fx * fz * fgb

                     + 12.0 * pa_x * fx * fx * fx * fz

                     + 10.0 * pa_xxx * fz * fx * fx

                     - 0.75 * fx * fx * pb_x * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_x

                     - 3.0 * pa_xx * fx * pb_x * fz * fgb

                     - 3.0 * pa_xx * fz * fga * pb_x * fx

                     - 6.0 * pa_x * fx * fz * fga * pb_zz

                     + 3.0 * fx * fx * fx * fz * pb_x

                     - pa_xxxx * pb_x * fz * fgb

                     + 15.0 * pa_xx * fz * fx * fx * pb_x

                     + 30.0 * pa_x * fx * fx * fz * pb_zz

                     + 6.0 * pa_xxxx * fz * pb_x * fx

                     + 24.0 * pa_xxx * fz * fx * pb_zz

                     - 3.0 * fx * fz * fga * pb_xzz

                     - 6.0 * pa_xx * fz * fga * pb_xzz

                     + 7.5 * fx * fx * fz * pb_xzz

                     + 36.0 * pa_xx * fz * fx * pb_xzz

                     + 14.0 * pa_xxxx * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_yyy_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxxx,
                                    double pb_y,
                                    double pb_yyy,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pb_y

                     + 4.5 * pa_xx * fx * fx * pb_y

                     + 1.5 * pa_xxxx * pb_y * fx

                     + 0.75 * fx * fx * pb_yyy

                     + 3.0 * pa_xx * fx * pb_yyy

                     + pa_xxxx * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_yyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxxx,
                                    double pb_y,
                                    double pb_yyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pb_y * fz * fgb

                     - 4.5 * fx * fx * fz * fga * pb_y

                     - 9.0 * pa_xx * fx * pb_y * fz * fgb

                     - 9.0 * pa_xx * fz * fga * pb_y * fx

                     + 9.0 * fx * fx * fx * fz * pb_y

                     - 3.0 * pa_xxxx * pb_y * fz * fgb

                     + 45.0 * pa_xx * fz * fx * fx * pb_y

                     + 18.0 * pa_xxxx * fz * pb_y * fx

                     - 3.0 * fx * fz * fga * pb_yyy

                     - 6.0 * pa_xx * fz * fga * pb_yyy

                     + 7.5 * fx * fx * fz * pb_yyy

                     + 36.0 * pa_xx * fz * fx * pb_yyy

                     + 14.0 * pa_xxxx * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_yyz_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxxx,
                                    double pb_yyz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_z

                     + 1.5 * pa_xx * fx * fx * pb_z

                     + 0.5 * pa_xxxx * fx * pb_z

                     + 0.75 * fx * fx * pb_yyz

                     + 3.0 * pa_xx * fx * pb_yyz

                     + pa_xxxx * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_yyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxxx,
                                    double pb_yyz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb * pb_z

                     - 1.5 * fx * fx * fz * fga * pb_z

                     - 3.0 * pa_xx * fx * fz * fgb * pb_z

                     - 3.0 * pa_xx * fz * fga * fx * pb_z

                     + 3.0 * fx * fx * fx * fz * pb_z

                     - pa_xxxx * fz * fgb * pb_z

                     + 15.0 * pa_xx * fz * fx * fx * pb_z

                     + 6.0 * pa_xxxx * fz * fx * pb_z

                     - 3.0 * fx * fz * fga * pb_yyz

                     - 6.0 * pa_xx * fz * fga * pb_yyz

                     + 7.5 * fx * fx * fz * pb_yyz

                     + 36.0 * pa_xx * fz * fx * pb_yyz

                     + 14.0 * pa_xxxx * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_yzz_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxxx,
                                    double pb_y,
                                    double pb_yzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_y

                     + 1.5 * pa_xx * fx * fx * pb_y

                     + 0.5 * pa_xxxx * pb_y * fx

                     + 0.75 * fx * fx * pb_yzz

                     + 3.0 * pa_xx * fx * pb_yzz

                     + pa_xxxx * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_yzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxxx,
                                    double pb_y,
                                    double pb_yzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pb_y * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_y

                     - 3.0 * pa_xx * fx * pb_y * fz * fgb

                     - 3.0 * pa_xx * fz * fga * pb_y * fx

                     + 3.0 * fx * fx * fx * fz * pb_y

                     - pa_xxxx * pb_y * fz * fgb

                     + 15.0 * pa_xx * fz * fx * fx * pb_y

                     + 6.0 * pa_xxxx * fz * pb_y * fx

                     - 3.0 * fx * fz * fga * pb_yzz

                     - 6.0 * pa_xx * fz * fga * pb_yzz

                     + 7.5 * fx * fx * fz * pb_yzz

                     + 36.0 * pa_xx * fz * fx * pb_yzz

                     + 14.0 * pa_xxxx * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_zzz_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxxx,
                                    double pb_z,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pb_z

                     + 4.5 * pa_xx * fx * fx * pb_z

                     + 1.5 * pa_xxxx * pb_z * fx

                     + 0.75 * fx * fx * pb_zzz

                     + 3.0 * pa_xx * fx * pb_zzz

                     + pa_xxxx * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_zzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxxx,
                                    double pb_z,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pb_z * fz * fgb

                     - 4.5 * fx * fx * fz * fga * pb_z

                     - 9.0 * pa_xx * fx * pb_z * fz * fgb

                     - 9.0 * pa_xx * fz * fga * pb_z * fx

                     + 9.0 * fx * fx * fx * fz * pb_z

                     - 3.0 * pa_xxxx * pb_z * fz * fgb

                     + 45.0 * pa_xx * fz * fx * fx * pb_z

                     + 18.0 * pa_xxxx * fz * pb_z * fx

                     - 3.0 * fx * fz * fga * pb_zzz

                     - 6.0 * pa_xx * fz * fga * pb_zzz

                     + 7.5 * fx * fx * fz * pb_zzz

                     + 36.0 * pa_xx * fz * fx * pb_zzz

                     + 14.0 * pa_xxxx * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xxx_s_0(double fx,
                                    double pa_xxxy,
                                    double pa_xxy,
                                    double pa_xy,
                                    double pa_y,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pa_y

                     + 2.25 * pa_xxy * fx * fx

                     + 6.75 * pa_xy * fx * fx * pb_x

                     + 2.25 * fx * fx * pa_y * pb_xx

                     + 1.5 * pa_xxxy * pb_x * fx

                     + 4.5 * pa_xxy * fx * pb_xx

                     + 1.5 * pa_xy * fx * pb_xxx

                     + pa_xxxy * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxxy,
                                    double pa_xxy,
                                    double pa_xy,
                                    double pa_y,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double r_0_0)
    {
        return r_0_0 * (15.0 * fx * fx * fx * fz * pa_y

                     - 2.25 * fx * fx * pa_y * fz * fgb

                     - 2.25 * fx * fx * fz * fga * pa_y

                     - 4.5 * pa_xxy * fx * fz * fgb

                     + 22.5 * pa_xxy * fx * fx * fz

                     + 67.5 * pa_xy * fx * fx * fz * pb_x

                     - 4.5 * pa_xy * fx * pb_x * fz * fgb

                     - 4.5 * pa_xy * fz * fga * pb_x * fx

                     - 4.5 * fx * fz * fga * pa_y * pb_xx

                     - 3.0 * pa_xxxy * pb_x * fz * fgb

                     + 22.5 * fx * fx * fz * pa_y * pb_xx

                     + 18.0 * pa_xxxy * fz * pb_x * fx

                     + 54.0 * pa_xxy * fx * fz * pb_xx

                     - 3.0 * pa_xy * fz * fga * pb_xxx

                     + 18.0 * pa_xy * fx * fz * pb_xxx

                     + 14.0 * pa_xxxy * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xxy_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxy,
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
        return s_0_0 * (1.125 * pa_x * fx * fx * fx

                     + 0.75 * fx * fx * fx * pb_x

                     + 0.25 * pa_xxx * fx * fx

                     + 1.5 * pa_xx * fx * fx * pb_x

                     + 2.25 * pa_xy * fx * fx * pb_y

                     + 0.75 * pa_x * fx * fx * pb_xx

                     + 1.5 * fx * fx * pa_y * pb_xy

                     + 0.5 * pa_xxxy * fx * pb_y

                     + 0.5 * pa_xxx * fx * pb_xx

                     + 3.0 * pa_xxy * fx * pb_xy

                     + 1.5 * pa_xy * fx * pb_xxy

                     + pa_xxxy * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxy,
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
        return r_0_0 * (9.0 * pa_x * fx * fx * fx * fz

                     - 0.75 * pa_x * fx * fx * fz * fgb

                     - 0.75 * pa_x * fz * fga * fx * fx

                     - 1.5 * fx * fx * fz * fga * pb_x

                     - 0.5 * pa_xxx * fx * fz * fgb

                     + 6.0 * fx * fx * fx * fz * pb_x

                     + 2.5 * pa_xxx * fz * fx * fx

                     + 15.0 * pa_xx * fx * fx * fz * pb_x

                     + 22.5 * pa_xy * fx * fx * fz * pb_y

                     - 1.5 * pa_xy * fx * fz * fgb * pb_y

                     - 1.5 * pa_xy * fz * fga * fx * pb_y

                     - 1.5 * pa_x * fz * fga * fx * pb_xx

                     - 3.0 * fx * fz * fga * pa_y * pb_xy

                     - pa_xxxy * fz * fgb * pb_y

                     + 7.5 * pa_x * fx * fx * fz * pb_xx

                     + 15.0 * fx * fx * fz * pa_y * pb_xy

                     + 6.0 * pa_xxxy * fz * fx * pb_y

                     + 6.0 * pa_xxx * fz * fx * pb_xx

                     + 36.0 * pa_xxy * fx * fz * pb_xy

                     - 3.0 * pa_xy * fz * fga * pb_xxy

                     + 18.0 * pa_xy * fx * fz * pb_xxy

                     + 14.0 * pa_xxxy * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xxz_s_0(double fx,
                                    double pa_xxxy,
                                    double pa_xxy,
                                    double pa_xy,
                                    double pa_y,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_xy * fx * fx * pb_z

                     + 1.5 * fx * fx * pa_y * pb_xz

                     + 0.5 * pa_xxxy * fx * pb_z

                     + 3.0 * pa_xxy * fx * pb_xz

                     + 1.5 * pa_xy * fx * pb_xxz

                     + pa_xxxy * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxxy,
                                    double pa_xxy,
                                    double pa_xy,
                                    double pa_y,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (22.5 * pa_xy * fx * fx * fz * pb_z

                     - 1.5 * pa_xy * fx * fz * fgb * pb_z

                     - 1.5 * pa_xy * fz * fga * fx * pb_z

                     - 3.0 * fx * fz * fga * pa_y * pb_xz

                     - pa_xxxy * fz * fgb * pb_z

                     + 15.0 * fx * fx * fz * pa_y * pb_xz

                     + 6.0 * pa_xxxy * fz * fx * pb_z

                     + 36.0 * pa_xxy * fx * fz * pb_xz

                     - 3.0 * pa_xy * fz * fga * pb_xxz

                     + 18.0 * pa_xy * fx * fz * pb_xxz

                     + 14.0 * pa_xxxy * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xyy_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxy,
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
        return s_0_0 * (0.375 * fx * fx * fx * pa_y

                     + 0.75 * fx * fx * fx * pb_y

                     + 0.75 * pa_xxy * fx * fx

                     + 1.5 * pa_xx * fx * fx * pb_y

                     + 0.75 * pa_xy * fx * fx * pb_x

                     + 1.5 * pa_x * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_y * pb_yy

                     + 0.5 * pa_xxxy * pb_x * fx

                     + pa_xxx * fx * pb_xy

                     + 1.5 * pa_xxy * fx * pb_yy

                     + 1.5 * pa_xy * fx * pb_xyy

                     + pa_xxxy * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxy,
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
        return r_0_0 * (- 0.75 * fx * fx * pa_y * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pa_y

                     - 1.5 * fx * fx * fz * fga * pb_y

                     - 1.5 * pa_xxy * fx * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_y

                     + 6.0 * fx * fx * fx * fz * pb_y

                     + 7.5 * pa_xxy * fx * fx * fz

                     + 15.0 * pa_xx * fx * fx * fz * pb_y

                     - 1.5 * pa_xy * fx * pb_x * fz * fgb

                     - 1.5 * pa_xy * fz * fga * pb_x * fx

                     - 3.0 * pa_x * fz * fga * fx * pb_xy

                     - 1.5 * fx * fz * fga * pa_y * pb_yy

                     - pa_xxxy * pb_x * fz * fgb

                     + 7.5 * pa_xy * fx * fx * fz * pb_x

                     + 15.0 * pa_x * fx * fx * fz * pb_xy

                     + 7.5 * fx * fx * fz * pa_y * pb_yy

                     + 6.0 * pa_xxxy * fz * pb_x * fx

                     + 12.0 * pa_xxx * fz * fx * pb_xy

                     + 18.0 * pa_xxy * fx * fz * pb_yy

                     - 3.0 * pa_xy * fz * fga * pb_xyy

                     + 18.0 * pa_xy * fx * fz * pb_xyy

                     + 14.0 * pa_xxxy * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxy,
                                    double pa_xxy,
                                    double pa_xy,
                                    double pa_y,
                                    double pb_xyz,
                                    double pb_xz,
                                    double pb_yz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_z

                     + 0.75 * pa_xx * fx * fx * pb_z

                     + 0.75 * pa_x * fx * fx * pb_xz

                     + 0.75 * fx * fx * pa_y * pb_yz

                     + 0.5 * pa_xxx * fx * pb_xz

                     + 1.5 * pa_xxy * fx * pb_yz

                     + 1.5 * pa_xy * fx * pb_xyz

                     + pa_xxxy * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xyz_r_0(double fga,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxy,
                                    double pa_xxy,
                                    double pa_xy,
                                    double pa_y,
                                    double pb_xyz,
                                    double pb_xz,
                                    double pb_yz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fga * pb_z

                     + 3.0 * fx * fx * fx * fz * pb_z

                     + 7.5 * pa_xx * fx * fx * fz * pb_z

                     - 1.5 * pa_x * fz * fga * fx * pb_xz

                     - 1.5 * fx * fz * fga * pa_y * pb_yz

                     + 7.5 * pa_x * fx * fx * fz * pb_xz

                     + 7.5 * fx * fx * fz * pa_y * pb_yz

                     + 6.0 * pa_xxx * fz * fx * pb_xz

                     + 18.0 * pa_xxy * fx * fz * pb_yz

                     - 3.0 * pa_xy * fz * fga * pb_xyz

                     + 18.0 * pa_xy * fx * fz * pb_xyz

                     + 14.0 * pa_xxxy * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xzz_s_0(double fx,
                                    double pa_xxxy,
                                    double pa_xxy,
                                    double pa_xy,
                                    double pa_y,
                                    double pb_x,
                                    double pb_xzz,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_y

                     + 0.75 * pa_xxy * fx * fx

                     + 0.75 * pa_xy * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_y * pb_zz

                     + 0.5 * pa_xxxy * pb_x * fx

                     + 1.5 * pa_xxy * fx * pb_zz

                     + 1.5 * pa_xy * fx * pb_xzz

                     + pa_xxxy * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxxy,
                                    double pa_xxy,
                                    double pa_xy,
                                    double pa_y,
                                    double pb_x,
                                    double pb_xzz,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pa_y * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pa_y

                     - 1.5 * pa_xxy * fx * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_y

                     + 7.5 * pa_xxy * fx * fx * fz

                     - 1.5 * pa_xy * fx * pb_x * fz * fgb

                     - 1.5 * pa_xy * fz * fga * pb_x * fx

                     - 1.5 * fx * fz * fga * pa_y * pb_zz

                     - pa_xxxy * pb_x * fz * fgb

                     + 7.5 * pa_xy * fx * fx * fz * pb_x

                     + 7.5 * fx * fx * fz * pa_y * pb_zz

                     + 6.0 * pa_xxxy * fz * pb_x * fx

                     + 18.0 * pa_xxy * fx * fz * pb_zz

                     - 3.0 * pa_xy * fz * fga * pb_xzz

                     + 18.0 * pa_xy * fx * fz * pb_xzz

                     + 14.0 * pa_xxxy * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_yyy_s_0(double fx,
                                    double pa_x,
                                    double pa_xxx,
                                    double pa_xxxy,
                                    double pa_xy,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * pa_x * fx * fx * fx

                     + 0.75 * pa_xxx * fx * fx

                     + 2.25 * pa_xy * fx * fx * pb_y

                     + 2.25 * pa_x * fx * fx * pb_yy

                     + 1.5 * pa_xxxy * pb_y * fx

                     + 1.5 * pa_xxx * fx * pb_yy

                     + 1.5 * pa_xy * fx * pb_yyy

                     + pa_xxxy * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_yyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xxx,
                                    double pa_xxxy,
                                    double pa_xy,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 2.25 * pa_x * fx * fx * fz * fgb

                     - 2.25 * pa_x * fz * fga * fx * fx

                     - 1.5 * pa_xxx * fx * fz * fgb

                     + 9.0 * pa_x * fx * fx * fx * fz

                     + 7.5 * pa_xxx * fz * fx * fx

                     - 4.5 * pa_xy * fx * pb_y * fz * fgb

                     - 4.5 * pa_xy * fz * fga * pb_y * fx

                     - 4.5 * pa_x * fz * fga * fx * pb_yy

                     - 3.0 * pa_xxxy * pb_y * fz * fgb

                     + 22.5 * pa_xy * fx * fx * fz * pb_y

                     + 22.5 * pa_x * fx * fx * fz * pb_yy

                     + 18.0 * pa_xxxy * fz * pb_y * fx

                     + 18.0 * pa_xxx * fz * fx * pb_yy

                     - 3.0 * pa_xy * fz * fga * pb_yyy

                     + 18.0 * pa_xy * fx * fz * pb_yyy

                     + 14.0 * pa_xxxy * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_yyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xxx,
                                    double pa_xxxy,
                                    double pa_xy,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xy * fx * fx * pb_z

                     + 1.5 * pa_x * fx * fx * pb_yz

                     + 0.5 * pa_xxxy * fx * pb_z

                     + pa_xxx * fx * pb_yz

                     + 1.5 * pa_xy * fx * pb_yyz

                     + pa_xxxy * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_yyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xxx,
                                    double pa_xxxy,
                                    double pa_xy,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xy * fx * fz * fgb * pb_z

                     - 1.5 * pa_xy * fz * fga * fx * pb_z

                     - 3.0 * pa_x * fz * fga * fx * pb_yz

                     - pa_xxxy * fz * fgb * pb_z

                     + 7.5 * pa_xy * fx * fx * fz * pb_z

                     + 15.0 * pa_x * fx * fx * fz * pb_yz

                     + 6.0 * pa_xxxy * fz * fx * pb_z

                     + 12.0 * pa_xxx * fz * fx * pb_yz

                     - 3.0 * pa_xy * fz * fga * pb_yyz

                     + 18.0 * pa_xy * fx * fz * pb_yyz

                     + 14.0 * pa_xxxy * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_yzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xxx,
                                    double pa_xxxy,
                                    double pa_xy,
                                    double pb_y,
                                    double pb_yzz,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 0.25 * pa_xxx * fx * fx

                     + 0.75 * pa_xy * fx * fx * pb_y

                     + 0.75 * pa_x * fx * fx * pb_zz

                     + 0.5 * pa_xxxy * pb_y * fx

                     + 0.5 * pa_xxx * fx * pb_zz

                     + 1.5 * pa_xy * fx * pb_yzz

                     + pa_xxxy * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_yzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xxx,
                                    double pa_xxxy,
                                    double pa_xy,
                                    double pb_y,
                                    double pb_yzz,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * pa_x * fx * fx * fz * fgb

                     - 0.75 * pa_x * fz * fga * fx * fx

                     - 0.5 * pa_xxx * fx * fz * fgb

                     + 3.0 * pa_x * fx * fx * fx * fz

                     + 2.5 * pa_xxx * fz * fx * fx

                     - 1.5 * pa_xy * fx * pb_y * fz * fgb

                     - 1.5 * pa_xy * fz * fga * pb_y * fx

                     - 1.5 * pa_x * fz * fga * fx * pb_zz

                     - pa_xxxy * pb_y * fz * fgb

                     + 7.5 * pa_xy * fx * fx * fz * pb_y

                     + 7.5 * pa_x * fx * fx * fz * pb_zz

                     + 6.0 * pa_xxxy * fz * pb_y * fx

                     + 6.0 * pa_xxx * fz * fx * pb_zz

                     - 3.0 * pa_xy * fz * fga * pb_yzz

                     + 18.0 * pa_xy * fx * fz * pb_yzz

                     + 14.0 * pa_xxxy * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_zzz_s_0(double fx,
                                    double pa_xxxy,
                                    double pa_xy,
                                    double pb_z,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_xy * fx * fx * pb_z

                     + 1.5 * pa_xxxy * pb_z * fx

                     + 1.5 * pa_xy * fx * pb_zzz

                     + pa_xxxy * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_zzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxxy,
                                    double pa_xy,
                                    double pb_z,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_xy * fx * pb_z * fz * fgb

                     - 4.5 * pa_xy * fz * fga * pb_z * fx

                     - 3.0 * pa_xxxy * pb_z * fz * fgb

                     + 22.5 * pa_xy * fx * fx * fz * pb_z

                     + 18.0 * pa_xxxy * fz * pb_z * fx

                     - 3.0 * pa_xy * fz * fga * pb_zzz

                     + 18.0 * pa_xy * fx * fz * pb_zzz

                     + 14.0 * pa_xxxy * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xxx_s_0(double fx,
                                    double pa_xxxz,
                                    double pa_xxz,
                                    double pa_xz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pa_z

                     + 2.25 * pa_xxz * fx * fx

                     + 6.75 * pa_xz * fx * fx * pb_x

                     + 2.25 * fx * fx * pa_z * pb_xx

                     + 1.5 * pa_xxxz * pb_x * fx

                     + 4.5 * pa_xxz * fx * pb_xx

                     + 1.5 * pa_xz * fx * pb_xxx

                     + pa_xxxz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxxz,
                                    double pa_xxz,
                                    double pa_xz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double r_0_0)
    {
        return r_0_0 * (15.0 * fx * fx * fx * fz * pa_z

                     - 2.25 * fx * fx * pa_z * fz * fgb

                     - 2.25 * fx * fx * fz * fga * pa_z

                     - 4.5 * pa_xxz * fx * fz * fgb

                     + 22.5 * pa_xxz * fx * fx * fz

                     + 67.5 * pa_xz * fx * fx * fz * pb_x

                     - 4.5 * pa_xz * fx * pb_x * fz * fgb

                     - 4.5 * pa_xz * fz * fga * pb_x * fx

                     - 4.5 * fx * fz * fga * pa_z * pb_xx

                     - 3.0 * pa_xxxz * pb_x * fz * fgb

                     + 22.5 * fx * fx * fz * pa_z * pb_xx

                     + 18.0 * pa_xxxz * fz * pb_x * fx

                     + 54.0 * pa_xxz * fx * fz * pb_xx

                     - 3.0 * pa_xz * fz * fga * pb_xxx

                     + 18.0 * pa_xz * fx * fz * pb_xxx

                     + 14.0 * pa_xxxz * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xxy_s_0(double fx,
                                    double pa_xxxz,
                                    double pa_xxz,
                                    double pa_xz,
                                    double pa_z,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_xz * fx * fx * pb_y

                     + 1.5 * fx * fx * pa_z * pb_xy

                     + 0.5 * pa_xxxz * fx * pb_y

                     + 3.0 * pa_xxz * fx * pb_xy

                     + 1.5 * pa_xz * fx * pb_xxy

                     + pa_xxxz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxxz,
                                    double pa_xxz,
                                    double pa_xz,
                                    double pa_z,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double r_0_0)
    {
        return r_0_0 * (22.5 * pa_xz * fx * fx * fz * pb_y

                     - 1.5 * pa_xz * fx * fz * fgb * pb_y

                     - 1.5 * pa_xz * fz * fga * fx * pb_y

                     - 3.0 * fx * fz * fga * pa_z * pb_xy

                     - pa_xxxz * fz * fgb * pb_y

                     + 15.0 * fx * fx * fz * pa_z * pb_xy

                     + 6.0 * pa_xxxz * fz * fx * pb_y

                     + 36.0 * pa_xxz * fx * fz * pb_xy

                     - 3.0 * pa_xz * fz * fga * pb_xxy

                     + 18.0 * pa_xz * fx * fz * pb_xxy

                     + 14.0 * pa_xxxz * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xxz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxz,
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
        return s_0_0 * (1.125 * pa_x * fx * fx * fx

                     + 0.75 * fx * fx * fx * pb_x

                     + 0.25 * pa_xxx * fx * fx

                     + 1.5 * pa_xx * fx * fx * pb_x

                     + 2.25 * pa_xz * fx * fx * pb_z

                     + 0.75 * pa_x * fx * fx * pb_xx

                     + 1.5 * fx * fx * pa_z * pb_xz

                     + 0.5 * pa_xxxz * fx * pb_z

                     + 0.5 * pa_xxx * fx * pb_xx

                     + 3.0 * pa_xxz * fx * pb_xz

                     + 1.5 * pa_xz * fx * pb_xxz

                     + pa_xxxz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxz,
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
        return r_0_0 * (9.0 * pa_x * fx * fx * fx * fz

                     - 0.75 * pa_x * fx * fx * fz * fgb

                     - 0.75 * pa_x * fz * fga * fx * fx

                     - 1.5 * fx * fx * fz * fga * pb_x

                     - 0.5 * pa_xxx * fx * fz * fgb

                     + 6.0 * fx * fx * fx * fz * pb_x

                     + 2.5 * pa_xxx * fz * fx * fx

                     + 15.0 * pa_xx * fx * fx * fz * pb_x

                     + 22.5 * pa_xz * fx * fx * fz * pb_z

                     - 1.5 * pa_xz * fx * fz * fgb * pb_z

                     - 1.5 * pa_xz * fz * fga * fx * pb_z

                     - 1.5 * pa_x * fz * fga * fx * pb_xx

                     - 3.0 * fx * fz * fga * pa_z * pb_xz

                     - pa_xxxz * fz * fgb * pb_z

                     + 7.5 * pa_x * fx * fx * fz * pb_xx

                     + 15.0 * fx * fx * fz * pa_z * pb_xz

                     + 6.0 * pa_xxxz * fz * fx * pb_z

                     + 6.0 * pa_xxx * fz * fx * pb_xx

                     + 36.0 * pa_xxz * fx * fz * pb_xz

                     - 3.0 * pa_xz * fz * fga * pb_xxz

                     + 18.0 * pa_xz * fx * fz * pb_xxz

                     + 14.0 * pa_xxxz * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xyy_s_0(double fx,
                                    double pa_xxxz,
                                    double pa_xxz,
                                    double pa_xz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xyy,
                                    double pb_yy,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_z

                     + 0.75 * pa_xxz * fx * fx

                     + 0.75 * pa_xz * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_z * pb_yy

                     + 0.5 * pa_xxxz * pb_x * fx

                     + 1.5 * pa_xxz * fx * pb_yy

                     + 1.5 * pa_xz * fx * pb_xyy

                     + pa_xxxz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxxz,
                                    double pa_xxz,
                                    double pa_xz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xyy,
                                    double pb_yy,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pa_z * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pa_z

                     - 1.5 * pa_xxz * fx * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_z

                     + 7.5 * pa_xxz * fx * fx * fz

                     - 1.5 * pa_xz * fx * pb_x * fz * fgb

                     - 1.5 * pa_xz * fz * fga * pb_x * fx

                     - 1.5 * fx * fz * fga * pa_z * pb_yy

                     - pa_xxxz * pb_x * fz * fgb

                     + 7.5 * pa_xz * fx * fx * fz * pb_x

                     + 7.5 * fx * fx * fz * pa_z * pb_yy

                     + 6.0 * pa_xxxz * fz * pb_x * fx

                     + 18.0 * pa_xxz * fx * fz * pb_yy

                     - 3.0 * pa_xz * fz * fga * pb_xyy

                     + 18.0 * pa_xz * fx * fz * pb_xyy

                     + 14.0 * pa_xxxz * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxz,
                                    double pa_xxz,
                                    double pa_xz,
                                    double pa_z,
                                    double pb_xy,
                                    double pb_xyz,
                                    double pb_y,
                                    double pb_yz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_y

                     + 0.75 * pa_xx * fx * fx * pb_y

                     + 0.75 * pa_x * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_z * pb_yz

                     + 0.5 * pa_xxx * fx * pb_xy

                     + 1.5 * pa_xxz * fx * pb_yz

                     + 1.5 * pa_xz * fx * pb_xyz

                     + pa_xxxz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xyz_r_0(double fga,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxz,
                                    double pa_xxz,
                                    double pa_xz,
                                    double pa_z,
                                    double pb_xy,
                                    double pb_xyz,
                                    double pb_y,
                                    double pb_yz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fga * pb_y

                     + 3.0 * fx * fx * fx * fz * pb_y

                     + 7.5 * pa_xx * fx * fx * fz * pb_y

                     - 1.5 * pa_x * fz * fga * fx * pb_xy

                     - 1.5 * fx * fz * fga * pa_z * pb_yz

                     + 7.5 * pa_x * fx * fx * fz * pb_xy

                     + 7.5 * fx * fx * fz * pa_z * pb_yz

                     + 6.0 * pa_xxx * fz * fx * pb_xy

                     + 18.0 * pa_xxz * fx * fz * pb_yz

                     - 3.0 * pa_xz * fz * fga * pb_xyz

                     + 18.0 * pa_xz * fx * fz * pb_xyz

                     + 14.0 * pa_xxxz * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxz,
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
        return s_0_0 * (0.375 * fx * fx * fx * pa_z

                     + 0.75 * fx * fx * fx * pb_z

                     + 0.75 * pa_xxz * fx * fx

                     + 1.5 * pa_xx * fx * fx * pb_z

                     + 0.75 * pa_xz * fx * fx * pb_x

                     + 1.5 * pa_x * fx * fx * pb_xz

                     + 0.75 * fx * fx * pa_z * pb_zz

                     + 0.5 * pa_xxxz * pb_x * fx

                     + pa_xxx * fx * pb_xz

                     + 1.5 * pa_xxz * fx * pb_zz

                     + 1.5 * pa_xz * fx * pb_xzz

                     + pa_xxxz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pa_xxxz,
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
        return r_0_0 * (- 0.75 * fx * fx * pa_z * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pa_z

                     - 1.5 * fx * fx * fz * fga * pb_z

                     - 1.5 * pa_xxz * fx * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_z

                     + 6.0 * fx * fx * fx * fz * pb_z

                     + 7.5 * pa_xxz * fx * fx * fz

                     + 15.0 * pa_xx * fx * fx * fz * pb_z

                     - 1.5 * pa_xz * fx * pb_x * fz * fgb

                     - 1.5 * pa_xz * fz * fga * pb_x * fx

                     - 3.0 * pa_x * fz * fga * fx * pb_xz

                     - 1.5 * fx * fz * fga * pa_z * pb_zz

                     - pa_xxxz * pb_x * fz * fgb

                     + 7.5 * pa_xz * fx * fx * fz * pb_x

                     + 15.0 * pa_x * fx * fx * fz * pb_xz

                     + 7.5 * fx * fx * fz * pa_z * pb_zz

                     + 6.0 * pa_xxxz * fz * pb_x * fx

                     + 12.0 * pa_xxx * fz * fx * pb_xz

                     + 18.0 * pa_xxz * fx * fz * pb_zz

                     - 3.0 * pa_xz * fz * fga * pb_xzz

                     + 18.0 * pa_xz * fx * fz * pb_xzz

                     + 14.0 * pa_xxxz * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_yyy_s_0(double fx,
                                    double pa_xxxz,
                                    double pa_xz,
                                    double pb_y,
                                    double pb_yyy,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_xz * fx * fx * pb_y

                     + 1.5 * pa_xxxz * pb_y * fx

                     + 1.5 * pa_xz * fx * pb_yyy

                     + pa_xxxz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_yyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxxz,
                                    double pa_xz,
                                    double pb_y,
                                    double pb_yyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_xz * fx * pb_y * fz * fgb

                     - 4.5 * pa_xz * fz * fga * pb_y * fx

                     - 3.0 * pa_xxxz * pb_y * fz * fgb

                     + 22.5 * pa_xz * fx * fx * fz * pb_y

                     + 18.0 * pa_xxxz * fz * pb_y * fx

                     - 3.0 * pa_xz * fz * fga * pb_yyy

                     + 18.0 * pa_xz * fx * fz * pb_yyy

                     + 14.0 * pa_xxxz * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_yyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xxx,
                                    double pa_xxxz,
                                    double pa_xz,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 0.25 * pa_xxx * fx * fx

                     + 0.75 * pa_xz * fx * fx * pb_z

                     + 0.75 * pa_x * fx * fx * pb_yy

                     + 0.5 * pa_xxxz * fx * pb_z

                     + 0.5 * pa_xxx * fx * pb_yy

                     + 1.5 * pa_xz * fx * pb_yyz

                     + pa_xxxz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_yyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xxx,
                                    double pa_xxxz,
                                    double pa_xz,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * pa_x * fx * fx * fz * fgb

                     - 0.75 * pa_x * fz * fga * fx * fx

                     - 0.5 * pa_xxx * fx * fz * fgb

                     + 3.0 * pa_x * fx * fx * fx * fz

                     + 2.5 * pa_xxx * fz * fx * fx

                     - 1.5 * pa_xz * fx * fz * fgb * pb_z

                     - 1.5 * pa_xz * fz * fga * fx * pb_z

                     - 1.5 * pa_x * fz * fga * fx * pb_yy

                     - pa_xxxz * fz * fgb * pb_z

                     + 7.5 * pa_xz * fx * fx * fz * pb_z

                     + 7.5 * pa_x * fx * fx * fz * pb_yy

                     + 6.0 * pa_xxxz * fz * fx * pb_z

                     + 6.0 * pa_xxx * fz * fx * pb_yy

                     - 3.0 * pa_xz * fz * fga * pb_yyz

                     + 18.0 * pa_xz * fx * fz * pb_yyz

                     + 14.0 * pa_xxxz * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_yzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xxx,
                                    double pa_xxxz,
                                    double pa_xz,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xz * fx * fx * pb_y

                     + 1.5 * pa_x * fx * fx * pb_yz

                     + 0.5 * pa_xxxz * pb_y * fx

                     + pa_xxx * fx * pb_yz

                     + 1.5 * pa_xz * fx * pb_yzz

                     + pa_xxxz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_yzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xxx,
                                    double pa_xxxz,
                                    double pa_xz,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xz * fx * pb_y * fz * fgb

                     - 1.5 * pa_xz * fz * fga * pb_y * fx

                     - 3.0 * pa_x * fz * fga * fx * pb_yz

                     - pa_xxxz * pb_y * fz * fgb

                     + 7.5 * pa_xz * fx * fx * fz * pb_y

                     + 15.0 * pa_x * fx * fx * fz * pb_yz

                     + 6.0 * pa_xxxz * fz * pb_y * fx

                     + 12.0 * pa_xxx * fz * fx * pb_yz

                     - 3.0 * pa_xz * fz * fga * pb_yzz

                     + 18.0 * pa_xz * fx * fz * pb_yzz

                     + 14.0 * pa_xxxz * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_zzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xxx,
                                    double pa_xxxz,
                                    double pa_xz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * pa_x * fx * fx * fx

                     + 0.75 * pa_xxx * fx * fx

                     + 2.25 * pa_xz * fx * fx * pb_z

                     + 2.25 * pa_x * fx * fx * pb_zz

                     + 1.5 * pa_xxxz * pb_z * fx

                     + 1.5 * pa_xxx * fx * pb_zz

                     + 1.5 * pa_xz * fx * pb_zzz

                     + pa_xxxz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_zzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xxx,
                                    double pa_xxxz,
                                    double pa_xz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 2.25 * pa_x * fx * fx * fz * fgb

                     - 2.25 * pa_x * fz * fga * fx * fx

                     - 1.5 * pa_xxx * fx * fz * fgb

                     + 9.0 * pa_x * fx * fx * fx * fz

                     + 7.5 * pa_xxx * fz * fx * fx

                     - 4.5 * pa_xz * fx * pb_z * fz * fgb

                     - 4.5 * pa_xz * fz * fga * pb_z * fx

                     - 4.5 * pa_x * fz * fga * fx * pb_zz

                     - 3.0 * pa_xxxz * pb_z * fz * fgb

                     + 22.5 * pa_xz * fx * fx * fz * pb_z

                     + 22.5 * pa_x * fx * fx * fz * pb_zz

                     + 18.0 * pa_xxxz * fz * pb_z * fx

                     + 18.0 * pa_xxx * fz * fx * pb_zz

                     - 3.0 * pa_xz * fz * fga * pb_zzz

                     + 18.0 * pa_xz * fx * fz * pb_zzz

                     + 14.0 * pa_xxxz * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xxx_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxyy,
                                    double pa_xyy,
                                    double pa_yy,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * fx

                     + 1.125 * fx * fx * fx * pb_x

                     + 1.5 * pa_xyy * fx * fx

                     + 2.25 * fx * fx * pa_yy * pb_x

                     + 0.75 * pa_xx * fx * fx * pb_x

                     + 1.5 * pa_x * fx * fx * pb_xx

                     + 1.5 * pa_xxyy * pb_x * fx

                     + 3.0 * pa_xyy * fx * pb_xx

                     + 0.25 * fx * fx * pb_xxx

                     + 0.5 * pa_xx * fx * pb_xxx

                     + 0.5 * fx * pa_yy * pb_xxx

                     + pa_xxyy * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxyy,
                                    double pa_xyy,
                                    double pa_yy,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fx * fz * fgb

                     - 1.5 * pa_x * fx * fx * fz * fga

                     - 3.0 * fx * fx * fz * fga * pb_x

                     - 3.0 * pa_xyy * fx * fz * fgb

                     + 6.0 * pa_x * fx * fx * fx * fz

                     + 9.0 * fx * fx * fx * fz * pb_x

                     + 15.0 * pa_xyy * fx * fx * fz

                     + 22.5 * fx * fx * pa_yy * fz * pb_x

                     - 0.75 * fx * fx * pb_x * fz * fgb

                     - 1.5 * pa_xx * fx * pb_x * fz * fgb

                     - 1.5 * pa_xx * fz * fga * pb_x * fx

                     - 3.0 * pa_x * fx * fz * fga * pb_xx

                     - 1.5 * fx * pa_yy * pb_x * fz * fgb

                     - 1.5 * fz * fga * pa_yy * pb_x * fx

                     - 3.0 * pa_xxyy * pb_x * fz * fgb

                     + 7.5 * pa_xx * fz * fx * fx * pb_x

                     + 15.0 * pa_x * fx * fx * fz * pb_xx

                     + 18.0 * pa_xxyy * fz * pb_x * fx

                     + 36.0 * pa_xyy * fx * fz * pb_xx

                     - fx * fz * fga * pb_xxx

                     - pa_xx * fz * fga * pb_xxx

                     + 2.5 * fx * fx * fz * pb_xxx

                     - fz * fga * pa_yy * pb_xxx

                     + 6.0 * pa_xx * fz * fx * pb_xxx

                     + 6.0 * fx * pa_yy * fz * pb_xxx

                     + 14.0 * pa_xxyy * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xxy_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_xxyy,
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
        return s_0_0 * (0.75 * fx * fx * fx * pa_y

                     + 0.375 * fx * fx * fx * pb_y

                     + 0.5 * pa_xxy * fx * fx

                     + 2.0 * pa_xy * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_yy * pb_y

                     + 0.25 * pa_xx * fx * fx * pb_y

                     + pa_x * fx * fx * pb_xy

                     + 0.5 * fx * fx * pa_y * pb_xx

                     + 0.5 * pa_xxyy * fx * pb_y

                     + pa_xxy * fx * pb_xx

                     + 2.0 * pa_xyy * fx * pb_xy

                     + 0.25 * fx * fx * pb_xxy

                     + 0.5 * pa_xx * fx * pb_xxy

                     + 0.5 * fx * pa_yy * pb_xxy

                     + pa_xxyy * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_xxyy,
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
        return r_0_0 * (6.0 * fx * fx * fx * pa_y * fz

                     - 0.5 * fx * fx * pa_y * fz * fgb

                     - fx * fx * fz * fga * pb_y

                     - 0.5 * fz * fga * pa_y * fx * fx

                     - pa_xxy * fx * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_y

                     + 5.0 * pa_xxy * fz * fx * fx

                     + 20.0 * pa_xy * fx * fx * fz * pb_x

                     + 7.5 * fx * fx * pa_yy * fz * pb_y

                     - 0.25 * fx * fx * fz * fgb * pb_y

                     - 0.5 * pa_xx * fx * fz * fgb * pb_y

                     - 0.5 * pa_xx * fz * fga * fx * pb_y

                     - 2.0 * pa_x * fx * fz * fga * pb_xy

                     - 0.5 * fx * pa_yy * fz * fgb * pb_y

                     - 0.5 * fz * fga * pa_yy * fx * pb_y

                     - fz * fga * pa_y * fx * pb_xx

                     - pa_xxyy * fz * fgb * pb_y

                     + 2.5 * pa_xx * fz * fx * fx * pb_y

                     + 10.0 * pa_x * fx * fx * fz * pb_xy

                     + 5.0 * fx * fx * pa_y * fz * pb_xx

                     + 6.0 * pa_xxyy * fz * fx * pb_y

                     + 12.0 * pa_xxy * fz * fx * pb_xx

                     + 24.0 * pa_xyy * fx * fz * pb_xy

                     - fx * fz * fga * pb_xxy

                     - pa_xx * fz * fga * pb_xxy

                     + 2.5 * fx * fx * fz * pb_xxy

                     - fz * fga * pa_yy * pb_xxy

                     + 6.0 * pa_xx * fz * fx * pb_xxy

                     + 6.0 * fx * pa_yy * fz * pb_xxy

                     + 14.0 * pa_xxyy * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xxz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxyy,
                                    double pa_xyy,
                                    double pa_yy,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_z

                     + 0.75 * fx * fx * pa_yy * pb_z

                     + 0.25 * pa_xx * fx * fx * pb_z

                     + pa_x * fx * fx * pb_xz

                     + 0.5 * pa_xxyy * fx * pb_z

                     + 2.0 * pa_xyy * fx * pb_xz

                     + 0.25 * fx * fx * pb_xxz

                     + 0.5 * pa_xx * fx * pb_xxz

                     + 0.5 * fx * pa_yy * pb_xxz

                     + pa_xxyy * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxyy,
                                    double pa_xyy,
                                    double pa_yy,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- fx * fx * fz * fga * pb_z

                     + 3.0 * fx * fx * fx * fz * pb_z

                     + 7.5 * fx * fx * pa_yy * fz * pb_z

                     - 0.25 * fx * fx * fz * fgb * pb_z

                     - 0.5 * pa_xx * fx * fz * fgb * pb_z

                     - 0.5 * pa_xx * fz * fga * fx * pb_z

                     - 2.0 * pa_x * fx * fz * fga * pb_xz

                     - 0.5 * fx * pa_yy * fz * fgb * pb_z

                     - 0.5 * fz * fga * pa_yy * fx * pb_z

                     - pa_xxyy * fz * fgb * pb_z

                     + 2.5 * pa_xx * fz * fx * fx * pb_z

                     + 10.0 * pa_x * fx * fx * fz * pb_xz

                     + 6.0 * pa_xxyy * fz * fx * pb_z

                     + 24.0 * pa_xyy * fx * fz * pb_xz

                     - fx * fz * fga * pb_xxz

                     - pa_xx * fz * fga * pb_xxz

                     + 2.5 * fx * fx * fz * pb_xxz

                     - fz * fga * pa_yy * pb_xxz

                     + 6.0 * pa_xx * fz * fx * pb_xxz

                     + 6.0 * fx * pa_yy * fz * pb_xxz

                     + 14.0 * pa_xxyy * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xyy_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_xxyy,
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
        return s_0_0 * (0.75 * pa_x * fx * fx * fx

                     + 0.375 * fx * fx * fx * pb_x

                     + 0.75 * pa_xx * fx * fx * pb_x

                     + 0.5 * pa_xyy * fx * fx

                     + 2.0 * pa_xy * fx * fx * pb_y

                     + 0.5 * pa_x * fx * fx * pb_yy

                     + 0.25 * fx * fx * pa_yy * pb_x

                     + fx * fx * pa_y * pb_xy

                     + 0.5 * pa_xxyy * pb_x * fx

                     + 2.0 * pa_xxy * fx * pb_xy

                     + pa_xyy * fx * pb_yy

                     + 0.25 * fx * fx * pb_xyy

                     + 0.5 * pa_xx * fx * pb_xyy

                     + 0.5 * fx * pa_yy * pb_xyy

                     + pa_xxyy * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_xxyy,
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
        return r_0_0 * (6.0 * pa_x * fx * fx * fx * fz

                     - 0.5 * pa_x * fx * fx * fz * fgb

                     - 0.5 * pa_x * fx * fx * fz * fga

                     - fz * fga * fx * fx * pb_x

                     - pa_xyy * fx * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_x

                     + 7.5 * pa_xx * fx * fx * fz * pb_x

                     + 5.0 * pa_xyy * fx * fx * fz

                     + 20.0 * pa_xy * fx * fx * fz * pb_y

                     - 0.25 * fx * fx * pb_x * fz * fgb

                     - 0.5 * pa_xx * fx * pb_x * fz * fgb

                     - 0.5 * pa_xx * fz * fga * pb_x * fx

                     - pa_x * fx * fz * fga * pb_yy

                     - 0.5 * fx * pa_yy * pb_x * fz * fgb

                     - 0.5 * fz * fga * pa_yy * pb_x * fx

                     - 2.0 * fz * fga * pa_y * fx * pb_xy

                     - pa_xxyy * pb_x * fz * fgb

                     + 5.0 * pa_x * fx * fx * fz * pb_yy

                     + 2.5 * fx * fx * pa_yy * fz * pb_x

                     + 10.0 * fx * fx * pa_y * fz * pb_xy

                     + 6.0 * pa_xxyy * fz * pb_x * fx

                     + 24.0 * pa_xxy * fz * fx * pb_xy

                     + 12.0 * pa_xyy * fx * fz * pb_yy

                     - fx * fz * fga * pb_xyy

                     - pa_xx * fz * fga * pb_xyy

                     + 2.5 * fx * fx * fz * pb_xyy

                     - fz * fga * pa_yy * pb_xyy

                     + 6.0 * pa_xx * fz * fx * pb_xyy

                     + 6.0 * fx * pa_yy * fz * pb_xyy

                     + 14.0 * pa_xxyy * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_xxyy,
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
        return s_0_0 * (pa_xy * fx * fx * pb_z

                     + 0.5 * pa_x * fx * fx * pb_yz

                     + 0.5 * fx * fx * pa_y * pb_xz

                     + pa_xxy * fx * pb_xz

                     + pa_xyy * fx * pb_yz

                     + 0.25 * fx * fx * pb_xyz

                     + 0.5 * pa_xx * fx * pb_xyz

                     + 0.5 * fx * pa_yy * pb_xyz

                     + pa_xxyy * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xyz_r_0(double fga,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_xxyy,
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
        return r_0_0 * (10.0 * pa_xy * fx * fx * fz * pb_z

                     - pa_x * fx * fz * fga * pb_yz

                     - fz * fga * pa_y * fx * pb_xz

                     + 5.0 * pa_x * fx * fx * fz * pb_yz

                     + 5.0 * fx * fx * pa_y * fz * pb_xz

                     + 12.0 * pa_xxy * fz * fx * pb_xz

                     + 12.0 * pa_xyy * fx * fz * pb_yz

                     - fx * fz * fga * pb_xyz

                     - pa_xx * fz * fga * pb_xyz

                     + 2.5 * fx * fx * fz * pb_xyz

                     - fz * fga * pa_yy * pb_xyz

                     + 6.0 * pa_xx * fz * fx * pb_xyz

                     + 6.0 * fx * pa_yy * fz * pb_xyz

                     + 14.0 * pa_xxyy * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxyy,
                                    double pa_xyy,
                                    double pa_yy,
                                    double pb_x,
                                    double pb_xzz,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (0.25 * pa_x * fx * fx * fx

                     + 0.5 * pa_xyy * fx * fx

                     + 0.125 * fx * fx * fx * pb_x

                     + 0.25 * pa_xx * fx * fx * pb_x

                     + 0.5 * pa_x * fx * fx * pb_zz

                     + 0.25 * fx * fx * pa_yy * pb_x

                     + 0.5 * pa_xxyy * pb_x * fx

                     + pa_xyy * fx * pb_zz

                     + 0.25 * fx * fx * pb_xzz

                     + 0.5 * pa_xx * fx * pb_xzz

                     + 0.5 * fx * pa_yy * pb_xzz

                     + pa_xxyy * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxyy,
                                    double pa_xyy,
                                    double pa_yy,
                                    double pb_x,
                                    double pb_xzz,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_x * fx * fx * fz * fgb

                     - 0.5 * pa_x * fx * fx * fz * fga

                     - pa_xyy * fx * fz * fgb

                     + 2.0 * pa_x * fx * fx * fx * fz

                     + 5.0 * pa_xyy * fx * fx * fz

                     - 0.25 * fx * fx * pb_x * fz * fgb

                     - 0.5 * fx * fx * fz * fga * pb_x

                     - 0.5 * pa_xx * fx * pb_x * fz * fgb

                     - 0.5 * pa_xx * fz * fga * pb_x * fx

                     - pa_x * fx * fz * fga * pb_zz

                     - 0.5 * fx * pa_yy * pb_x * fz * fgb

                     + fx * fx * fx * fz * pb_x

                     - 0.5 * fz * fga * pa_yy * pb_x * fx

                     - pa_xxyy * pb_x * fz * fgb

                     + 2.5 * pa_xx * fz * fx * fx * pb_x

                     + 5.0 * pa_x * fx * fx * fz * pb_zz

                     + 2.5 * fx * fx * pa_yy * fz * pb_x

                     + 6.0 * pa_xxyy * fz * pb_x * fx

                     + 12.0 * pa_xyy * fx * fz * pb_zz

                     - fx * fz * fga * pb_xzz

                     - pa_xx * fz * fga * pb_xzz

                     + 2.5 * fx * fx * fz * pb_xzz

                     - fz * fga * pa_yy * pb_xzz

                     + 6.0 * pa_xx * fz * fx * pb_xzz

                     + 6.0 * fx * pa_yy * fz * pb_xzz

                     + 14.0 * pa_xxyy * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_yyy_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_xxyy,
                                    double pa_y,
                                    double pa_yy,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * fx * pa_y

                     + 1.125 * fx * fx * fx * pb_y

                     + 1.5 * pa_xxy * fx * fx

                     + 2.25 * pa_xx * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_yy * pb_y

                     + 1.5 * fx * fx * pa_y * pb_yy

                     + 1.5 * pa_xxyy * pb_y * fx

                     + 3.0 * pa_xxy * fx * pb_yy

                     + 0.25 * fx * fx * pb_yyy

                     + 0.5 * pa_xx * fx * pb_yyy

                     + 0.5 * fx * pa_yy * pb_yyy

                     + pa_xxyy * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_yyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_xxyy,
                                    double pa_y,
                                    double pa_yy,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_y * fz * fgb

                     - 1.5 * fz * fga * pa_y * fx * fx

                     - 3.0 * fz * fga * fx * fx * pb_y

                     - 3.0 * pa_xxy * fx * fz * fgb

                     + 6.0 * fx * fx * fx * pa_y * fz

                     + 9.0 * fx * fx * fx * fz * pb_y

                     + 15.0 * pa_xxy * fz * fx * fx

                     + 22.5 * pa_xx * fx * fx * fz * pb_y

                     - 0.75 * fx * fx * pb_y * fz * fgb

                     - 1.5 * pa_xx * fx * pb_y * fz * fgb

                     - 1.5 * pa_xx * fz * fga * pb_y * fx

                     - 1.5 * fx * pa_yy * pb_y * fz * fgb

                     - 1.5 * fz * fga * pa_yy * pb_y * fx

                     - 3.0 * fz * fga * pa_y * fx * pb_yy

                     - 3.0 * pa_xxyy * pb_y * fz * fgb

                     + 7.5 * fx * fx * pa_yy * fz * pb_y

                     + 15.0 * fx * fx * pa_y * fz * pb_yy

                     + 18.0 * pa_xxyy * fz * pb_y * fx

                     + 36.0 * pa_xxy * fz * fx * pb_yy

                     - fx * fz * fga * pb_yyy

                     - pa_xx * fz * fga * pb_yyy

                     + 2.5 * fx * fx * fz * pb_yyy

                     - fz * fga * pa_yy * pb_yyy

                     + 6.0 * pa_xx * fz * fx * pb_yyy

                     + 6.0 * fx * pa_yy * fz * pb_yyy

                     + 14.0 * pa_xxyy * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_yyz_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_xxyy,
                                    double pa_y,
                                    double pa_yy,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_z

                     + 0.75 * pa_xx * fx * fx * pb_z

                     + 0.25 * fx * fx * pa_yy * pb_z

                     + fx * fx * pa_y * pb_yz

                     + 0.5 * pa_xxyy * fx * pb_z

                     + 2.0 * pa_xxy * fx * pb_yz

                     + 0.25 * fx * fx * pb_yyz

                     + 0.5 * pa_xx * fx * pb_yyz

                     + 0.5 * fx * pa_yy * pb_yyz

                     + pa_xxyy * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_yyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_xxyy,
                                    double pa_y,
                                    double pa_yy,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- fz * fga * fx * fx * pb_z

                     + 3.0 * fx * fx * fx * fz * pb_z

                     + 7.5 * pa_xx * fx * fx * fz * pb_z

                     - 0.25 * fx * fx * fz * fgb * pb_z

                     - 0.5 * pa_xx * fx * fz * fgb * pb_z

                     - 0.5 * pa_xx * fz * fga * fx * pb_z

                     - 0.5 * fx * pa_yy * fz * fgb * pb_z

                     - 0.5 * fz * fga * pa_yy * fx * pb_z

                     - 2.0 * fz * fga * pa_y * fx * pb_yz

                     - pa_xxyy * fz * fgb * pb_z

                     + 2.5 * fx * fx * pa_yy * fz * pb_z

                     + 10.0 * fx * fx * pa_y * fz * pb_yz

                     + 6.0 * pa_xxyy * fz * fx * pb_z

                     + 24.0 * pa_xxy * fz * fx * pb_yz

                     - fx * fz * fga * pb_yyz

                     - pa_xx * fz * fga * pb_yyz

                     + 2.5 * fx * fx * fz * pb_yyz

                     - fz * fga * pa_yy * pb_yyz

                     + 6.0 * pa_xx * fz * fx * pb_yyz

                     + 6.0 * fx * pa_yy * fz * pb_yyz

                     + 14.0 * pa_xxyy * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_yzz_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_xxyy,
                                    double pa_y,
                                    double pa_yy,
                                    double pb_y,
                                    double pb_yzz,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * fx * pa_y

                     + 0.5 * pa_xxy * fx * fx

                     + 0.125 * fx * fx * fx * pb_y

                     + 0.25 * pa_xx * fx * fx * pb_y

                     + 0.25 * fx * fx * pa_yy * pb_y

                     + 0.5 * fx * fx * pa_y * pb_zz

                     + 0.5 * pa_xxyy * pb_y * fx

                     + pa_xxy * fx * pb_zz

                     + 0.25 * fx * fx * pb_yzz

                     + 0.5 * pa_xx * fx * pb_yzz

                     + 0.5 * fx * pa_yy * pb_yzz

                     + pa_xxyy * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_yzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_xxyy,
                                    double pa_y,
                                    double pa_yy,
                                    double pb_y,
                                    double pb_yzz,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * pa_y * fz * fgb

                     - 0.5 * fz * fga * pa_y * fx * fx

                     - pa_xxy * fx * fz * fgb

                     + 2.0 * fx * fx * fx * pa_y * fz

                     + 5.0 * pa_xxy * fz * fx * fx

                     - 0.25 * fx * fx * pb_y * fz * fgb

                     - 0.5 * fx * fx * fz * fga * pb_y

                     - 0.5 * pa_xx * fx * pb_y * fz * fgb

                     - 0.5 * pa_xx * fz * fga * pb_y * fx

                     - 0.5 * fx * pa_yy * pb_y * fz * fgb

                     + fx * fx * fx * fz * pb_y

                     - 0.5 * fz * fga * pa_yy * pb_y * fx

                     - fz * fga * pa_y * fx * pb_zz

                     - pa_xxyy * pb_y * fz * fgb

                     + 2.5 * pa_xx * fz * fx * fx * pb_y

                     + 2.5 * fx * fx * pa_yy * fz * pb_y

                     + 5.0 * fx * fx * pa_y * fz * pb_zz

                     + 6.0 * pa_xxyy * fz * pb_y * fx

                     + 12.0 * pa_xxy * fz * fx * pb_zz

                     - fx * fz * fga * pb_yzz

                     - pa_xx * fz * fga * pb_yzz

                     + 2.5 * fx * fx * fz * pb_yzz

                     - fz * fga * pa_yy * pb_yzz

                     + 6.0 * pa_xx * fz * fx * pb_yzz

                     + 6.0 * fx * pa_yy * fz * pb_yzz

                     + 14.0 * pa_xxyy * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_zzz_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxyy,
                                    double pa_yy,
                                    double pb_z,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_z

                     + 0.75 * pa_xx * fx * fx * pb_z

                     + 0.75 * fx * fx * pa_yy * pb_z

                     + 1.5 * pa_xxyy * pb_z * fx

                     + 0.25 * fx * fx * pb_zzz

                     + 0.5 * pa_xx * fx * pb_zzz

                     + 0.5 * fx * pa_yy * pb_zzz

                     + pa_xxyy * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_zzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxyy,
                                    double pa_yy,
                                    double pb_z,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pb_z * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_z

                     - 1.5 * pa_xx * fx * pb_z * fz * fgb

                     - 1.5 * pa_xx * fz * fga * pb_z * fx

                     - 1.5 * fx * pa_yy * pb_z * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_z

                     - 1.5 * fz * fga * pa_yy * pb_z * fx

                     - 3.0 * pa_xxyy * pb_z * fz * fgb

                     + 7.5 * pa_xx * fz * fx * fx * pb_z

                     + 7.5 * fx * fx * pa_yy * fz * pb_z

                     + 18.0 * pa_xxyy * fz * pb_z * fx

                     - fx * fz * fga * pb_zzz

                     - pa_xx * fz * fga * pb_zzz

                     + 2.5 * fx * fx * fz * pb_zzz

                     - fz * fga * pa_yy * pb_zzz

                     + 6.0 * pa_xx * fz * fx * pb_zzz

                     + 6.0 * fx * pa_yy * fz * pb_zzz

                     + 14.0 * pa_xxyy * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xxx_s_0(double fx,
                                    double pa_xxyz,
                                    double pa_xyz,
                                    double pa_yz,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xyz * fx * fx

                     + 2.25 * fx * fx * pa_yz * pb_x

                     + 1.5 * pa_xxyz * pb_x * fx

                     + 3.0 * pa_xyz * fx * pb_xx

                     + 0.5 * fx * pa_yz * pb_xxx

                     + pa_xxyz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxyz,
                                    double pa_xyz,
                                    double pa_yz,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xyz * fx * fz * fgb

                     + 15.0 * pa_xyz * fx * fx * fz

                     + 22.5 * fx * fx * pa_yz * fz * pb_x

                     - 1.5 * fx * pa_yz * pb_x * fz * fgb

                     - 1.5 * fz * fga * pa_yz * pb_x * fx

                     - 3.0 * pa_xxyz * pb_x * fz * fgb

                     + 18.0 * pa_xxyz * fz * pb_x * fx

                     + 36.0 * pa_xyz * fx * fz * pb_xx

                     - fz * fga * pa_yz * pb_xxx

                     + 6.0 * fx * pa_yz * fz * pb_xxx

                     + 14.0 * pa_xxyz * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xxy_s_0(double fx,
                                    double pa_xxyz,
                                    double pa_xxz,
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
        return s_0_0 * (0.375 * fx * fx * fx * pa_z

                     + 0.25 * pa_xxz * fx * fx

                     + pa_xz * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_yz * pb_y

                     + 0.25 * fx * fx * pa_z * pb_xx

                     + 0.5 * pa_xxyz * fx * pb_y

                     + 0.5 * pa_xxz * fx * pb_xx

                     + 2.0 * pa_xyz * fx * pb_xy

                     + 0.5 * fx * pa_yz * pb_xxy

                     + pa_xxyz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxyz,
                                    double pa_xxz,
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
        return r_0_0 * (3.0 * fx * fx * fx * fz * pa_z

                     - 0.25 * fx * fx * pa_z * fz * fgb

                     - 0.25 * fz * fga * fx * fx * pa_z

                     - 0.5 * pa_xxz * fx * fz * fgb

                     + 2.5 * pa_xxz * fx * fx * fz

                     + 10.0 * pa_xz * fx * fx * fz * pb_x

                     + 7.5 * fx * fx * pa_yz * fz * pb_y

                     - 0.5 * fx * pa_yz * fz * fgb * pb_y

                     - 0.5 * fz * fga * pa_yz * fx * pb_y

                     - 0.5 * fz * fga * fx * pa_z * pb_xx

                     - pa_xxyz * fz * fgb * pb_y

                     + 2.5 * fx * fx * fz * pa_z * pb_xx

                     + 6.0 * pa_xxyz * fz * fx * pb_y

                     + 6.0 * pa_xxz * fx * fz * pb_xx

                     + 24.0 * pa_xyz * fx * fz * pb_xy

                     - fz * fga * pa_yz * pb_xxy

                     + 6.0 * fx * pa_yz * fz * pb_xxy

                     + 14.0 * pa_xxyz * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xxz_s_0(double fx,
                                    double pa_xxy,
                                    double pa_xxyz,
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
        return s_0_0 * (0.375 * fx * fx * fx * pa_y

                     + 0.25 * pa_xxy * fx * fx

                     + pa_xy * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_yz * pb_z

                     + 0.25 * fx * fx * pa_y * pb_xx

                     + 0.5 * pa_xxyz * fx * pb_z

                     + 0.5 * pa_xxy * fx * pb_xx

                     + 2.0 * pa_xyz * fx * pb_xz

                     + 0.5 * fx * pa_yz * pb_xxz

                     + pa_xxyz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxy,
                                    double pa_xxyz,
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
        return r_0_0 * (3.0 * fx * fx * fx * pa_y * fz

                     - 0.25 * fx * fx * pa_y * fz * fgb

                     - 0.25 * fz * fga * pa_y * fx * fx

                     - 0.5 * pa_xxy * fx * fz * fgb

                     + 2.5 * pa_xxy * fz * fx * fx

                     + 10.0 * pa_xy * fx * fx * fz * pb_x

                     + 7.5 * fx * fx * pa_yz * fz * pb_z

                     - 0.5 * fx * pa_yz * fz * fgb * pb_z

                     - 0.5 * fz * fga * pa_yz * fx * pb_z

                     - 0.5 * fz * fga * pa_y * fx * pb_xx

                     - pa_xxyz * fz * fgb * pb_z

                     + 2.5 * fx * fx * pa_y * fz * pb_xx

                     + 6.0 * pa_xxyz * fz * fx * pb_z

                     + 6.0 * pa_xxy * fz * fx * pb_xx

                     + 24.0 * pa_xyz * fx * fz * pb_xz

                     - fz * fga * pa_yz * pb_xxz

                     + 6.0 * fx * pa_yz * fz * pb_xxz

                     + 14.0 * pa_xxyz * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xyy_s_0(double fx,
                                    double pa_xxyz,
                                    double pa_xxz,
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
        return s_0_0 * (0.5 * pa_xyz * fx * fx

                     + pa_xz * fx * fx * pb_y

                     + 0.25 * fx * fx * pa_yz * pb_x

                     + 0.5 * fx * fx * pa_z * pb_xy

                     + 0.5 * pa_xxyz * pb_x * fx

                     + pa_xxz * fx * pb_xy

                     + pa_xyz * fx * pb_yy

                     + 0.5 * fx * pa_yz * pb_xyy

                     + pa_xxyz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxyz,
                                    double pa_xxz,
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
        return r_0_0 * (- pa_xyz * fx * fz * fgb

                     + 5.0 * pa_xyz * fx * fx * fz

                     + 10.0 * pa_xz * fx * fx * fz * pb_y

                     - 0.5 * fx * pa_yz * pb_x * fz * fgb

                     - 0.5 * fz * fga * pa_yz * pb_x * fx

                     - fz * fga * fx * pa_z * pb_xy

                     - pa_xxyz * pb_x * fz * fgb

                     + 2.5 * fx * fx * pa_yz * fz * pb_x

                     + 5.0 * fx * fx * fz * pa_z * pb_xy

                     + 6.0 * pa_xxyz * fz * pb_x * fx

                     + 12.0 * pa_xxz * fx * fz * pb_xy

                     + 12.0 * pa_xyz * fx * fz * pb_yy

                     - fz * fga * pa_yz * pb_xyy

                     + 6.0 * fx * pa_yz * fz * pb_xyy

                     + 14.0 * pa_xxyz * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_xxyz,
                                    double pa_xxz,
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
        return s_0_0 * (0.25 * pa_x * fx * fx * fx

                     + 0.125 * fx * fx * fx * pb_x

                     + 0.25 * pa_xx * fx * fx * pb_x

                     + 0.5 * pa_xy * fx * fx * pb_y

                     + 0.5 * pa_xz * fx * fx * pb_z

                     + 0.25 * fx * fx * pa_y * pb_xy

                     + 0.25 * fx * fx * pa_z * pb_xz

                     + 0.5 * pa_xxy * fx * pb_xy

                     + 0.5 * pa_xxz * fx * pb_xz

                     + pa_xyz * fx * pb_yz

                     + 0.5 * fx * pa_yz * pb_xyz

                     + pa_xxyz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xyz_r_0(double fga,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_xxyz,
                                    double pa_xxz,
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
        return r_0_0 * (2.0 * pa_x * fx * fx * fx * fz

                     - 0.25 * fz * fga * fx * fx * pb_x

                     + fx * fx * fx * fz * pb_x

                     + 2.5 * pa_xx * fx * fx * fz * pb_x

                     + 5.0 * pa_xy * fx * fx * fz * pb_y

                     + 5.0 * pa_xz * fx * fx * fz * pb_z

                     - 0.5 * fz * fga * pa_y * fx * pb_xy

                     - 0.5 * fz * fga * fx * pa_z * pb_xz

                     + 2.5 * fx * fx * pa_y * fz * pb_xy

                     + 2.5 * fx * fx * fz * pa_z * pb_xz

                     + 6.0 * pa_xxy * fz * fx * pb_xy

                     + 6.0 * pa_xxz * fx * fz * pb_xz

                     + 12.0 * pa_xyz * fx * fz * pb_yz

                     - fz * fga * pa_yz * pb_xyz

                     + 6.0 * fx * pa_yz * fz * pb_xyz

                     + 14.0 * pa_xxyz * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xzz_s_0(double fx,
                                    double pa_xxy,
                                    double pa_xxyz,
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
        return s_0_0 * (0.5 * pa_xyz * fx * fx

                     + pa_xy * fx * fx * pb_z

                     + 0.25 * fx * fx * pa_yz * pb_x

                     + 0.5 * fx * fx * pa_y * pb_xz

                     + 0.5 * pa_xxyz * pb_x * fx

                     + pa_xxy * fx * pb_xz

                     + pa_xyz * fx * pb_zz

                     + 0.5 * fx * pa_yz * pb_xzz

                     + pa_xxyz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxy,
                                    double pa_xxyz,
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
        return r_0_0 * (- pa_xyz * fx * fz * fgb

                     + 5.0 * pa_xyz * fx * fx * fz

                     + 10.0 * pa_xy * fx * fx * fz * pb_z

                     - 0.5 * fx * pa_yz * pb_x * fz * fgb

                     - 0.5 * fz * fga * pa_yz * pb_x * fx

                     - fz * fga * pa_y * fx * pb_xz

                     - pa_xxyz * pb_x * fz * fgb

                     + 2.5 * fx * fx * pa_yz * fz * pb_x

                     + 5.0 * fx * fx * pa_y * fz * pb_xz

                     + 6.0 * pa_xxyz * fz * pb_x * fx

                     + 12.0 * pa_xxy * fz * fx * pb_xz

                     + 12.0 * pa_xyz * fx * fz * pb_zz

                     - fz * fga * pa_yz * pb_xzz

                     + 6.0 * fx * pa_yz * fz * pb_xzz

                     + 14.0 * pa_xxyz * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_yyy_s_0(double fx,
                                    double pa_xxyz,
                                    double pa_xxz,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_z

                     + 0.75 * pa_xxz * fx * fx

                     + 0.75 * fx * fx * pa_yz * pb_y

                     + 0.75 * fx * fx * pa_z * pb_yy

                     + 1.5 * pa_xxyz * pb_y * fx

                     + 1.5 * pa_xxz * fx * pb_yy

                     + 0.5 * fx * pa_yz * pb_yyy

                     + pa_xxyz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_yyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxyz,
                                    double pa_xxz,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pa_z * fz * fgb

                     - 0.75 * fz * fga * fx * fx * pa_z

                     - 1.5 * pa_xxz * fx * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_z

                     + 7.5 * pa_xxz * fx * fx * fz

                     - 1.5 * fx * pa_yz * pb_y * fz * fgb

                     - 1.5 * fz * fga * pa_yz * pb_y * fx

                     - 1.5 * fz * fga * fx * pa_z * pb_yy

                     - 3.0 * pa_xxyz * pb_y * fz * fgb

                     + 7.5 * fx * fx * pa_yz * fz * pb_y

                     + 7.5 * fx * fx * fz * pa_z * pb_yy

                     + 18.0 * pa_xxyz * fz * pb_y * fx

                     + 18.0 * pa_xxz * fx * fz * pb_yy

                     - fz * fga * pa_yz * pb_yyy

                     + 6.0 * fx * pa_yz * fz * pb_yyy

                     + 14.0 * pa_xxyz * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_yyz_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_xxyz,
                                    double pa_xxz,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx * pa_y

                     + 0.25 * fx * fx * fx * pb_y

                     + 0.25 * pa_xxy * fx * fx

                     + 0.5 * pa_xx * fx * fx * pb_y

                     + 0.25 * fx * fx * pa_yz * pb_z

                     + 0.25 * fx * fx * pa_y * pb_yy

                     + 0.5 * fx * fx * pa_z * pb_yz

                     + 0.5 * pa_xxyz * fx * pb_z

                     + 0.5 * pa_xxy * fx * pb_yy

                     + pa_xxz * fx * pb_yz

                     + 0.5 * fx * pa_yz * pb_yyz

                     + pa_xxyz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_yyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_xxyz,
                                    double pa_xxz,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.25 * fx * fx * pa_y * fz * fgb

                     - 0.25 * fz * fga * pa_y * fx * fx

                     - 0.5 * fz * fga * fx * fx * pb_y

                     - 0.5 * pa_xxy * fx * fz * fgb

                     + fx * fx * fx * pa_y * fz

                     + 2.0 * fx * fx * fx * fz * pb_y

                     + 2.5 * pa_xxy * fz * fx * fx

                     + 5.0 * pa_xx * fx * fx * fz * pb_y

                     - 0.5 * fx * pa_yz * fz * fgb * pb_z

                     - 0.5 * fz * fga * pa_yz * fx * pb_z

                     - 0.5 * fz * fga * pa_y * fx * pb_yy

                     - fz * fga * fx * pa_z * pb_yz

                     - pa_xxyz * fz * fgb * pb_z

                     + 2.5 * fx * fx * pa_yz * fz * pb_z

                     + 2.5 * fx * fx * pa_y * fz * pb_yy

                     + 5.0 * fx * fx * fz * pa_z * pb_yz

                     + 6.0 * pa_xxyz * fz * fx * pb_z

                     + 6.0 * pa_xxy * fz * fx * pb_yy

                     + 12.0 * pa_xxz * fx * fz * pb_yz

                     - fz * fga * pa_yz * pb_yyz

                     + 6.0 * fx * pa_yz * fz * pb_yyz

                     + 14.0 * pa_xxyz * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_yzz_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_xxyz,
                                    double pa_xxz,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double pb_z,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx * pa_z

                     + 0.25 * fx * fx * fx * pb_z

                     + 0.25 * pa_xxz * fx * fx

                     + 0.5 * pa_xx * fx * fx * pb_z

                     + 0.25 * fx * fx * pa_yz * pb_y

                     + 0.5 * fx * fx * pa_y * pb_yz

                     + 0.25 * fx * fx * pa_z * pb_zz

                     + 0.5 * pa_xxyz * pb_y * fx

                     + pa_xxy * fx * pb_yz

                     + 0.5 * pa_xxz * fx * pb_zz

                     + 0.5 * fx * pa_yz * pb_yzz

                     + pa_xxyz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_yzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_xxyz,
                                    double pa_xxz,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double pb_z,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.25 * fx * fx * pa_z * fz * fgb

                     - 0.25 * fz * fga * fx * fx * pa_z

                     - 0.5 * fz * fga * fx * fx * pb_z

                     - 0.5 * pa_xxz * fx * fz * fgb

                     + fx * fx * fx * fz * pa_z

                     + 2.0 * fx * fx * fx * fz * pb_z

                     + 2.5 * pa_xxz * fx * fx * fz

                     + 5.0 * pa_xx * fx * fx * fz * pb_z

                     - 0.5 * fx * pa_yz * pb_y * fz * fgb

                     - 0.5 * fz * fga * pa_yz * pb_y * fx

                     - fz * fga * pa_y * fx * pb_yz

                     - 0.5 * fz * fga * fx * pa_z * pb_zz

                     - pa_xxyz * pb_y * fz * fgb

                     + 2.5 * fx * fx * pa_yz * fz * pb_y

                     + 5.0 * fx * fx * pa_y * fz * pb_yz

                     + 2.5 * fx * fx * fz * pa_z * pb_zz

                     + 6.0 * pa_xxyz * fz * pb_y * fx

                     + 12.0 * pa_xxy * fz * fx * pb_yz

                     + 6.0 * pa_xxz * fx * fz * pb_zz

                     - fz * fga * pa_yz * pb_yzz

                     + 6.0 * fx * pa_yz * fz * pb_yzz

                     + 14.0 * pa_xxyz * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_zzz_s_0(double fx,
                                    double pa_xxy,
                                    double pa_xxyz,
                                    double pa_y,
                                    double pa_yz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_y

                     + 0.75 * pa_xxy * fx * fx

                     + 0.75 * fx * fx * pa_yz * pb_z

                     + 0.75 * fx * fx * pa_y * pb_zz

                     + 1.5 * pa_xxyz * pb_z * fx

                     + 1.5 * pa_xxy * fx * pb_zz

                     + 0.5 * fx * pa_yz * pb_zzz

                     + pa_xxyz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_zzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxy,
                                    double pa_xxyz,
                                    double pa_y,
                                    double pa_yz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pa_y * fz * fgb

                     - 0.75 * fz * fga * pa_y * fx * fx

                     - 1.5 * pa_xxy * fx * fz * fgb

                     + 3.0 * fx * fx * fx * pa_y * fz

                     + 7.5 * pa_xxy * fz * fx * fx

                     - 1.5 * fx * pa_yz * pb_z * fz * fgb

                     - 1.5 * fz * fga * pa_yz * pb_z * fx

                     - 1.5 * fz * fga * pa_y * fx * pb_zz

                     - 3.0 * pa_xxyz * pb_z * fz * fgb

                     + 7.5 * fx * fx * pa_yz * fz * pb_z

                     + 7.5 * fx * fx * pa_y * fz * pb_zz

                     + 18.0 * pa_xxyz * fz * pb_z * fx

                     + 18.0 * pa_xxy * fz * fx * pb_zz

                     - fz * fga * pa_yz * pb_zzz

                     + 6.0 * fx * pa_yz * fz * pb_zzz

                     + 14.0 * pa_xxyz * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xxx_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxzz,
                                    double pa_xzz,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * fx

                     + 1.125 * fx * fx * fx * pb_x

                     + 1.5 * pa_xzz * fx * fx

                     + 2.25 * fx * fx * pa_zz * pb_x

                     + 0.75 * pa_xx * fx * fx * pb_x

                     + 1.5 * pa_x * fx * fx * pb_xx

                     + 1.5 * pa_xxzz * pb_x * fx

                     + 3.0 * pa_xzz * fx * pb_xx

                     + 0.25 * fx * fx * pb_xxx

                     + 0.5 * pa_xx * fx * pb_xxx

                     + 0.5 * fx * pa_zz * pb_xxx

                     + pa_xxzz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxzz,
                                    double pa_xzz,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fx * fz * fgb

                     - 1.5 * pa_x * fx * fx * fz * fga

                     - 3.0 * fx * fx * fz * fga * pb_x

                     - 3.0 * pa_xzz * fx * fz * fgb

                     + 6.0 * pa_x * fx * fx * fx * fz

                     + 9.0 * fx * fx * fx * fz * pb_x

                     + 15.0 * pa_xzz * fx * fx * fz

                     + 22.5 * fx * fx * pa_zz * fz * pb_x

                     - 0.75 * fx * fx * pb_x * fz * fgb

                     - 1.5 * pa_xx * fx * pb_x * fz * fgb

                     - 1.5 * pa_xx * fz * fga * pb_x * fx

                     - 3.0 * pa_x * fx * fz * fga * pb_xx

                     - 1.5 * fx * pa_zz * pb_x * fz * fgb

                     - 1.5 * fz * fga * pa_zz * pb_x * fx

                     - 3.0 * pa_xxzz * pb_x * fz * fgb

                     + 7.5 * pa_xx * fz * fx * fx * pb_x

                     + 15.0 * pa_x * fx * fx * fz * pb_xx

                     + 18.0 * pa_xxzz * fz * pb_x * fx

                     + 36.0 * pa_xzz * fx * fz * pb_xx

                     - fx * fz * fga * pb_xxx

                     - pa_xx * fz * fga * pb_xxx

                     + 2.5 * fx * fx * fz * pb_xxx

                     - fz * fga * pa_zz * pb_xxx

                     + 6.0 * pa_xx * fz * fx * pb_xxx

                     + 6.0 * fx * pa_zz * fz * pb_xxx

                     + 14.0 * pa_xxzz * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xxy_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxzz,
                                    double pa_xzz,
                                    double pa_zz,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_zz * pb_y

                     + 0.25 * pa_xx * fx * fx * pb_y

                     + pa_x * fx * fx * pb_xy

                     + 0.5 * pa_xxzz * fx * pb_y

                     + 2.0 * pa_xzz * fx * pb_xy

                     + 0.25 * fx * fx * pb_xxy

                     + 0.5 * pa_xx * fx * pb_xxy

                     + 0.5 * fx * pa_zz * pb_xxy

                     + pa_xxzz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxzz,
                                    double pa_xzz,
                                    double pa_zz,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double r_0_0)
    {
        return r_0_0 * (- fx * fx * fz * fga * pb_y

                     + 3.0 * fx * fx * fx * fz * pb_y

                     + 7.5 * fx * fx * pa_zz * fz * pb_y

                     - 0.25 * fx * fx * fz * fgb * pb_y

                     - 0.5 * pa_xx * fx * fz * fgb * pb_y

                     - 0.5 * pa_xx * fz * fga * fx * pb_y

                     - 2.0 * pa_x * fx * fz * fga * pb_xy

                     - 0.5 * fx * pa_zz * fz * fgb * pb_y

                     - 0.5 * fz * fga * pa_zz * fx * pb_y

                     - pa_xxzz * fz * fgb * pb_y

                     + 2.5 * pa_xx * fz * fx * fx * pb_y

                     + 10.0 * pa_x * fx * fx * fz * pb_xy

                     + 6.0 * pa_xxzz * fz * fx * pb_y

                     + 24.0 * pa_xzz * fx * fz * pb_xy

                     - fx * fz * fga * pb_xxy

                     - pa_xx * fz * fga * pb_xxy

                     + 2.5 * fx * fx * fz * pb_xxy

                     - fz * fga * pa_zz * pb_xxy

                     + 6.0 * pa_xx * fz * fx * pb_xxy

                     + 6.0 * fx * pa_zz * fz * pb_xxy

                     + 14.0 * pa_xxzz * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xxz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxz,
                                    double pa_xxzz,
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
        return s_0_0 * (0.75 * fx * fx * fx * pa_z

                     + 0.375 * fx * fx * fx * pb_z

                     + 0.5 * pa_xxz * fx * fx

                     + 2.0 * pa_xz * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_zz * pb_z

                     + 0.25 * pa_xx * fx * fx * pb_z

                     + pa_x * fx * fx * pb_xz

                     + 0.5 * fx * fx * pa_z * pb_xx

                     + 0.5 * pa_xxzz * fx * pb_z

                     + pa_xxz * fx * pb_xx

                     + 2.0 * pa_xzz * fx * pb_xz

                     + 0.25 * fx * fx * pb_xxz

                     + 0.5 * pa_xx * fx * pb_xxz

                     + 0.5 * fx * pa_zz * pb_xxz

                     + pa_xxzz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxz,
                                    double pa_xxzz,
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
        return r_0_0 * (6.0 * fx * fx * fx * pa_z * fz

                     - 0.5 * fx * fx * pa_z * fz * fgb

                     - fx * fx * fz * fga * pb_z

                     - 0.5 * fz * fga * pa_z * fx * fx

                     - pa_xxz * fx * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_z

                     + 5.0 * pa_xxz * fz * fx * fx

                     + 20.0 * pa_xz * fx * fx * fz * pb_x

                     + 7.5 * fx * fx * pa_zz * fz * pb_z

                     - 0.25 * fx * fx * fz * fgb * pb_z

                     - 0.5 * pa_xx * fx * fz * fgb * pb_z

                     - 0.5 * pa_xx * fz * fga * fx * pb_z

                     - 2.0 * pa_x * fx * fz * fga * pb_xz

                     - 0.5 * fx * pa_zz * fz * fgb * pb_z

                     - 0.5 * fz * fga * pa_zz * fx * pb_z

                     - fz * fga * pa_z * fx * pb_xx

                     - pa_xxzz * fz * fgb * pb_z

                     + 2.5 * pa_xx * fz * fx * fx * pb_z

                     + 10.0 * pa_x * fx * fx * fz * pb_xz

                     + 5.0 * fx * fx * pa_z * fz * pb_xx

                     + 6.0 * pa_xxzz * fz * fx * pb_z

                     + 12.0 * pa_xxz * fz * fx * pb_xx

                     + 24.0 * pa_xzz * fx * fz * pb_xz

                     - fx * fz * fga * pb_xxz

                     - pa_xx * fz * fga * pb_xxz

                     + 2.5 * fx * fx * fz * pb_xxz

                     - fz * fga * pa_zz * pb_xxz

                     + 6.0 * pa_xx * fz * fx * pb_xxz

                     + 6.0 * fx * pa_zz * fz * pb_xxz

                     + 14.0 * pa_xxzz * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xyy_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxzz,
                                    double pa_xzz,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xyy,
                                    double pb_yy,
                                    double s_0_0)
    {
        return s_0_0 * (0.25 * pa_x * fx * fx * fx

                     + 0.5 * pa_xzz * fx * fx

                     + 0.125 * fx * fx * fx * pb_x

                     + 0.25 * pa_xx * fx * fx * pb_x

                     + 0.5 * pa_x * fx * fx * pb_yy

                     + 0.25 * fx * fx * pa_zz * pb_x

                     + 0.5 * pa_xxzz * pb_x * fx

                     + pa_xzz * fx * pb_yy

                     + 0.25 * fx * fx * pb_xyy

                     + 0.5 * pa_xx * fx * pb_xyy

                     + 0.5 * fx * pa_zz * pb_xyy

                     + pa_xxzz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxzz,
                                    double pa_xzz,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xyy,
                                    double pb_yy,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_x * fx * fx * fz * fgb

                     - 0.5 * pa_x * fx * fx * fz * fga

                     - pa_xzz * fx * fz * fgb

                     + 2.0 * pa_x * fx * fx * fx * fz

                     + 5.0 * pa_xzz * fx * fx * fz

                     - 0.25 * fx * fx * pb_x * fz * fgb

                     - 0.5 * fx * fx * fz * fga * pb_x

                     - 0.5 * pa_xx * fx * pb_x * fz * fgb

                     - 0.5 * pa_xx * fz * fga * pb_x * fx

                     - pa_x * fx * fz * fga * pb_yy

                     - 0.5 * fx * pa_zz * pb_x * fz * fgb

                     + fx * fx * fx * fz * pb_x

                     - 0.5 * fz * fga * pa_zz * pb_x * fx

                     - pa_xxzz * pb_x * fz * fgb

                     + 2.5 * pa_xx * fz * fx * fx * pb_x

                     + 5.0 * pa_x * fx * fx * fz * pb_yy

                     + 2.5 * fx * fx * pa_zz * fz * pb_x

                     + 6.0 * pa_xxzz * fz * pb_x * fx

                     + 12.0 * pa_xzz * fx * fz * pb_yy

                     - fx * fz * fga * pb_xyy

                     - pa_xx * fz * fga * pb_xyy

                     + 2.5 * fx * fx * fz * pb_xyy

                     - fz * fga * pa_zz * pb_xyy

                     + 6.0 * pa_xx * fz * fx * pb_xyy

                     + 6.0 * fx * pa_zz * fz * pb_xyy

                     + 14.0 * pa_xxzz * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxz,
                                    double pa_xxzz,
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
        return s_0_0 * (pa_xz * fx * fx * pb_y

                     + 0.5 * pa_x * fx * fx * pb_yz

                     + 0.5 * fx * fx * pa_z * pb_xy

                     + pa_xxz * fx * pb_xy

                     + pa_xzz * fx * pb_yz

                     + 0.25 * fx * fx * pb_xyz

                     + 0.5 * pa_xx * fx * pb_xyz

                     + 0.5 * fx * pa_zz * pb_xyz

                     + pa_xxzz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xyz_r_0(double fga,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxz,
                                    double pa_xxzz,
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
        return r_0_0 * (10.0 * pa_xz * fx * fx * fz * pb_y

                     - pa_x * fx * fz * fga * pb_yz

                     - fz * fga * pa_z * fx * pb_xy

                     + 5.0 * pa_x * fx * fx * fz * pb_yz

                     + 5.0 * fx * fx * pa_z * fz * pb_xy

                     + 12.0 * pa_xxz * fz * fx * pb_xy

                     + 12.0 * pa_xzz * fx * fz * pb_yz

                     - fx * fz * fga * pb_xyz

                     - pa_xx * fz * fga * pb_xyz

                     + 2.5 * fx * fx * fz * pb_xyz

                     - fz * fga * pa_zz * pb_xyz

                     + 6.0 * pa_xx * fz * fx * pb_xyz

                     + 6.0 * fx * pa_zz * fz * pb_xyz

                     + 14.0 * pa_xxzz * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxz,
                                    double pa_xxzz,
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
        return s_0_0 * (0.75 * pa_x * fx * fx * fx

                     + 0.375 * fx * fx * fx * pb_x

                     + 0.75 * pa_xx * fx * fx * pb_x

                     + 0.5 * pa_xzz * fx * fx

                     + 2.0 * pa_xz * fx * fx * pb_z

                     + 0.5 * pa_x * fx * fx * pb_zz

                     + 0.25 * fx * fx * pa_zz * pb_x

                     + fx * fx * pa_z * pb_xz

                     + 0.5 * pa_xxzz * pb_x * fx

                     + 2.0 * pa_xxz * fx * pb_xz

                     + pa_xzz * fx * pb_zz

                     + 0.25 * fx * fx * pb_xzz

                     + 0.5 * pa_xx * fx * pb_xzz

                     + 0.5 * fx * pa_zz * pb_xzz

                     + pa_xxzz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxz,
                                    double pa_xxzz,
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
        return r_0_0 * (6.0 * pa_x * fx * fx * fx * fz

                     - 0.5 * pa_x * fx * fx * fz * fgb

                     - 0.5 * pa_x * fx * fx * fz * fga

                     - fz * fga * fx * fx * pb_x

                     - pa_xzz * fx * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_x

                     + 7.5 * pa_xx * fx * fx * fz * pb_x

                     + 5.0 * pa_xzz * fx * fx * fz

                     + 20.0 * pa_xz * fx * fx * fz * pb_z

                     - 0.25 * fx * fx * pb_x * fz * fgb

                     - 0.5 * pa_xx * fx * pb_x * fz * fgb

                     - 0.5 * pa_xx * fz * fga * pb_x * fx

                     - pa_x * fx * fz * fga * pb_zz

                     - 0.5 * fx * pa_zz * pb_x * fz * fgb

                     - 0.5 * fz * fga * pa_zz * pb_x * fx

                     - 2.0 * fz * fga * pa_z * fx * pb_xz

                     - pa_xxzz * pb_x * fz * fgb

                     + 5.0 * pa_x * fx * fx * fz * pb_zz

                     + 2.5 * fx * fx * pa_zz * fz * pb_x

                     + 10.0 * fx * fx * pa_z * fz * pb_xz

                     + 6.0 * pa_xxzz * fz * pb_x * fx

                     + 24.0 * pa_xxz * fz * fx * pb_xz

                     + 12.0 * pa_xzz * fx * fz * pb_zz

                     - fx * fz * fga * pb_xzz

                     - pa_xx * fz * fga * pb_xzz

                     + 2.5 * fx * fx * fz * pb_xzz

                     - fz * fga * pa_zz * pb_xzz

                     + 6.0 * pa_xx * fz * fx * pb_xzz

                     + 6.0 * fx * pa_zz * fz * pb_xzz

                     + 14.0 * pa_xxzz * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_yyy_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxzz,
                                    double pa_zz,
                                    double pb_y,
                                    double pb_yyy,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_y

                     + 0.75 * pa_xx * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_zz * pb_y

                     + 1.5 * pa_xxzz * pb_y * fx

                     + 0.25 * fx * fx * pb_yyy

                     + 0.5 * pa_xx * fx * pb_yyy

                     + 0.5 * fx * pa_zz * pb_yyy

                     + pa_xxzz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_yyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxzz,
                                    double pa_zz,
                                    double pb_y,
                                    double pb_yyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pb_y * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_y

                     - 1.5 * pa_xx * fx * pb_y * fz * fgb

                     - 1.5 * pa_xx * fz * fga * pb_y * fx

                     - 1.5 * fx * pa_zz * pb_y * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_y

                     - 1.5 * fz * fga * pa_zz * pb_y * fx

                     - 3.0 * pa_xxzz * pb_y * fz * fgb

                     + 7.5 * pa_xx * fz * fx * fx * pb_y

                     + 7.5 * fx * fx * pa_zz * fz * pb_y

                     + 18.0 * pa_xxzz * fz * pb_y * fx

                     - fx * fz * fga * pb_yyy

                     - pa_xx * fz * fga * pb_yyy

                     + 2.5 * fx * fx * fz * pb_yyy

                     - fz * fga * pa_zz * pb_yyy

                     + 6.0 * pa_xx * fz * fx * pb_yyy

                     + 6.0 * fx * pa_zz * fz * pb_yyy

                     + 14.0 * pa_xxzz * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_yyz_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxz,
                                    double pa_xxzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * fx * pa_z

                     + 0.5 * pa_xxz * fx * fx

                     + 0.125 * fx * fx * fx * pb_z

                     + 0.25 * pa_xx * fx * fx * pb_z

                     + 0.25 * fx * fx * pa_zz * pb_z

                     + 0.5 * fx * fx * pa_z * pb_yy

                     + 0.5 * pa_xxzz * fx * pb_z

                     + pa_xxz * fx * pb_yy

                     + 0.25 * fx * fx * pb_yyz

                     + 0.5 * pa_xx * fx * pb_yyz

                     + 0.5 * fx * pa_zz * pb_yyz

                     + pa_xxzz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_yyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxz,
                                    double pa_xxzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * pa_z * fz * fgb

                     - 0.5 * fz * fga * pa_z * fx * fx

                     - pa_xxz * fx * fz * fgb

                     + 2.0 * fx * fx * fx * pa_z * fz

                     + 5.0 * pa_xxz * fz * fx * fx

                     - 0.25 * fx * fx * fz * fgb * pb_z

                     - 0.5 * fx * fx * fz * fga * pb_z

                     - 0.5 * pa_xx * fx * fz * fgb * pb_z

                     - 0.5 * pa_xx * fz * fga * fx * pb_z

                     - 0.5 * fx * pa_zz * fz * fgb * pb_z

                     + fx * fx * fx * fz * pb_z

                     - 0.5 * fz * fga * pa_zz * fx * pb_z

                     - fz * fga * pa_z * fx * pb_yy

                     - pa_xxzz * fz * fgb * pb_z

                     + 2.5 * pa_xx * fz * fx * fx * pb_z

                     + 2.5 * fx * fx * pa_zz * fz * pb_z

                     + 5.0 * fx * fx * pa_z * fz * pb_yy

                     + 6.0 * pa_xxzz * fz * fx * pb_z

                     + 12.0 * pa_xxz * fz * fx * pb_yy

                     - fx * fz * fga * pb_yyz

                     - pa_xx * fz * fga * pb_yyz

                     + 2.5 * fx * fx * fz * pb_yyz

                     - fz * fga * pa_zz * pb_yyz

                     + 6.0 * pa_xx * fz * fx * pb_yyz

                     + 6.0 * fx * pa_zz * fz * pb_yyz

                     + 14.0 * pa_xxzz * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_yzz_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxz,
                                    double pa_xxzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_y

                     + 0.75 * pa_xx * fx * fx * pb_y

                     + 0.25 * fx * fx * pa_zz * pb_y

                     + fx * fx * pa_z * pb_yz

                     + 0.5 * pa_xxzz * pb_y * fx

                     + 2.0 * pa_xxz * fx * pb_yz

                     + 0.25 * fx * fx * pb_yzz

                     + 0.5 * pa_xx * fx * pb_yzz

                     + 0.5 * fx * pa_zz * pb_yzz

                     + pa_xxzz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_yzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxz,
                                    double pa_xxzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double r_0_0)
    {
        return r_0_0 * (- fz * fga * fx * fx * pb_y

                     + 3.0 * fx * fx * fx * fz * pb_y

                     + 7.5 * pa_xx * fx * fx * fz * pb_y

                     - 0.25 * fx * fx * pb_y * fz * fgb

                     - 0.5 * pa_xx * fx * pb_y * fz * fgb

                     - 0.5 * pa_xx * fz * fga * pb_y * fx

                     - 0.5 * fx * pa_zz * pb_y * fz * fgb

                     - 0.5 * fz * fga * pa_zz * pb_y * fx

                     - 2.0 * fz * fga * pa_z * fx * pb_yz

                     - pa_xxzz * pb_y * fz * fgb

                     + 2.5 * fx * fx * pa_zz * fz * pb_y

                     + 10.0 * fx * fx * pa_z * fz * pb_yz

                     + 6.0 * pa_xxzz * fz * pb_y * fx

                     + 24.0 * pa_xxz * fz * fx * pb_yz

                     - fx * fz * fga * pb_yzz

                     - pa_xx * fz * fga * pb_yzz

                     + 2.5 * fx * fx * fz * pb_yzz

                     - fz * fga * pa_zz * pb_yzz

                     + 6.0 * pa_xx * fz * fx * pb_yzz

                     + 6.0 * fx * pa_zz * fz * pb_yzz

                     + 14.0 * pa_xxzz * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_zzz_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxz,
                                    double pa_xxzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * fx * pa_z

                     + 1.125 * fx * fx * fx * pb_z

                     + 1.5 * pa_xxz * fx * fx

                     + 2.25 * pa_xx * fx * fx * pb_z

                     + 0.75 * fx * fx * pa_zz * pb_z

                     + 1.5 * fx * fx * pa_z * pb_zz

                     + 1.5 * pa_xxzz * pb_z * fx

                     + 3.0 * pa_xxz * fx * pb_zz

                     + 0.25 * fx * fx * pb_zzz

                     + 0.5 * pa_xx * fx * pb_zzz

                     + 0.5 * fx * pa_zz * pb_zzz

                     + pa_xxzz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_zzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxz,
                                    double pa_xxzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_z * fz * fgb

                     - 1.5 * fz * fga * pa_z * fx * fx

                     - 3.0 * fz * fga * fx * fx * pb_z

                     - 3.0 * pa_xxz * fx * fz * fgb

                     + 6.0 * fx * fx * fx * pa_z * fz

                     + 9.0 * fx * fx * fx * fz * pb_z

                     + 15.0 * pa_xxz * fz * fx * fx

                     + 22.5 * pa_xx * fx * fx * fz * pb_z

                     - 0.75 * fx * fx * pb_z * fz * fgb

                     - 1.5 * pa_xx * fx * pb_z * fz * fgb

                     - 1.5 * pa_xx * fz * fga * pb_z * fx

                     - 1.5 * fx * pa_zz * pb_z * fz * fgb

                     - 1.5 * fz * fga * pa_zz * pb_z * fx

                     - 3.0 * fz * fga * pa_z * fx * pb_zz

                     - 3.0 * pa_xxzz * pb_z * fz * fgb

                     + 7.5 * fx * fx * pa_zz * fz * pb_z

                     + 15.0 * fx * fx * pa_z * fz * pb_zz

                     + 18.0 * pa_xxzz * fz * pb_z * fx

                     + 36.0 * pa_xxz * fz * fx * pb_zz

                     - fx * fz * fga * pb_zzz

                     - pa_xx * fz * fga * pb_zzz

                     + 2.5 * fx * fx * fz * pb_zzz

                     - fz * fga * pa_zz * pb_zzz

                     + 6.0 * pa_xx * fz * fx * pb_zzz

                     + 6.0 * fx * pa_zz * fz * pb_zzz

                     + 14.0 * pa_xxzz * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xxx_s_0(double fx,
                                    double pa_xy,
                                    double pa_xyyy,
                                    double pa_y,
                                    double pa_yyy,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pa_y

                     + 0.75 * fx * fx * pa_yyy

                     + 2.25 * pa_xy * fx * fx * pb_x

                     + 2.25 * fx * fx * pa_y * pb_xx

                     + 1.5 * pa_xyyy * pb_x * fx

                     + 1.5 * fx * pa_yyy * pb_xx

                     + 1.5 * pa_xy * fx * pb_xxx

                     + pa_xyyy * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xy,
                                    double pa_xyyy,
                                    double pa_y,
                                    double pa_yyy,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pa_y * fz * fgb

                     - 2.25 * fx * fx * pa_y * fz * fga

                     - 1.5 * fx * pa_yyy * fz * fgb

                     + 9.0 * fx * fx * fx * pa_y * fz

                     + 7.5 * fx * fx * pa_yyy * fz

                     - 4.5 * pa_xy * fx * pb_x * fz * fgb

                     - 4.5 * pa_xy * fz * fga * pb_x * fx

                     - 4.5 * fx * pa_y * fz * fga * pb_xx

                     - 3.0 * pa_xyyy * pb_x * fz * fgb

                     + 22.5 * pa_xy * fz * fx * fx * pb_x

                     + 22.5 * fx * fx * pa_y * fz * pb_xx

                     + 18.0 * pa_xyyy * fz * pb_x * fx

                     + 18.0 * fx * pa_yyy * fz * pb_xx

                     - 3.0 * pa_xy * fz * fga * pb_xxx

                     + 18.0 * pa_xy * fz * fx * pb_xxx

                     + 14.0 * pa_xyyy * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xxy_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pa_xyyy,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 0.75 * fx * fx * fx * pb_x

                     + 0.75 * pa_xyy * fx * fx

                     + 1.5 * fx * fx * pa_yy * pb_x

                     + 0.75 * pa_xy * fx * fx * pb_y

                     + 0.75 * pa_x * fx * fx * pb_xx

                     + 1.5 * fx * fx * pa_y * pb_xy

                     + 0.5 * pa_xyyy * fx * pb_y

                     + 1.5 * pa_xyy * fx * pb_xx

                     + fx * pa_yyy * pb_xy

                     + 1.5 * pa_xy * fx * pb_xxy

                     + pa_xyyy * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pa_xyyy,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * pa_x * fx * fx * fz * fgb

                     - 0.75 * pa_x * fx * fx * fz * fga

                     - 1.5 * fx * fx * fz * fga * pb_x

                     - 1.5 * pa_xyy * fx * fz * fgb

                     + 3.0 * pa_x * fx * fx * fx * fz

                     + 6.0 * fx * fx * fx * fz * pb_x

                     + 7.5 * pa_xyy * fz * fx * fx

                     + 15.0 * fx * fx * pa_yy * fz * pb_x

                     - 1.5 * pa_xy * fx * fz * fgb * pb_y

                     - 1.5 * pa_xy * fz * fga * fx * pb_y

                     - 1.5 * pa_x * fx * fz * fga * pb_xx

                     - 3.0 * fx * pa_y * fz * fga * pb_xy

                     - pa_xyyy * fz * fgb * pb_y

                     + 7.5 * pa_xy * fz * fx * fx * pb_y

                     + 7.5 * pa_x * fx * fx * fz * pb_xx

                     + 15.0 * fx * fx * pa_y * fz * pb_xy

                     + 6.0 * pa_xyyy * fz * fx * pb_y

                     + 18.0 * pa_xyy * fz * fx * pb_xx

                     + 12.0 * fx * pa_yyy * fz * pb_xy

                     - 3.0 * pa_xy * fz * fga * pb_xxy

                     + 18.0 * pa_xy * fz * fx * pb_xxy

                     + 14.0 * pa_xyyy * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xxz_s_0(double fx,
                                    double pa_xy,
                                    double pa_xyyy,
                                    double pa_y,
                                    double pa_yyy,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xy * fx * fx * pb_z

                     + 1.5 * fx * fx * pa_y * pb_xz

                     + 0.5 * pa_xyyy * fx * pb_z

                     + fx * pa_yyy * pb_xz

                     + 1.5 * pa_xy * fx * pb_xxz

                     + pa_xyyy * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xy,
                                    double pa_xyyy,
                                    double pa_y,
                                    double pa_yyy,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xy * fx * fz * fgb * pb_z

                     - 1.5 * pa_xy * fz * fga * fx * pb_z

                     - 3.0 * fx * pa_y * fz * fga * pb_xz

                     - pa_xyyy * fz * fgb * pb_z

                     + 7.5 * pa_xy * fz * fx * fx * pb_z

                     + 15.0 * fx * fx * pa_y * fz * pb_xz

                     + 6.0 * pa_xyyy * fz * fx * pb_z

                     + 12.0 * fx * pa_yyy * fz * pb_xz

                     - 3.0 * pa_xy * fz * fga * pb_xxz

                     + 18.0 * pa_xy * fz * fx * pb_xxz

                     + 14.0 * pa_xyyy * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xyy_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pa_xyyy,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyy,
                                    double pb_y,
                                    double pb_yy,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pa_y

                     + 0.75 * fx * fx * fx * pb_y

                     + 2.25 * pa_xy * fx * fx * pb_x

                     + 0.25 * fx * fx * pa_yyy

                     + 1.5 * fx * fx * pa_yy * pb_y

                     + 1.5 * pa_x * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_y * pb_yy

                     + 0.5 * pa_xyyy * pb_x * fx

                     + 3.0 * pa_xyy * fx * pb_xy

                     + 0.5 * fx * pa_yyy * pb_yy

                     + 1.5 * pa_xy * fx * pb_xyy

                     + pa_xyyy * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pa_xyyy,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyy,
                                    double pb_y,
                                    double pb_yy,
                                    double r_0_0)
    {
        return r_0_0 * (9.0 * fx * fx * fx * pa_y * fz

                     - 0.75 * fx * fx * pa_y * fz * fgb

                     - 0.75 * fx * fx * pa_y * fz * fga

                     - 1.5 * fx * fx * fz * fga * pb_y

                     - 0.5 * fx * pa_yyy * fz * fgb

                     + 6.0 * fx * fx * fx * fz * pb_y

                     + 22.5 * pa_xy * fx * fx * fz * pb_x

                     + 2.5 * fx * fx * pa_yyy * fz

                     + 15.0 * fx * fx * pa_yy * fz * pb_y

                     - 1.5 * pa_xy * fx * pb_x * fz * fgb

                     - 1.5 * pa_xy * fz * fga * pb_x * fx

                     - 3.0 * pa_x * fx * fz * fga * pb_xy

                     - 1.5 * fx * pa_y * fz * fga * pb_yy

                     - pa_xyyy * pb_x * fz * fgb

                     + 15.0 * pa_x * fx * fx * fz * pb_xy

                     + 7.5 * fx * fx * pa_y * fz * pb_yy

                     + 6.0 * pa_xyyy * fz * pb_x * fx

                     + 36.0 * pa_xyy * fz * fx * pb_xy

                     + 6.0 * fx * pa_yyy * fz * pb_yy

                     - 3.0 * pa_xy * fz * fga * pb_xyy

                     + 18.0 * pa_xy * fz * fx * pb_xyy

                     + 14.0 * pa_xyyy * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pa_xyyy,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_xyz,
                                    double pb_xz,
                                    double pb_yz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_z

                     + 0.75 * fx * fx * pa_yy * pb_z

                     + 0.75 * pa_x * fx * fx * pb_xz

                     + 0.75 * fx * fx * pa_y * pb_yz

                     + 1.5 * pa_xyy * fx * pb_xz

                     + 0.5 * fx * pa_yyy * pb_yz

                     + 1.5 * pa_xy * fx * pb_xyz

                     + pa_xyyy * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xyz_r_0(double fga,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pa_xyyy,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_xyz,
                                    double pb_xz,
                                    double pb_yz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fga * pb_z

                     + 3.0 * fx * fx * fx * fz * pb_z

                     + 7.5 * fx * fx * pa_yy * fz * pb_z

                     - 1.5 * pa_x * fx * fz * fga * pb_xz

                     - 1.5 * fx * pa_y * fz * fga * pb_yz

                     + 7.5 * pa_x * fx * fx * fz * pb_xz

                     + 7.5 * fx * fx * pa_y * fz * pb_yz

                     + 18.0 * pa_xyy * fz * fx * pb_xz

                     + 6.0 * fx * pa_yyy * fz * pb_yz

                     - 3.0 * pa_xy * fz * fga * pb_xyz

                     + 18.0 * pa_xy * fz * fx * pb_xyz

                     + 14.0 * pa_xyyy * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xzz_s_0(double fx,
                                    double pa_xy,
                                    double pa_xyyy,
                                    double pa_y,
                                    double pa_yyy,
                                    double pb_x,
                                    double pb_xzz,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_y

                     + 0.25 * fx * fx * pa_yyy

                     + 0.75 * pa_xy * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_y * pb_zz

                     + 0.5 * pa_xyyy * pb_x * fx

                     + 0.5 * fx * pa_yyy * pb_zz

                     + 1.5 * pa_xy * fx * pb_xzz

                     + pa_xyyy * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xy,
                                    double pa_xyyy,
                                    double pa_y,
                                    double pa_yyy,
                                    double pb_x,
                                    double pb_xzz,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pa_y * fz * fgb

                     - 0.75 * fx * fx * pa_y * fz * fga

                     - 0.5 * fx * pa_yyy * fz * fgb

                     + 3.0 * fx * fx * fx * pa_y * fz

                     + 2.5 * fx * fx * pa_yyy * fz

                     - 1.5 * pa_xy * fx * pb_x * fz * fgb

                     - 1.5 * pa_xy * fz * fga * pb_x * fx

                     - 1.5 * fx * pa_y * fz * fga * pb_zz

                     - pa_xyyy * pb_x * fz * fgb

                     + 7.5 * pa_xy * fz * fx * fx * pb_x

                     + 7.5 * fx * fx * pa_y * fz * pb_zz

                     + 6.0 * pa_xyyy * fz * pb_x * fx

                     + 6.0 * fx * pa_yyy * fz * pb_zz

                     - 3.0 * pa_xy * fz * fga * pb_xzz

                     + 18.0 * pa_xy * fz * fx * pb_xzz

                     + 14.0 * pa_xyyy * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_yyy_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pa_xyyy,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * pa_x * fx * fx * fx

                     + 2.25 * pa_xyy * fx * fx

                     + 6.75 * pa_xy * fx * fx * pb_y

                     + 2.25 * pa_x * fx * fx * pb_yy

                     + 1.5 * pa_xyyy * pb_y * fx

                     + 4.5 * pa_xyy * fx * pb_yy

                     + 1.5 * pa_xy * fx * pb_yyy

                     + pa_xyyy * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_yyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pa_xyyy,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double r_0_0)
    {
        return r_0_0 * (15.0 * pa_x * fx * fx * fx * fz

                     - 2.25 * pa_x * fx * fx * fz * fgb

                     - 2.25 * pa_x * fx * fx * fz * fga

                     - 4.5 * pa_xyy * fx * fz * fgb

                     + 22.5 * pa_xyy * fz * fx * fx

                     + 67.5 * pa_xy * fx * fx * fz * pb_y

                     - 4.5 * pa_xy * fx * pb_y * fz * fgb

                     - 4.5 * pa_xy * fz * fga * pb_y * fx

                     - 4.5 * pa_x * fx * fz * fga * pb_yy

                     - 3.0 * pa_xyyy * pb_y * fz * fgb

                     + 22.5 * pa_x * fx * fx * fz * pb_yy

                     + 18.0 * pa_xyyy * fz * pb_y * fx

                     + 54.0 * pa_xyy * fz * fx * pb_yy

                     - 3.0 * pa_xy * fz * fga * pb_yyy

                     + 18.0 * pa_xy * fz * fx * pb_yyy

                     + 14.0 * pa_xyyy * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_yyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pa_xyyy,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_xy * fx * fx * pb_z

                     + 1.5 * pa_x * fx * fx * pb_yz

                     + 0.5 * pa_xyyy * fx * pb_z

                     + 3.0 * pa_xyy * fx * pb_yz

                     + 1.5 * pa_xy * fx * pb_yyz

                     + pa_xyyy * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_yyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pa_xyyy,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (22.5 * pa_xy * fx * fx * fz * pb_z

                     - 1.5 * pa_xy * fx * fz * fgb * pb_z

                     - 1.5 * pa_xy * fz * fga * fx * pb_z

                     - 3.0 * pa_x * fx * fz * fga * pb_yz

                     - pa_xyyy * fz * fgb * pb_z

                     + 15.0 * pa_x * fx * fx * fz * pb_yz

                     + 6.0 * pa_xyyy * fz * fx * pb_z

                     + 36.0 * pa_xyy * fz * fx * pb_yz

                     - 3.0 * pa_xy * fz * fga * pb_yyz

                     + 18.0 * pa_xy * fz * fx * pb_yyz

                     + 14.0 * pa_xyyy * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_yzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pa_xyyy,
                                    double pb_y,
                                    double pb_yzz,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 0.75 * pa_xyy * fx * fx

                     + 0.75 * pa_xy * fx * fx * pb_y

                     + 0.75 * pa_x * fx * fx * pb_zz

                     + 0.5 * pa_xyyy * pb_y * fx

                     + 1.5 * pa_xyy * fx * pb_zz

                     + 1.5 * pa_xy * fx * pb_yzz

                     + pa_xyyy * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_yzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pa_xyyy,
                                    double pb_y,
                                    double pb_yzz,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * pa_x * fx * fx * fz * fgb

                     - 0.75 * pa_x * fx * fx * fz * fga

                     - 1.5 * pa_xyy * fx * fz * fgb

                     + 3.0 * pa_x * fx * fx * fx * fz

                     + 7.5 * pa_xyy * fz * fx * fx

                     - 1.5 * pa_xy * fx * pb_y * fz * fgb

                     - 1.5 * pa_xy * fz * fga * pb_y * fx

                     - 1.5 * pa_x * fx * fz * fga * pb_zz

                     - pa_xyyy * pb_y * fz * fgb

                     + 7.5 * pa_xy * fz * fx * fx * pb_y

                     + 7.5 * pa_x * fx * fx * fz * pb_zz

                     + 6.0 * pa_xyyy * fz * pb_y * fx

                     + 18.0 * pa_xyy * fz * fx * pb_zz

                     - 3.0 * pa_xy * fz * fga * pb_yzz

                     + 18.0 * pa_xy * fz * fx * pb_yzz

                     + 14.0 * pa_xyyy * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_zzz_s_0(double fx,
                                    double pa_xy,
                                    double pa_xyyy,
                                    double pb_z,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_xy * fx * fx * pb_z

                     + 1.5 * pa_xyyy * pb_z * fx

                     + 1.5 * pa_xy * fx * pb_zzz

                     + pa_xyyy * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_zzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xy,
                                    double pa_xyyy,
                                    double pb_z,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_xy * fx * pb_z * fz * fgb

                     - 4.5 * pa_xy * fz * fga * pb_z * fx

                     - 3.0 * pa_xyyy * pb_z * fz * fgb

                     + 22.5 * pa_xy * fz * fx * fx * pb_z

                     + 18.0 * pa_xyyy * fz * pb_z * fx

                     - 3.0 * pa_xy * fz * fga * pb_zzz

                     + 18.0 * pa_xy * fz * fx * pb_zzz

                     + 14.0 * pa_xyyy * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xxx_s_0(double fx,
                                    double pa_xyyz,
                                    double pa_xz,
                                    double pa_yyz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_z

                     + 0.75 * fx * fx * pa_yyz

                     + 0.75 * pa_xz * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_z * pb_xx

                     + 1.5 * pa_xyyz * pb_x * fx

                     + 1.5 * fx * pa_yyz * pb_xx

                     + 0.5 * pa_xz * fx * pb_xxx

                     + pa_xyyz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xyyz,
                                    double pa_xz,
                                    double pa_yyz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pa_z * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pa_z

                     - 1.5 * fx * pa_yyz * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_z

                     + 7.5 * fx * fx * pa_yyz * fz

                     - 1.5 * pa_xz * fx * pb_x * fz * fgb

                     - 1.5 * pa_xz * fz * fga * pb_x * fx

                     - 1.5 * fx * fz * fga * pa_z * pb_xx

                     - 3.0 * pa_xyyz * pb_x * fz * fgb

                     + 7.5 * pa_xz * fx * fx * fz * pb_x

                     + 7.5 * fx * fx * fz * pa_z * pb_xx

                     + 18.0 * pa_xyyz * fz * pb_x * fx

                     + 18.0 * fx * pa_yyz * fz * pb_xx

                     - pa_xz * fz * fga * pb_xxx

                     + 6.0 * pa_xz * fx * fz * pb_xxx

                     + 14.0 * pa_xyyz * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xxy_s_0(double fx,
                                    double pa_xyyz,
                                    double pa_xyz,
                                    double pa_xz,
                                    double pa_yyz,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xyz * fx * fx

                     + fx * fx * pa_yz * pb_x

                     + 0.25 * pa_xz * fx * fx * pb_y

                     + 0.5 * fx * fx * pa_z * pb_xy

                     + 0.5 * pa_xyyz * fx * pb_y

                     + pa_xyz * fx * pb_xx

                     + fx * pa_yyz * pb_xy

                     + 0.5 * pa_xz * fx * pb_xxy

                     + pa_xyyz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xyyz,
                                    double pa_xyz,
                                    double pa_xz,
                                    double pa_yyz,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double r_0_0)
    {
        return r_0_0 * (- pa_xyz * fx * fz * fgb

                     + 5.0 * pa_xyz * fx * fx * fz

                     + 10.0 * fx * fx * pa_yz * fz * pb_x

                     - 0.5 * pa_xz * fx * fz * fgb * pb_y

                     - 0.5 * pa_xz * fz * fga * fx * pb_y

                     - fx * fz * fga * pa_z * pb_xy

                     - pa_xyyz * fz * fgb * pb_y

                     + 2.5 * pa_xz * fx * fx * fz * pb_y

                     + 5.0 * fx * fx * fz * pa_z * pb_xy

                     + 6.0 * pa_xyyz * fz * fx * pb_y

                     + 12.0 * pa_xyz * fx * fz * pb_xx

                     + 12.0 * fx * pa_yyz * fz * pb_xy

                     - pa_xz * fz * fga * pb_xxy

                     + 6.0 * pa_xz * fx * fz * pb_xxy

                     + 14.0 * pa_xyyz * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xxz_s_0(double fx,
                                    double pa_x,
                                    double pa_xyy,
                                    double pa_xyyz,
                                    double pa_xz,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.125 * pa_x * fx * fx * fx

                     + 0.25 * fx * fx * fx * pb_x

                     + 0.25 * pa_xyy * fx * fx

                     + 0.5 * fx * fx * pa_yy * pb_x

                     + 0.25 * pa_xz * fx * fx * pb_z

                     + 0.25 * pa_x * fx * fx * pb_xx

                     + 0.5 * fx * fx * pa_z * pb_xz

                     + 0.5 * pa_xyyz * fx * pb_z

                     + 0.5 * pa_xyy * fx * pb_xx

                     + fx * pa_yyz * pb_xz

                     + 0.5 * pa_xz * fx * pb_xxz

                     + pa_xyyz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xyy,
                                    double pa_xyyz,
                                    double pa_xz,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.25 * pa_x * fx * fx * fz * fgb

                     - 0.25 * pa_x * fz * fga * fx * fx

                     - 0.5 * fx * fx * fz * fga * pb_x

                     - 0.5 * pa_xyy * fx * fz * fgb

                     + pa_x * fx * fx * fx * fz

                     + 2.0 * fx * fx * fx * fz * pb_x

                     + 2.5 * pa_xyy * fz * fx * fx

                     + 5.0 * fx * fx * pa_yy * fz * pb_x

                     - 0.5 * pa_xz * fx * fz * fgb * pb_z

                     - 0.5 * pa_xz * fz * fga * fx * pb_z

                     - 0.5 * pa_x * fz * fga * fx * pb_xx

                     - fx * fz * fga * pa_z * pb_xz

                     - pa_xyyz * fz * fgb * pb_z

                     + 2.5 * pa_xz * fx * fx * fz * pb_z

                     + 2.5 * pa_x * fx * fx * fz * pb_xx

                     + 5.0 * fx * fx * fz * pa_z * pb_xz

                     + 6.0 * pa_xyyz * fz * fx * pb_z

                     + 6.0 * pa_xyy * fz * fx * pb_xx

                     + 12.0 * fx * pa_yyz * fz * pb_xz

                     - pa_xz * fz * fga * pb_xxz

                     + 6.0 * pa_xz * fx * fz * pb_xxz

                     + 14.0 * pa_xyyz * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xyy_s_0(double fx,
                                    double pa_xyyz,
                                    double pa_xyz,
                                    double pa_xz,
                                    double pa_yyz,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyy,
                                    double pb_y,
                                    double pb_yy,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_z

                     + 0.75 * pa_xz * fx * fx * pb_x

                     + 0.25 * fx * fx * pa_yyz

                     + fx * fx * pa_yz * pb_y

                     + 0.25 * fx * fx * pa_z * pb_yy

                     + 0.5 * pa_xyyz * pb_x * fx

                     + 2.0 * pa_xyz * fx * pb_xy

                     + 0.5 * fx * pa_yyz * pb_yy

                     + 0.5 * pa_xz * fx * pb_xyy

                     + pa_xyyz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xyyz,
                                    double pa_xyz,
                                    double pa_xz,
                                    double pa_yyz,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyy,
                                    double pb_y,
                                    double pb_yy,
                                    double r_0_0)
    {
        return r_0_0 * (3.0 * fx * fx * fx * fz * pa_z

                     - 0.25 * fx * fx * pa_z * fz * fgb

                     - 0.25 * fx * fx * fz * fga * pa_z

                     - 0.5 * fx * pa_yyz * fz * fgb

                     + 7.5 * pa_xz * fx * fx * fz * pb_x

                     + 2.5 * fx * fx * pa_yyz * fz

                     + 10.0 * fx * fx * pa_yz * fz * pb_y

                     - 0.5 * pa_xz * fx * pb_x * fz * fgb

                     - 0.5 * pa_xz * fz * fga * pb_x * fx

                     - 0.5 * fx * fz * fga * pa_z * pb_yy

                     - pa_xyyz * pb_x * fz * fgb

                     + 2.5 * fx * fx * fz * pa_z * pb_yy

                     + 6.0 * pa_xyyz * fz * pb_x * fx

                     + 24.0 * pa_xyz * fx * fz * pb_xy

                     + 6.0 * fx * pa_yyz * fz * pb_yy

                     - pa_xz * fz * fga * pb_xyy

                     + 6.0 * pa_xz * fx * fz * pb_xyy

                     + 14.0 * pa_xyyz * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pa_xyyz,
                                    double pa_xyz,
                                    double pa_xz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyz,
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
        return s_0_0 * (0.25 * fx * fx * fx * pa_y

                     + 0.125 * fx * fx * fx * pb_y

                     + 0.5 * pa_xy * fx * fx * pb_x

                     + 0.25 * fx * fx * pa_yy * pb_y

                     + 0.5 * fx * fx * pa_yz * pb_z

                     + 0.25 * pa_x * fx * fx * pb_xy

                     + 0.25 * fx * fx * pa_z * pb_yz

                     + 0.5 * pa_xyy * fx * pb_xy

                     + pa_xyz * fx * pb_xz

                     + 0.5 * fx * pa_yyz * pb_yz

                     + 0.5 * pa_xz * fx * pb_xyz

                     + pa_xyyz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xyz_r_0(double fga,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pa_xyyz,
                                    double pa_xyz,
                                    double pa_xz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyz,
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
        return r_0_0 * (2.0 * fx * fx * fx * pa_y * fz

                     - 0.25 * fx * fx * fz * fga * pb_y

                     + fx * fx * fx * fz * pb_y

                     + 5.0 * pa_xy * fx * fx * fz * pb_x

                     + 2.5 * fx * fx * pa_yy * fz * pb_y

                     + 5.0 * fx * fx * pa_yz * fz * pb_z

                     - 0.5 * pa_x * fz * fga * fx * pb_xy

                     - 0.5 * fx * fz * fga * pa_z * pb_yz

                     + 2.5 * pa_x * fx * fx * fz * pb_xy

                     + 2.5 * fx * fx * fz * pa_z * pb_yz

                     + 6.0 * pa_xyy * fz * fx * pb_xy

                     + 12.0 * pa_xyz * fx * fz * pb_xz

                     + 6.0 * fx * pa_yyz * fz * pb_yz

                     - pa_xz * fz * fga * pb_xyz

                     + 6.0 * pa_xz * fx * fz * pb_xyz

                     + 14.0 * pa_xyyz * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xyy,
                                    double pa_xyyz,
                                    double pa_xz,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xz,
                                    double pb_xzz,
                                    double pb_z,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx * pa_z

                     + 0.25 * fx * fx * fx * pb_z

                     + 0.25 * fx * fx * pa_yyz

                     + 0.5 * fx * fx * pa_yy * pb_z

                     + 0.25 * pa_xz * fx * fx * pb_x

                     + 0.5 * pa_x * fx * fx * pb_xz

                     + 0.25 * fx * fx * pa_z * pb_zz

                     + 0.5 * pa_xyyz * pb_x * fx

                     + pa_xyy * fx * pb_xz

                     + 0.5 * fx * pa_yyz * pb_zz

                     + 0.5 * pa_xz * fx * pb_xzz

                     + pa_xyyz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xyy,
                                    double pa_xyyz,
                                    double pa_xz,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xz,
                                    double pb_xzz,
                                    double pb_z,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.25 * fx * fx * pa_z * fz * fgb

                     - 0.25 * fx * fx * fz * fga * pa_z

                     - 0.5 * fx * fx * fz * fga * pb_z

                     - 0.5 * fx * pa_yyz * fz * fgb

                     + fx * fx * fx * fz * pa_z

                     + 2.0 * fx * fx * fx * fz * pb_z

                     + 2.5 * fx * fx * pa_yyz * fz

                     + 5.0 * fx * fx * pa_yy * fz * pb_z

                     - 0.5 * pa_xz * fx * pb_x * fz * fgb

                     - 0.5 * pa_xz * fz * fga * pb_x * fx

                     - pa_x * fz * fga * fx * pb_xz

                     - 0.5 * fx * fz * fga * pa_z * pb_zz

                     - pa_xyyz * pb_x * fz * fgb

                     + 2.5 * pa_xz * fx * fx * fz * pb_x

                     + 5.0 * pa_x * fx * fx * fz * pb_xz

                     + 2.5 * fx * fx * fz * pa_z * pb_zz

                     + 6.0 * pa_xyyz * fz * pb_x * fx

                     + 12.0 * pa_xyy * fz * fx * pb_xz

                     + 6.0 * fx * pa_yyz * fz * pb_zz

                     - pa_xz * fz * fga * pb_xzz

                     + 6.0 * pa_xz * fx * fz * pb_xzz

                     + 14.0 * pa_xyyz * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_yyy_s_0(double fx,
                                    double pa_xyyz,
                                    double pa_xyz,
                                    double pa_xz,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xyz * fx * fx

                     + 2.25 * pa_xz * fx * fx * pb_y

                     + 1.5 * pa_xyyz * pb_y * fx

                     + 3.0 * pa_xyz * fx * pb_yy

                     + 0.5 * pa_xz * fx * pb_yyy

                     + pa_xyyz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_yyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xyyz,
                                    double pa_xyz,
                                    double pa_xz,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xyz * fx * fz * fgb

                     + 15.0 * pa_xyz * fx * fx * fz

                     + 22.5 * pa_xz * fx * fx * fz * pb_y

                     - 1.5 * pa_xz * fx * pb_y * fz * fgb

                     - 1.5 * pa_xz * fz * fga * pb_y * fx

                     - 3.0 * pa_xyyz * pb_y * fz * fgb

                     + 18.0 * pa_xyyz * fz * pb_y * fx

                     + 36.0 * pa_xyz * fx * fz * pb_yy

                     - pa_xz * fz * fga * pb_yyy

                     + 6.0 * pa_xz * fx * fz * pb_yyy

                     + 14.0 * pa_xyyz * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_yyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pa_xyyz,
                                    double pa_xyz,
                                    double pa_xz,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 0.25 * pa_xyy * fx * fx

                     + pa_xy * fx * fx * pb_y

                     + 0.75 * pa_xz * fx * fx * pb_z

                     + 0.25 * pa_x * fx * fx * pb_yy

                     + 0.5 * pa_xyyz * fx * pb_z

                     + 0.5 * pa_xyy * fx * pb_yy

                     + 2.0 * pa_xyz * fx * pb_yz

                     + 0.5 * pa_xz * fx * pb_yyz

                     + pa_xyyz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_yyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pa_xyyz,
                                    double pa_xyz,
                                    double pa_xz,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (3.0 * pa_x * fx * fx * fx * fz

                     - 0.25 * pa_x * fx * fx * fz * fgb

                     - 0.25 * pa_x * fz * fga * fx * fx

                     - 0.5 * pa_xyy * fx * fz * fgb

                     + 2.5 * pa_xyy * fz * fx * fx

                     + 10.0 * pa_xy * fx * fx * fz * pb_y

                     + 7.5 * pa_xz * fx * fx * fz * pb_z

                     - 0.5 * pa_xz * fx * fz * fgb * pb_z

                     - 0.5 * pa_xz * fz * fga * fx * pb_z

                     - 0.5 * pa_x * fz * fga * fx * pb_yy

                     - pa_xyyz * fz * fgb * pb_z

                     + 2.5 * pa_x * fx * fx * fz * pb_yy

                     + 6.0 * pa_xyyz * fz * fx * pb_z

                     + 6.0 * pa_xyy * fz * fx * pb_yy

                     + 24.0 * pa_xyz * fx * fz * pb_yz

                     - pa_xz * fz * fga * pb_yyz

                     + 6.0 * pa_xz * fx * fz * pb_yyz

                     + 14.0 * pa_xyyz * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_yzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pa_xyyz,
                                    double pa_xyz,
                                    double pa_xz,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double pb_z,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xyz * fx * fx

                     + pa_xy * fx * fx * pb_z

                     + 0.25 * pa_xz * fx * fx * pb_y

                     + 0.5 * pa_x * fx * fx * pb_yz

                     + 0.5 * pa_xyyz * pb_y * fx

                     + pa_xyy * fx * pb_yz

                     + pa_xyz * fx * pb_zz

                     + 0.5 * pa_xz * fx * pb_yzz

                     + pa_xyyz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_yzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pa_xyyz,
                                    double pa_xyz,
                                    double pa_xz,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double pb_z,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- pa_xyz * fx * fz * fgb

                     + 5.0 * pa_xyz * fx * fx * fz

                     + 10.0 * pa_xy * fx * fx * fz * pb_z

                     - 0.5 * pa_xz * fx * pb_y * fz * fgb

                     - 0.5 * pa_xz * fz * fga * pb_y * fx

                     - pa_x * fz * fga * fx * pb_yz

                     - pa_xyyz * pb_y * fz * fgb

                     + 2.5 * pa_xz * fx * fx * fz * pb_y

                     + 5.0 * pa_x * fx * fx * fz * pb_yz

                     + 6.0 * pa_xyyz * fz * pb_y * fx

                     + 12.0 * pa_xyy * fz * fx * pb_yz

                     + 12.0 * pa_xyz * fx * fz * pb_zz

                     - pa_xz * fz * fga * pb_yzz

                     + 6.0 * pa_xz * fx * fz * pb_yzz

                     + 14.0 * pa_xyyz * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_zzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xyy,
                                    double pa_xyyz,
                                    double pa_xz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 0.75 * pa_xyy * fx * fx

                     + 0.75 * pa_xz * fx * fx * pb_z

                     + 0.75 * pa_x * fx * fx * pb_zz

                     + 1.5 * pa_xyyz * pb_z * fx

                     + 1.5 * pa_xyy * fx * pb_zz

                     + 0.5 * pa_xz * fx * pb_zzz

                     + pa_xyyz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_zzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xyy,
                                    double pa_xyyz,
                                    double pa_xz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * pa_x * fx * fx * fz * fgb

                     - 0.75 * pa_x * fz * fga * fx * fx

                     - 1.5 * pa_xyy * fx * fz * fgb

                     + 3.0 * pa_x * fx * fx * fx * fz

                     + 7.5 * pa_xyy * fz * fx * fx

                     - 1.5 * pa_xz * fx * pb_z * fz * fgb

                     - 1.5 * pa_xz * fz * fga * pb_z * fx

                     - 1.5 * pa_x * fz * fga * fx * pb_zz

                     - 3.0 * pa_xyyz * pb_z * fz * fgb

                     + 7.5 * pa_xz * fx * fx * fz * pb_z

                     + 7.5 * pa_x * fx * fx * fz * pb_zz

                     + 18.0 * pa_xyyz * fz * pb_z * fx

                     + 18.0 * pa_xyy * fz * fx * pb_zz

                     - pa_xz * fz * fga * pb_zzz

                     + 6.0 * pa_xz * fx * fz * pb_zzz

                     + 14.0 * pa_xyyz * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xxx_s_0(double fx,
                                    double pa_xy,
                                    double pa_xyzz,
                                    double pa_y,
                                    double pa_yzz,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_y

                     + 0.75 * fx * fx * pa_yzz

                     + 0.75 * pa_xy * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_y * pb_xx

                     + 1.5 * pa_xyzz * pb_x * fx

                     + 1.5 * fx * pa_yzz * pb_xx

                     + 0.5 * pa_xy * fx * pb_xxx

                     + pa_xyzz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xy,
                                    double pa_xyzz,
                                    double pa_y,
                                    double pa_yzz,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pa_y * fz * fgb

                     - 0.75 * fx * fx * pa_y * fz * fga

                     - 1.5 * fx * pa_yzz * fz * fgb

                     + 3.0 * fx * fx * fx * pa_y * fz

                     + 7.5 * fx * fx * pa_yzz * fz

                     - 1.5 * pa_xy * fx * pb_x * fz * fgb

                     - 1.5 * pa_xy * fz * fga * pb_x * fx

                     - 1.5 * fx * pa_y * fz * fga * pb_xx

                     - 3.0 * pa_xyzz * pb_x * fz * fgb

                     + 7.5 * pa_xy * fz * fx * fx * pb_x

                     + 7.5 * fx * fx * pa_y * fz * pb_xx

                     + 18.0 * pa_xyzz * fz * pb_x * fx

                     + 18.0 * fx * pa_yzz * fz * pb_xx

                     - pa_xy * fz * fga * pb_xxx

                     + 6.0 * pa_xy * fz * fx * pb_xxx

                     + 14.0 * pa_xyzz * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xxy_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyzz,
                                    double pa_xzz,
                                    double pa_y,
                                    double pa_yzz,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double s_0_0)
    {
        return s_0_0 * (0.125 * pa_x * fx * fx * fx

                     + 0.25 * fx * fx * fx * pb_x

                     + 0.25 * pa_xzz * fx * fx

                     + 0.5 * fx * fx * pa_zz * pb_x

                     + 0.25 * pa_xy * fx * fx * pb_y

                     + 0.25 * pa_x * fx * fx * pb_xx

                     + 0.5 * fx * fx * pa_y * pb_xy

                     + 0.5 * pa_xyzz * fx * pb_y

                     + 0.5 * pa_xzz * fx * pb_xx

                     + fx * pa_yzz * pb_xy

                     + 0.5 * pa_xy * fx * pb_xxy

                     + pa_xyzz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyzz,
                                    double pa_xzz,
                                    double pa_y,
                                    double pa_yzz,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.25 * pa_x * fx * fx * fz * fgb

                     - 0.25 * pa_x * fx * fx * fz * fga

                     - 0.5 * fx * fx * fz * fga * pb_x

                     - 0.5 * pa_xzz * fx * fz * fgb

                     + pa_x * fx * fx * fx * fz

                     + 2.0 * fx * fx * fx * fz * pb_x

                     + 2.5 * pa_xzz * fx * fx * fz

                     + 5.0 * fx * fx * pa_zz * fz * pb_x

                     - 0.5 * pa_xy * fx * fz * fgb * pb_y

                     - 0.5 * pa_xy * fz * fga * fx * pb_y

                     - 0.5 * pa_x * fx * fz * fga * pb_xx

                     - fx * pa_y * fz * fga * pb_xy

                     - pa_xyzz * fz * fgb * pb_y

                     + 2.5 * pa_xy * fz * fx * fx * pb_y

                     + 2.5 * pa_x * fx * fx * fz * pb_xx

                     + 5.0 * fx * fx * pa_y * fz * pb_xy

                     + 6.0 * pa_xyzz * fz * fx * pb_y

                     + 6.0 * pa_xzz * fx * fz * pb_xx

                     + 12.0 * fx * pa_yzz * fz * pb_xy

                     - pa_xy * fz * fga * pb_xxy

                     + 6.0 * pa_xy * fz * fx * pb_xxy

                     + 14.0 * pa_xyzz * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xxz_s_0(double fx,
                                    double pa_xy,
                                    double pa_xyz,
                                    double pa_xyzz,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xyz * fx * fx

                     + fx * fx * pa_yz * pb_x

                     + 0.25 * pa_xy * fx * fx * pb_z

                     + 0.5 * fx * fx * pa_y * pb_xz

                     + 0.5 * pa_xyzz * fx * pb_z

                     + pa_xyz * fx * pb_xx

                     + fx * pa_yzz * pb_xz

                     + 0.5 * pa_xy * fx * pb_xxz

                     + pa_xyzz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xy,
                                    double pa_xyz,
                                    double pa_xyzz,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- pa_xyz * fx * fz * fgb

                     + 5.0 * pa_xyz * fz * fx * fx

                     + 10.0 * fx * fx * pa_yz * fz * pb_x

                     - 0.5 * pa_xy * fx * fz * fgb * pb_z

                     - 0.5 * pa_xy * fz * fga * fx * pb_z

                     - fx * pa_y * fz * fga * pb_xz

                     - pa_xyzz * fz * fgb * pb_z

                     + 2.5 * pa_xy * fz * fx * fx * pb_z

                     + 5.0 * fx * fx * pa_y * fz * pb_xz

                     + 6.0 * pa_xyzz * fz * fx * pb_z

                     + 12.0 * pa_xyz * fz * fx * pb_xx

                     + 12.0 * fx * pa_yzz * fz * pb_xz

                     - pa_xy * fz * fga * pb_xxz

                     + 6.0 * pa_xy * fz * fx * pb_xxz

                     + 14.0 * pa_xyzz * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xyy_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyzz,
                                    double pa_xzz,
                                    double pa_y,
                                    double pa_yzz,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyy,
                                    double pb_y,
                                    double pb_yy,
                                    double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx * pa_y

                     + 0.25 * fx * fx * fx * pb_y

                     + 0.25 * fx * fx * pa_yzz

                     + 0.5 * fx * fx * pa_zz * pb_y

                     + 0.25 * pa_xy * fx * fx * pb_x

                     + 0.5 * pa_x * fx * fx * pb_xy

                     + 0.25 * fx * fx * pa_y * pb_yy

                     + 0.5 * pa_xyzz * pb_x * fx

                     + pa_xzz * fx * pb_xy

                     + 0.5 * fx * pa_yzz * pb_yy

                     + 0.5 * pa_xy * fx * pb_xyy

                     + pa_xyzz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyzz,
                                    double pa_xzz,
                                    double pa_y,
                                    double pa_yzz,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyy,
                                    double pb_y,
                                    double pb_yy,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.25 * fx * fx * pa_y * fz * fgb

                     - 0.25 * fx * fx * pa_y * fz * fga

                     - 0.5 * fx * fx * fz * fga * pb_y

                     - 0.5 * fx * pa_yzz * fz * fgb

                     + fx * fx * fx * pa_y * fz

                     + 2.0 * fx * fx * fx * fz * pb_y

                     + 2.5 * fx * fx * pa_yzz * fz

                     + 5.0 * fx * fx * pa_zz * fz * pb_y

                     - 0.5 * pa_xy * fx * pb_x * fz * fgb

                     - 0.5 * pa_xy * fz * fga * pb_x * fx

                     - pa_x * fx * fz * fga * pb_xy

                     - 0.5 * fx * pa_y * fz * fga * pb_yy

                     - pa_xyzz * pb_x * fz * fgb

                     + 2.5 * pa_xy * fz * fx * fx * pb_x

                     + 5.0 * pa_x * fx * fx * fz * pb_xy

                     + 2.5 * fx * fx * pa_y * fz * pb_yy

                     + 6.0 * pa_xyzz * fz * pb_x * fx

                     + 12.0 * pa_xzz * fx * fz * pb_xy

                     + 6.0 * fx * pa_yzz * fz * pb_yy

                     - pa_xy * fz * fga * pb_xyy

                     + 6.0 * pa_xy * fz * fx * pb_xyy

                     + 14.0 * pa_xyzz * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyz,
                                    double pa_xyzz,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyz,
                                    double pb_xz,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * fx * pa_z

                     + 0.125 * fx * fx * fx * pb_z

                     + 0.5 * pa_xz * fx * fx * pb_x

                     + 0.5 * fx * fx * pa_yz * pb_y

                     + 0.25 * fx * fx * pa_zz * pb_z

                     + 0.25 * pa_x * fx * fx * pb_xz

                     + 0.25 * fx * fx * pa_y * pb_yz

                     + pa_xyz * fx * pb_xy

                     + 0.5 * pa_xzz * fx * pb_xz

                     + 0.5 * fx * pa_yzz * pb_yz

                     + 0.5 * pa_xy * fx * pb_xyz

                     + pa_xyzz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xyz_r_0(double fga,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyz,
                                    double pa_xyzz,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyz,
                                    double pb_xz,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (2.0 * fx * fx * fx * pa_z * fz

                     - 0.25 * fx * fx * fz * fga * pb_z

                     + fx * fx * fx * fz * pb_z

                     + 5.0 * pa_xz * fx * fx * fz * pb_x

                     + 5.0 * fx * fx * pa_yz * fz * pb_y

                     + 2.5 * fx * fx * pa_zz * fz * pb_z

                     - 0.5 * pa_x * fx * fz * fga * pb_xz

                     - 0.5 * fx * pa_y * fz * fga * pb_yz

                     + 2.5 * pa_x * fx * fx * fz * pb_xz

                     + 2.5 * fx * fx * pa_y * fz * pb_yz

                     + 12.0 * pa_xyz * fz * fx * pb_xy

                     + 6.0 * pa_xzz * fx * fz * pb_xz

                     + 6.0 * fx * pa_yzz * fz * pb_yz

                     - pa_xy * fz * fga * pb_xyz

                     + 6.0 * pa_xy * fz * fx * pb_xyz

                     + 14.0 * pa_xyzz * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xzz_s_0(double fx,
                                    double pa_xy,
                                    double pa_xyz,
                                    double pa_xyzz,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pb_x,
                                    double pb_xz,
                                    double pb_xzz,
                                    double pb_z,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_y

                     + 0.75 * pa_xy * fx * fx * pb_x

                     + 0.25 * fx * fx * pa_yzz

                     + fx * fx * pa_yz * pb_z

                     + 0.25 * fx * fx * pa_y * pb_zz

                     + 0.5 * pa_xyzz * pb_x * fx

                     + 2.0 * pa_xyz * fx * pb_xz

                     + 0.5 * fx * pa_yzz * pb_zz

                     + 0.5 * pa_xy * fx * pb_xzz

                     + pa_xyzz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xy,
                                    double pa_xyz,
                                    double pa_xyzz,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pb_x,
                                    double pb_xz,
                                    double pb_xzz,
                                    double pb_z,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (3.0 * fx * fx * fx * pa_y * fz

                     - 0.25 * fx * fx * pa_y * fz * fgb

                     - 0.25 * fx * fx * pa_y * fz * fga

                     - 0.5 * fx * pa_yzz * fz * fgb

                     + 7.5 * pa_xy * fx * fx * fz * pb_x

                     + 2.5 * fx * fx * pa_yzz * fz

                     + 10.0 * fx * fx * pa_yz * fz * pb_z

                     - 0.5 * pa_xy * fx * pb_x * fz * fgb

                     - 0.5 * pa_xy * fz * fga * pb_x * fx

                     - 0.5 * fx * pa_y * fz * fga * pb_zz

                     - pa_xyzz * pb_x * fz * fgb

                     + 2.5 * fx * fx * pa_y * fz * pb_zz

                     + 6.0 * pa_xyzz * fz * pb_x * fx

                     + 24.0 * pa_xyz * fz * fx * pb_xz

                     + 6.0 * fx * pa_yzz * fz * pb_zz

                     - pa_xy * fz * fga * pb_xzz

                     + 6.0 * pa_xy * fz * fx * pb_xzz

                     + 14.0 * pa_xyzz * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_yyy_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyzz,
                                    double pa_xzz,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 0.75 * pa_xzz * fx * fx

                     + 0.75 * pa_xy * fx * fx * pb_y

                     + 0.75 * pa_x * fx * fx * pb_yy

                     + 1.5 * pa_xyzz * pb_y * fx

                     + 1.5 * pa_xzz * fx * pb_yy

                     + 0.5 * pa_xy * fx * pb_yyy

                     + pa_xyzz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_yyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyzz,
                                    double pa_xzz,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * pa_x * fx * fx * fz * fgb

                     - 0.75 * pa_x * fx * fx * fz * fga

                     - 1.5 * pa_xzz * fx * fz * fgb

                     + 3.0 * pa_x * fx * fx * fx * fz

                     + 7.5 * pa_xzz * fx * fx * fz

                     - 1.5 * pa_xy * fx * pb_y * fz * fgb

                     - 1.5 * pa_xy * fz * fga * pb_y * fx

                     - 1.5 * pa_x * fx * fz * fga * pb_yy

                     - 3.0 * pa_xyzz * pb_y * fz * fgb

                     + 7.5 * pa_xy * fz * fx * fx * pb_y

                     + 7.5 * pa_x * fx * fx * fz * pb_yy

                     + 18.0 * pa_xyzz * fz * pb_y * fx

                     + 18.0 * pa_xzz * fx * fz * pb_yy

                     - pa_xy * fz * fga * pb_yyy

                     + 6.0 * pa_xy * fz * fx * pb_yyy

                     + 14.0 * pa_xyzz * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_yyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyz,
                                    double pa_xyzz,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xyz * fx * fx

                     + pa_xz * fx * fx * pb_y

                     + 0.25 * pa_xy * fx * fx * pb_z

                     + 0.5 * pa_x * fx * fx * pb_yz

                     + 0.5 * pa_xyzz * fx * pb_z

                     + pa_xyz * fx * pb_yy

                     + pa_xzz * fx * pb_yz

                     + 0.5 * pa_xy * fx * pb_yyz

                     + pa_xyzz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_yyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyz,
                                    double pa_xyzz,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- pa_xyz * fx * fz * fgb

                     + 5.0 * pa_xyz * fz * fx * fx

                     + 10.0 * pa_xz * fx * fx * fz * pb_y

                     - 0.5 * pa_xy * fx * fz * fgb * pb_z

                     - 0.5 * pa_xy * fz * fga * fx * pb_z

                     - pa_x * fx * fz * fga * pb_yz

                     - pa_xyzz * fz * fgb * pb_z

                     + 2.5 * pa_xy * fz * fx * fx * pb_z

                     + 5.0 * pa_x * fx * fx * fz * pb_yz

                     + 6.0 * pa_xyzz * fz * fx * pb_z

                     + 12.0 * pa_xyz * fz * fx * pb_yy

                     + 12.0 * pa_xzz * fx * fz * pb_yz

                     - pa_xy * fz * fga * pb_yyz

                     + 6.0 * pa_xy * fz * fx * pb_yyz

                     + 14.0 * pa_xyzz * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_yzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyz,
                                    double pa_xyzz,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double pb_z,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 0.75 * pa_xy * fx * fx * pb_y

                     + 0.25 * pa_xzz * fx * fx

                     + pa_xz * fx * fx * pb_z

                     + 0.25 * pa_x * fx * fx * pb_zz

                     + 0.5 * pa_xyzz * pb_y * fx

                     + 2.0 * pa_xyz * fx * pb_yz

                     + 0.5 * pa_xzz * fx * pb_zz

                     + 0.5 * pa_xy * fx * pb_yzz

                     + pa_xyzz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_yzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyz,
                                    double pa_xyzz,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double pb_z,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (3.0 * pa_x * fx * fx * fx * fz

                     - 0.25 * pa_x * fx * fx * fz * fgb

                     - 0.25 * pa_x * fx * fx * fz * fga

                     - 0.5 * pa_xzz * fx * fz * fgb

                     + 7.5 * pa_xy * fx * fx * fz * pb_y

                     + 2.5 * pa_xzz * fx * fx * fz

                     + 10.0 * pa_xz * fx * fx * fz * pb_z

                     - 0.5 * pa_xy * fx * pb_y * fz * fgb

                     - 0.5 * pa_xy * fz * fga * pb_y * fx

                     - 0.5 * pa_x * fx * fz * fga * pb_zz

                     - pa_xyzz * pb_y * fz * fgb

                     + 2.5 * pa_x * fx * fx * fz * pb_zz

                     + 6.0 * pa_xyzz * fz * pb_y * fx

                     + 24.0 * pa_xyz * fz * fx * pb_yz

                     + 6.0 * pa_xzz * fx * fz * pb_zz

                     - pa_xy * fz * fga * pb_yzz

                     + 6.0 * pa_xy * fz * fx * pb_yzz

                     + 14.0 * pa_xyzz * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_zzz_s_0(double fx,
                                    double pa_xy,
                                    double pa_xyz,
                                    double pa_xyzz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xyz * fx * fx

                     + 2.25 * pa_xy * fx * fx * pb_z

                     + 1.5 * pa_xyzz * pb_z * fx

                     + 3.0 * pa_xyz * fx * pb_zz

                     + 0.5 * pa_xy * fx * pb_zzz

                     + pa_xyzz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_zzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xy,
                                    double pa_xyz,
                                    double pa_xyzz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xyz * fx * fz * fgb

                     + 15.0 * pa_xyz * fz * fx * fx

                     + 22.5 * pa_xy * fx * fx * fz * pb_z

                     - 1.5 * pa_xy * fx * pb_z * fz * fgb

                     - 1.5 * pa_xy * fz * fga * pb_z * fx

                     - 3.0 * pa_xyzz * pb_z * fz * fgb

                     + 18.0 * pa_xyzz * fz * pb_z * fx

                     + 36.0 * pa_xyz * fz * fx * pb_zz

                     - pa_xy * fz * fga * pb_zzz

                     + 6.0 * pa_xy * fz * fx * pb_zzz

                     + 14.0 * pa_xyzz * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xxx_s_0(double fx,
                                    double pa_xz,
                                    double pa_xzzz,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pa_z

                     + 0.75 * fx * fx * pa_zzz

                     + 2.25 * pa_xz * fx * fx * pb_x

                     + 2.25 * fx * fx * pa_z * pb_xx

                     + 1.5 * pa_xzzz * pb_x * fx

                     + 1.5 * fx * pa_zzz * pb_xx

                     + 1.5 * pa_xz * fx * pb_xxx

                     + pa_xzzz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xz,
                                    double pa_xzzz,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pa_z * fz * fgb

                     - 2.25 * fx * fx * pa_z * fz * fga

                     - 1.5 * fx * pa_zzz * fz * fgb

                     + 9.0 * fx * fx * fx * pa_z * fz

                     + 7.5 * fx * fx * pa_zzz * fz

                     - 4.5 * pa_xz * fx * pb_x * fz * fgb

                     - 4.5 * pa_xz * fz * fga * pb_x * fx

                     - 4.5 * fx * pa_z * fz * fga * pb_xx

                     - 3.0 * pa_xzzz * pb_x * fz * fgb

                     + 22.5 * pa_xz * fz * fx * fx * pb_x

                     + 22.5 * fx * fx * pa_z * fz * pb_xx

                     + 18.0 * pa_xzzz * fz * pb_x * fx

                     + 18.0 * fx * pa_zzz * fz * pb_xx

                     - 3.0 * pa_xz * fz * fga * pb_xxx

                     + 18.0 * pa_xz * fz * fx * pb_xxx

                     + 14.0 * pa_xzzz * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xxy_s_0(double fx,
                                    double pa_xz,
                                    double pa_xzzz,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xz * fx * fx * pb_y

                     + 1.5 * fx * fx * pa_z * pb_xy

                     + 0.5 * pa_xzzz * fx * pb_y

                     + fx * pa_zzz * pb_xy

                     + 1.5 * pa_xz * fx * pb_xxy

                     + pa_xzzz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xz,
                                    double pa_xzzz,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xz * fx * fz * fgb * pb_y

                     - 1.5 * pa_xz * fz * fga * fx * pb_y

                     - 3.0 * fx * pa_z * fz * fga * pb_xy

                     - pa_xzzz * fz * fgb * pb_y

                     + 7.5 * pa_xz * fz * fx * fx * pb_y

                     + 15.0 * fx * fx * pa_z * fz * pb_xy

                     + 6.0 * pa_xzzz * fz * fx * pb_y

                     + 12.0 * fx * pa_zzz * fz * pb_xy

                     - 3.0 * pa_xz * fz * fga * pb_xxy

                     + 18.0 * pa_xz * fz * fx * pb_xxy

                     + 14.0 * pa_xzzz * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xxz_s_0(double fx,
                                    double pa_x,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pa_xzzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 0.75 * fx * fx * fx * pb_x

                     + 0.75 * pa_xzz * fx * fx

                     + 1.5 * fx * fx * pa_zz * pb_x

                     + 0.75 * pa_xz * fx * fx * pb_z

                     + 0.75 * pa_x * fx * fx * pb_xx

                     + 1.5 * fx * fx * pa_z * pb_xz

                     + 0.5 * pa_xzzz * fx * pb_z

                     + 1.5 * pa_xzz * fx * pb_xx

                     + fx * pa_zzz * pb_xz

                     + 1.5 * pa_xz * fx * pb_xxz

                     + pa_xzzz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pa_xzzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * pa_x * fx * fx * fz * fgb

                     - 0.75 * pa_x * fx * fx * fz * fga

                     - 1.5 * fx * fx * fz * fga * pb_x

                     - 1.5 * pa_xzz * fx * fz * fgb

                     + 3.0 * pa_x * fx * fx * fx * fz

                     + 6.0 * fx * fx * fx * fz * pb_x

                     + 7.5 * pa_xzz * fz * fx * fx

                     + 15.0 * fx * fx * pa_zz * fz * pb_x

                     - 1.5 * pa_xz * fx * fz * fgb * pb_z

                     - 1.5 * pa_xz * fz * fga * fx * pb_z

                     - 1.5 * pa_x * fx * fz * fga * pb_xx

                     - 3.0 * fx * pa_z * fz * fga * pb_xz

                     - pa_xzzz * fz * fgb * pb_z

                     + 7.5 * pa_xz * fz * fx * fx * pb_z

                     + 7.5 * pa_x * fx * fx * fz * pb_xx

                     + 15.0 * fx * fx * pa_z * fz * pb_xz

                     + 6.0 * pa_xzzz * fz * fx * pb_z

                     + 18.0 * pa_xzz * fz * fx * pb_xx

                     + 12.0 * fx * pa_zzz * fz * pb_xz

                     - 3.0 * pa_xz * fz * fga * pb_xxz

                     + 18.0 * pa_xz * fz * fx * pb_xxz

                     + 14.0 * pa_xzzz * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xyy_s_0(double fx,
                                    double pa_xz,
                                    double pa_xzzz,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_x,
                                    double pb_xyy,
                                    double pb_yy,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_z

                     + 0.25 * fx * fx * pa_zzz

                     + 0.75 * pa_xz * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_z * pb_yy

                     + 0.5 * pa_xzzz * pb_x * fx

                     + 0.5 * fx * pa_zzz * pb_yy

                     + 1.5 * pa_xz * fx * pb_xyy

                     + pa_xzzz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xz,
                                    double pa_xzzz,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_x,
                                    double pb_xyy,
                                    double pb_yy,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pa_z * fz * fgb

                     - 0.75 * fx * fx * pa_z * fz * fga

                     - 0.5 * fx * pa_zzz * fz * fgb

                     + 3.0 * fx * fx * fx * pa_z * fz

                     + 2.5 * fx * fx * pa_zzz * fz

                     - 1.5 * pa_xz * fx * pb_x * fz * fgb

                     - 1.5 * pa_xz * fz * fga * pb_x * fx

                     - 1.5 * fx * pa_z * fz * fga * pb_yy

                     - pa_xzzz * pb_x * fz * fgb

                     + 7.5 * pa_xz * fz * fx * fx * pb_x

                     + 7.5 * fx * fx * pa_z * fz * pb_yy

                     + 6.0 * pa_xzzz * fz * pb_x * fx

                     + 6.0 * fx * pa_zzz * fz * pb_yy

                     - 3.0 * pa_xz * fz * fga * pb_xyy

                     + 18.0 * pa_xz * fz * fx * pb_xyy

                     + 14.0 * pa_xzzz * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pa_xzzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_xy,
                                    double pb_xyz,
                                    double pb_y,
                                    double pb_yz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_zz * pb_y

                     + 0.75 * pa_x * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_z * pb_yz

                     + 1.5 * pa_xzz * fx * pb_xy

                     + 0.5 * fx * pa_zzz * pb_yz

                     + 1.5 * pa_xz * fx * pb_xyz

                     + pa_xzzz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xyz_r_0(double fga,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pa_xzzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_xy,
                                    double pb_xyz,
                                    double pb_y,
                                    double pb_yz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fga * pb_y

                     + 3.0 * fx * fx * fx * fz * pb_y

                     + 7.5 * fx * fx * pa_zz * fz * pb_y

                     - 1.5 * pa_x * fx * fz * fga * pb_xy

                     - 1.5 * fx * pa_z * fz * fga * pb_yz

                     + 7.5 * pa_x * fx * fx * fz * pb_xy

                     + 7.5 * fx * fx * pa_z * fz * pb_yz

                     + 18.0 * pa_xzz * fz * fx * pb_xy

                     + 6.0 * fx * pa_zzz * fz * pb_yz

                     - 3.0 * pa_xz * fz * fga * pb_xyz

                     + 18.0 * pa_xz * fz * fx * pb_xyz

                     + 14.0 * pa_xzzz * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pa_xzzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_x,
                                    double pb_xz,
                                    double pb_xzz,
                                    double pb_z,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pa_z

                     + 0.75 * fx * fx * fx * pb_z

                     + 2.25 * pa_xz * fx * fx * pb_x

                     + 0.25 * fx * fx * pa_zzz

                     + 1.5 * fx * fx * pa_zz * pb_z

                     + 1.5 * pa_x * fx * fx * pb_xz

                     + 0.75 * fx * fx * pa_z * pb_zz

                     + 0.5 * pa_xzzz * pb_x * fx

                     + 3.0 * pa_xzz * fx * pb_xz

                     + 0.5 * fx * pa_zzz * pb_zz

                     + 1.5 * pa_xz * fx * pb_xzz

                     + pa_xzzz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pa_xzzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_x,
                                    double pb_xz,
                                    double pb_xzz,
                                    double pb_z,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (9.0 * fx * fx * fx * pa_z * fz

                     - 0.75 * fx * fx * pa_z * fz * fgb

                     - 0.75 * fx * fx * pa_z * fz * fga

                     - 1.5 * fx * fx * fz * fga * pb_z

                     - 0.5 * fx * pa_zzz * fz * fgb

                     + 6.0 * fx * fx * fx * fz * pb_z

                     + 22.5 * pa_xz * fx * fx * fz * pb_x

                     + 2.5 * fx * fx * pa_zzz * fz

                     + 15.0 * fx * fx * pa_zz * fz * pb_z

                     - 1.5 * pa_xz * fx * pb_x * fz * fgb

                     - 1.5 * pa_xz * fz * fga * pb_x * fx

                     - 3.0 * pa_x * fx * fz * fga * pb_xz

                     - 1.5 * fx * pa_z * fz * fga * pb_zz

                     - pa_xzzz * pb_x * fz * fgb

                     + 15.0 * pa_x * fx * fx * fz * pb_xz

                     + 7.5 * fx * fx * pa_z * fz * pb_zz

                     + 6.0 * pa_xzzz * fz * pb_x * fx

                     + 36.0 * pa_xzz * fz * fx * pb_xz

                     + 6.0 * fx * pa_zzz * fz * pb_zz

                     - 3.0 * pa_xz * fz * fga * pb_xzz

                     + 18.0 * pa_xz * fz * fx * pb_xzz

                     + 14.0 * pa_xzzz * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_yyy_s_0(double fx,
                                    double pa_xz,
                                    double pa_xzzz,
                                    double pb_y,
                                    double pb_yyy,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_xz * fx * fx * pb_y

                     + 1.5 * pa_xzzz * pb_y * fx

                     + 1.5 * pa_xz * fx * pb_yyy

                     + pa_xzzz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_yyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xz,
                                    double pa_xzzz,
                                    double pb_y,
                                    double pb_yyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_xz * fx * pb_y * fz * fgb

                     - 4.5 * pa_xz * fz * fga * pb_y * fx

                     - 3.0 * pa_xzzz * pb_y * fz * fgb

                     + 22.5 * pa_xz * fz * fx * fx * pb_y

                     + 18.0 * pa_xzzz * fz * pb_y * fx

                     - 3.0 * pa_xz * fz * fga * pb_yyy

                     + 18.0 * pa_xz * fz * fx * pb_yyy

                     + 14.0 * pa_xzzz * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_yyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pa_xzzz,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 0.75 * pa_xzz * fx * fx

                     + 0.75 * pa_xz * fx * fx * pb_z

                     + 0.75 * pa_x * fx * fx * pb_yy

                     + 0.5 * pa_xzzz * fx * pb_z

                     + 1.5 * pa_xzz * fx * pb_yy

                     + 1.5 * pa_xz * fx * pb_yyz

                     + pa_xzzz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_yyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pa_xzzz,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * pa_x * fx * fx * fz * fgb

                     - 0.75 * pa_x * fx * fx * fz * fga

                     - 1.5 * pa_xzz * fx * fz * fgb

                     + 3.0 * pa_x * fx * fx * fx * fz

                     + 7.5 * pa_xzz * fz * fx * fx

                     - 1.5 * pa_xz * fx * fz * fgb * pb_z

                     - 1.5 * pa_xz * fz * fga * fx * pb_z

                     - 1.5 * pa_x * fx * fz * fga * pb_yy

                     - pa_xzzz * fz * fgb * pb_z

                     + 7.5 * pa_xz * fz * fx * fx * pb_z

                     + 7.5 * pa_x * fx * fx * fz * pb_yy

                     + 6.0 * pa_xzzz * fz * fx * pb_z

                     + 18.0 * pa_xzz * fz * fx * pb_yy

                     - 3.0 * pa_xz * fz * fga * pb_yyz

                     + 18.0 * pa_xz * fz * fx * pb_yyz

                     + 14.0 * pa_xzzz * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_yzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pa_xzzz,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_xz * fx * fx * pb_y

                     + 1.5 * pa_x * fx * fx * pb_yz

                     + 0.5 * pa_xzzz * pb_y * fx

                     + 3.0 * pa_xzz * fx * pb_yz

                     + 1.5 * pa_xz * fx * pb_yzz

                     + pa_xzzz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_yzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pa_xzzz,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double r_0_0)
    {
        return r_0_0 * (22.5 * pa_xz * fx * fx * fz * pb_y

                     - 1.5 * pa_xz * fx * pb_y * fz * fgb

                     - 1.5 * pa_xz * fz * fga * pb_y * fx

                     - 3.0 * pa_x * fx * fz * fga * pb_yz

                     - pa_xzzz * pb_y * fz * fgb

                     + 15.0 * pa_x * fx * fx * fz * pb_yz

                     + 6.0 * pa_xzzz * fz * pb_y * fx

                     + 36.0 * pa_xzz * fz * fx * pb_yz

                     - 3.0 * pa_xz * fz * fga * pb_yzz

                     + 18.0 * pa_xz * fz * fx * pb_yzz

                     + 14.0 * pa_xzzz * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_zzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pa_xzzz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * pa_x * fx * fx * fx

                     + 2.25 * pa_xzz * fx * fx

                     + 6.75 * pa_xz * fx * fx * pb_z

                     + 2.25 * pa_x * fx * fx * pb_zz

                     + 1.5 * pa_xzzz * pb_z * fx

                     + 4.5 * pa_xzz * fx * pb_zz

                     + 1.5 * pa_xz * fx * pb_zzz

                     + pa_xzzz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_zzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pa_xzzz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (15.0 * pa_x * fx * fx * fx * fz

                     - 2.25 * pa_x * fx * fx * fz * fgb

                     - 2.25 * pa_x * fx * fx * fz * fga

                     - 4.5 * pa_xzz * fx * fz * fgb

                     + 22.5 * pa_xzz * fz * fx * fx

                     + 67.5 * pa_xz * fx * fx * fz * pb_z

                     - 4.5 * pa_xz * fx * pb_z * fz * fgb

                     - 4.5 * pa_xz * fz * fga * pb_z * fx

                     - 4.5 * pa_x * fx * fz * fga * pb_zz

                     - 3.0 * pa_xzzz * pb_z * fz * fgb

                     + 22.5 * pa_x * fx * fx * fz * pb_zz

                     + 18.0 * pa_xzzz * fz * pb_z * fx

                     + 54.0 * pa_xzz * fz * fx * pb_zz

                     - 3.0 * pa_xz * fz * fga * pb_zzz

                     + 18.0 * pa_xz * fz * fx * pb_zzz

                     + 14.0 * pa_xzzz * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xxx_s_0(double fx,
                                    double pa_yy,
                                    double pa_yyyy,
                                    double pb_x,
                                    double pb_xxx,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pb_x

                     + 4.5 * pa_yy * fx * fx * pb_x

                     + 1.5 * pa_yyyy * pb_x * fx

                     + 0.75 * fx * fx * pb_xxx

                     + 3.0 * pa_yy * fx * pb_xxx

                     + pa_yyyy * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yy,
                                    double pa_yyyy,
                                    double pb_x,
                                    double pb_xxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pb_x * fz * fgb

                     - 4.5 * fx * fx * fz * fga * pb_x

                     - 9.0 * pa_yy * fx * pb_x * fz * fgb

                     - 9.0 * pa_yy * fz * fga * pb_x * fx

                     + 9.0 * fx * fx * fx * fz * pb_x

                     - 3.0 * pa_yyyy * pb_x * fz * fgb

                     + 45.0 * pa_yy * fz * fx * fx * pb_x

                     + 18.0 * pa_yyyy * fz * pb_x * fx

                     - 3.0 * fx * fz * fga * pb_xxx

                     - 6.0 * pa_yy * fz * fga * pb_xxx

                     + 7.5 * fx * fx * fz * pb_xxx

                     + 36.0 * pa_yy * fz * fx * pb_xxx

                     + 14.0 * pa_yyyy * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xxy_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pa_yyyy,
                                    double pb_xx,
                                    double pb_xxy,
                                    double pb_y,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_y * fx * fx * fx

                     + pa_yyy * fx * fx

                     + 0.375 * fx * fx * fx * pb_y

                     + 1.5 * pa_yy * fx * fx * pb_y

                     + 3.0 * pa_y * fx * fx * pb_xx

                     + 0.5 * pa_yyyy * fx * pb_y

                     + 2.0 * pa_yyy * fx * pb_xx

                     + 0.75 * fx * fx * pb_xxy

                     + 3.0 * pa_yy * fx * pb_xxy

                     + pa_yyyy * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pa_yyyy,
                                    double pb_xx,
                                    double pb_xxy,
                                    double pb_y,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_y * fx * fx * fz * fgb

                     - 3.0 * pa_y * fx * fx * fz * fga

                     - 2.0 * pa_yyy * fx * fz * fgb

                     + 12.0 * pa_y * fx * fx * fx * fz

                     + 10.0 * pa_yyy * fz * fx * fx

                     - 0.75 * fx * fx * fz * fgb * pb_y

                     - 1.5 * fx * fx * fz * fga * pb_y

                     - 3.0 * pa_yy * fx * fz * fgb * pb_y

                     - 3.0 * pa_yy * fz * fga * fx * pb_y

                     - 6.0 * pa_y * fx * fz * fga * pb_xx

                     + 3.0 * fx * fx * fx * fz * pb_y

                     - pa_yyyy * fz * fgb * pb_y

                     + 15.0 * pa_yy * fz * fx * fx * pb_y

                     + 30.0 * pa_y * fx * fx * fz * pb_xx

                     + 6.0 * pa_yyyy * fz * fx * pb_y

                     + 24.0 * pa_yyy * fz * fx * pb_xx

                     - 3.0 * fx * fz * fga * pb_xxy

                     - 6.0 * pa_yy * fz * fga * pb_xxy

                     + 7.5 * fx * fx * fz * pb_xxy

                     + 36.0 * pa_yy * fz * fx * pb_xxy

                     + 14.0 * pa_yyyy * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xxz_s_0(double fx,
                                    double pa_yy,
                                    double pa_yyyy,
                                    double pb_xxz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_z

                     + 1.5 * pa_yy * fx * fx * pb_z

                     + 0.5 * pa_yyyy * fx * pb_z

                     + 0.75 * fx * fx * pb_xxz

                     + 3.0 * pa_yy * fx * pb_xxz

                     + pa_yyyy * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yy,
                                    double pa_yyyy,
                                    double pb_xxz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb * pb_z

                     - 1.5 * fx * fx * fz * fga * pb_z

                     - 3.0 * pa_yy * fx * fz * fgb * pb_z

                     - 3.0 * pa_yy * fz * fga * fx * pb_z

                     + 3.0 * fx * fx * fx * fz * pb_z

                     - pa_yyyy * fz * fgb * pb_z

                     + 15.0 * pa_yy * fz * fx * fx * pb_z

                     + 6.0 * pa_yyyy * fz * fx * pb_z

                     - 3.0 * fx * fz * fga * pb_xxz

                     - 6.0 * pa_yy * fz * fga * pb_xxz

                     + 7.5 * fx * fx * fz * pb_xxz

                     + 36.0 * pa_yy * fz * fx * pb_xxz

                     + 14.0 * pa_yyyy * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xyy_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pa_yyyy,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyy,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pb_x

                     + 4.5 * pa_yy * fx * fx * pb_x

                     + 6.0 * pa_y * fx * fx * pb_xy

                     + 0.5 * pa_yyyy * pb_x * fx

                     + 4.0 * pa_yyy * fx * pb_xy

                     + 0.75 * fx * fx * pb_xyy

                     + 3.0 * pa_yy * fx * pb_xyy

                     + pa_yyyy * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pa_yyyy,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * fx * fx * fz * fga * pb_x

                     + 15.0 * fx * fx * fx * fz * pb_x

                     + 45.0 * pa_yy * fx * fx * fz * pb_x

                     - 0.75 * fx * fx * pb_x * fz * fgb

                     - 3.0 * pa_yy * fx * pb_x * fz * fgb

                     - 3.0 * pa_yy * fz * fga * pb_x * fx

                     - 12.0 * pa_y * fx * fz * fga * pb_xy

                     - pa_yyyy * pb_x * fz * fgb

                     + 60.0 * pa_y * fx * fx * fz * pb_xy

                     + 6.0 * pa_yyyy * fz * pb_x * fx

                     + 48.0 * pa_yyy * fz * fx * pb_xy

                     - 3.0 * fx * fz * fga * pb_xyy

                     - 6.0 * pa_yy * fz * fga * pb_xyy

                     + 7.5 * fx * fx * fz * pb_xyy

                     + 36.0 * pa_yy * fz * fx * pb_xyy

                     + 14.0 * pa_yyyy * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xyz_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pa_yyyy,
                                    double pb_xyz,
                                    double pb_xz,
                                    double s_0_0)
    {
        return s_0_0 * (3.0 * pa_y * fx * fx * pb_xz

                     + 2.0 * pa_yyy * fx * pb_xz

                     + 0.75 * fx * fx * pb_xyz

                     + 3.0 * pa_yy * fx * pb_xyz

                     + pa_yyyy * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xyz_r_0(double fga,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pa_yyyy,
                                    double pb_xyz,
                                    double pb_xz,
                                    double r_0_0)
    {
        return r_0_0 * (- 6.0 * pa_y * fx * fz * fga * pb_xz

                     + 30.0 * pa_y * fx * fx * fz * pb_xz

                     + 24.0 * pa_yyy * fz * fx * pb_xz

                     - 3.0 * fx * fz * fga * pb_xyz

                     - 6.0 * pa_yy * fz * fga * pb_xyz

                     + 7.5 * fx * fx * fz * pb_xyz

                     + 36.0 * pa_yy * fz * fx * pb_xyz

                     + 14.0 * pa_yyyy * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xzz_s_0(double fx,
                                    double pa_yy,
                                    double pa_yyyy,
                                    double pb_x,
                                    double pb_xzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_x

                     + 1.5 * pa_yy * fx * fx * pb_x

                     + 0.5 * pa_yyyy * pb_x * fx

                     + 0.75 * fx * fx * pb_xzz

                     + 3.0 * pa_yy * fx * pb_xzz

                     + pa_yyyy * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yy,
                                    double pa_yyyy,
                                    double pb_x,
                                    double pb_xzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pb_x * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_x

                     - 3.0 * pa_yy * fx * pb_x * fz * fgb

                     - 3.0 * pa_yy * fz * fga * pb_x * fx

                     + 3.0 * fx * fx * fx * fz * pb_x

                     - pa_yyyy * pb_x * fz * fgb

                     + 15.0 * pa_yy * fz * fx * fx * pb_x

                     + 6.0 * pa_yyyy * fz * pb_x * fx

                     - 3.0 * fx * fz * fga * pb_xzz

                     - 6.0 * pa_yy * fz * fga * pb_xzz

                     + 7.5 * fx * fx * fz * pb_xzz

                     + 36.0 * pa_yy * fz * fx * pb_xzz

                     + 14.0 * pa_yyyy * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_yyy_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pa_yyyy,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double s_0_0)
    {
        return s_0_0 * (7.5 * pa_y * fx * fx * fx

                     + 5.625 * fx * fx * fx * pb_y

                     + 3.0 * pa_yyy * fx * fx

                     + 13.5 * pa_yy * fx * fx * pb_y

                     + 9.0 * pa_y * fx * fx * pb_yy

                     + 1.5 * pa_yyyy * pb_y * fx

                     + 6.0 * pa_yyy * fx * pb_yy

                     + 0.75 * fx * fx * pb_yyy

                     + 3.0 * pa_yy * fx * pb_yyy

                     + pa_yyyy * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_yyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pa_yyyy,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double r_0_0)
    {
        return r_0_0 * (60.0 * pa_y * fx * fx * fx * fz

                     - 9.0 * pa_y * fx * fx * fz * fgb

                     - 9.0 * pa_y * fx * fx * fz * fga

                     - 13.5 * fx * fx * fz * fga * pb_y

                     - 6.0 * pa_yyy * fx * fz * fgb

                     + 45.0 * fx * fx * fx * fz * pb_y

                     + 30.0 * pa_yyy * fz * fx * fx

                     + 135.0 * pa_yy * fx * fx * fz * pb_y

                     - 2.25 * fx * fx * pb_y * fz * fgb

                     - 9.0 * pa_yy * fx * pb_y * fz * fgb

                     - 9.0 * pa_yy * fz * fga * pb_y * fx

                     - 18.0 * pa_y * fx * fz * fga * pb_yy

                     - 3.0 * pa_yyyy * pb_y * fz * fgb

                     + 90.0 * pa_y * fx * fx * fz * pb_yy

                     + 18.0 * pa_yyyy * fz * pb_y * fx

                     + 72.0 * pa_yyy * fz * fx * pb_yy

                     - 3.0 * fx * fz * fga * pb_yyy

                     - 6.0 * pa_yy * fz * fga * pb_yyy

                     + 7.5 * fx * fx * fz * pb_yyy

                     + 36.0 * pa_yy * fz * fx * pb_yyy

                     + 14.0 * pa_yyyy * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_yyz_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pa_yyyy,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pb_z

                     + 4.5 * pa_yy * fx * fx * pb_z

                     + 6.0 * pa_y * fx * fx * pb_yz

                     + 0.5 * pa_yyyy * fx * pb_z

                     + 4.0 * pa_yyy * fx * pb_yz

                     + 0.75 * fx * fx * pb_yyz

                     + 3.0 * pa_yy * fx * pb_yyz

                     + pa_yyyy * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_yyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pa_yyyy,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * fx * fx * fz * fga * pb_z

                     + 15.0 * fx * fx * fx * fz * pb_z

                     + 45.0 * pa_yy * fx * fx * fz * pb_z

                     - 0.75 * fx * fx * fz * fgb * pb_z

                     - 3.0 * pa_yy * fx * fz * fgb * pb_z

                     - 3.0 * pa_yy * fz * fga * fx * pb_z

                     - 12.0 * pa_y * fx * fz * fga * pb_yz

                     - pa_yyyy * fz * fgb * pb_z

                     + 60.0 * pa_y * fx * fx * fz * pb_yz

                     + 6.0 * pa_yyyy * fz * fx * pb_z

                     + 48.0 * pa_yyy * fz * fx * pb_yz

                     - 3.0 * fx * fz * fga * pb_yyz

                     - 6.0 * pa_yy * fz * fga * pb_yyz

                     + 7.5 * fx * fx * fz * pb_yyz

                     + 36.0 * pa_yy * fz * fx * pb_yyz

                     + 14.0 * pa_yyyy * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_yzz_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pa_yyyy,
                                    double pb_y,
                                    double pb_yzz,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_y * fx * fx * fx

                     + pa_yyy * fx * fx

                     + 0.375 * fx * fx * fx * pb_y

                     + 1.5 * pa_yy * fx * fx * pb_y

                     + 3.0 * pa_y * fx * fx * pb_zz

                     + 0.5 * pa_yyyy * pb_y * fx

                     + 2.0 * pa_yyy * fx * pb_zz

                     + 0.75 * fx * fx * pb_yzz

                     + 3.0 * pa_yy * fx * pb_yzz

                     + pa_yyyy * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_yzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pa_yyyy,
                                    double pb_y,
                                    double pb_yzz,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_y * fx * fx * fz * fgb

                     - 3.0 * pa_y * fx * fx * fz * fga

                     - 2.0 * pa_yyy * fx * fz * fgb

                     + 12.0 * pa_y * fx * fx * fx * fz

                     + 10.0 * pa_yyy * fz * fx * fx

                     - 0.75 * fx * fx * pb_y * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_y

                     - 3.0 * pa_yy * fx * pb_y * fz * fgb

                     - 3.0 * pa_yy * fz * fga * pb_y * fx

                     - 6.0 * pa_y * fx * fz * fga * pb_zz

                     + 3.0 * fx * fx * fx * fz * pb_y

                     - pa_yyyy * pb_y * fz * fgb

                     + 15.0 * pa_yy * fz * fx * fx * pb_y

                     + 30.0 * pa_y * fx * fx * fz * pb_zz

                     + 6.0 * pa_yyyy * fz * pb_y * fx

                     + 24.0 * pa_yyy * fz * fx * pb_zz

                     - 3.0 * fx * fz * fga * pb_yzz

                     - 6.0 * pa_yy * fz * fga * pb_yzz

                     + 7.5 * fx * fx * fz * pb_yzz

                     + 36.0 * pa_yy * fz * fx * pb_yzz

                     + 14.0 * pa_yyyy * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_zzz_s_0(double fx,
                                    double pa_yy,
                                    double pa_yyyy,
                                    double pb_z,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pb_z

                     + 4.5 * pa_yy * fx * fx * pb_z

                     + 1.5 * pa_yyyy * pb_z * fx

                     + 0.75 * fx * fx * pb_zzz

                     + 3.0 * pa_yy * fx * pb_zzz

                     + pa_yyyy * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_zzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yy,
                                    double pa_yyyy,
                                    double pb_z,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pb_z * fz * fgb

                     - 4.5 * fx * fx * fz * fga * pb_z

                     - 9.0 * pa_yy * fx * pb_z * fz * fgb

                     - 9.0 * pa_yy * fz * fga * pb_z * fx

                     + 9.0 * fx * fx * fx * fz * pb_z

                     - 3.0 * pa_yyyy * pb_z * fz * fgb

                     + 45.0 * pa_yy * fz * fx * fx * pb_z

                     + 18.0 * pa_yyyy * fz * pb_z * fx

                     - 3.0 * fx * fz * fga * pb_zzz

                     - 6.0 * pa_yy * fz * fga * pb_zzz

                     + 7.5 * fx * fx * fz * pb_zzz

                     + 36.0 * pa_yy * fz * fx * pb_zzz

                     + 14.0 * pa_yyyy * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xxx_s_0(double fx,
                                    double pa_yyyz,
                                    double pa_yz,
                                    double pb_x,
                                    double pb_xxx,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_yz * fx * fx * pb_x

                     + 1.5 * pa_yyyz * pb_x * fx

                     + 1.5 * pa_yz * fx * pb_xxx

                     + pa_yyyz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yyyz,
                                    double pa_yz,
                                    double pb_x,
                                    double pb_xxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_yz * fx * pb_x * fz * fgb

                     - 4.5 * pa_yz * fz * fga * pb_x * fx

                     - 3.0 * pa_yyyz * pb_x * fz * fgb

                     + 22.5 * pa_yz * fx * fx * fz * pb_x

                     + 18.0 * pa_yyyz * fz * pb_x * fx

                     - 3.0 * pa_yz * fz * fga * pb_xxx

                     + 18.0 * pa_yz * fx * fz * pb_xxx

                     + 14.0 * pa_yyyz * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xxy_s_0(double fx,
                                    double pa_yyyz,
                                    double pa_yyz,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_xx,
                                    double pb_xxy,
                                    double pb_y,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_z

                     + 0.75 * pa_yyz * fx * fx

                     + 0.75 * pa_yz * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_z * pb_xx

                     + 0.5 * pa_yyyz * fx * pb_y

                     + 1.5 * pa_yyz * fx * pb_xx

                     + 1.5 * pa_yz * fx * pb_xxy

                     + pa_yyyz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yyyz,
                                    double pa_yyz,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_xx,
                                    double pb_xxy,
                                    double pb_y,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pa_z * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pa_z

                     - 1.5 * pa_yyz * fx * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_z

                     + 7.5 * pa_yyz * fx * fx * fz

                     - 1.5 * pa_yz * fx * fz * fgb * pb_y

                     - 1.5 * pa_yz * fz * fga * fx * pb_y

                     - 1.5 * fx * fz * fga * pa_z * pb_xx

                     - pa_yyyz * fz * fgb * pb_y

                     + 7.5 * pa_yz * fx * fx * fz * pb_y

                     + 7.5 * fx * fx * fz * pa_z * pb_xx

                     + 6.0 * pa_yyyz * fz * fx * pb_y

                     + 18.0 * pa_yyz * fx * fz * pb_xx

                     - 3.0 * pa_yz * fz * fga * pb_xxy

                     + 18.0 * pa_yz * fx * fz * pb_xxy

                     + 14.0 * pa_yyyz * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xxz_s_0(double fx,
                                    double pa_y,
                                    double pa_yyy,
                                    double pa_yyyz,
                                    double pa_yz,
                                    double pb_xx,
                                    double pb_xxz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_y * fx * fx * fx

                     + 0.25 * pa_yyy * fx * fx

                     + 0.75 * pa_yz * fx * fx * pb_z

                     + 0.75 * pa_y * fx * fx * pb_xx

                     + 0.5 * pa_yyyz * fx * pb_z

                     + 0.5 * pa_yyy * fx * pb_xx

                     + 1.5 * pa_yz * fx * pb_xxz

                     + pa_yyyz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yyy,
                                    double pa_yyyz,
                                    double pa_yz,
                                    double pb_xx,
                                    double pb_xxz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * pa_y * fx * fx * fz * fgb

                     - 0.75 * pa_y * fz * fga * fx * fx

                     - 0.5 * pa_yyy * fx * fz * fgb

                     + 3.0 * pa_y * fx * fx * fx * fz

                     + 2.5 * pa_yyy * fz * fx * fx

                     - 1.5 * pa_yz * fx * fz * fgb * pb_z

                     - 1.5 * pa_yz * fz * fga * fx * pb_z

                     - 1.5 * pa_y * fz * fga * fx * pb_xx

                     - pa_yyyz * fz * fgb * pb_z

                     + 7.5 * pa_yz * fx * fx * fz * pb_z

                     + 7.5 * pa_y * fx * fx * fz * pb_xx

                     + 6.0 * pa_yyyz * fz * fx * pb_z

                     + 6.0 * pa_yyy * fz * fx * pb_xx

                     - 3.0 * pa_yz * fz * fga * pb_xxz

                     + 18.0 * pa_yz * fx * fz * pb_xxz

                     + 14.0 * pa_yyyz * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xyy_s_0(double fx,
                                    double pa_yyyz,
                                    double pa_yyz,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyy,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_yz * fx * fx * pb_x

                     + 1.5 * fx * fx * pa_z * pb_xy

                     + 0.5 * pa_yyyz * pb_x * fx

                     + 3.0 * pa_yyz * fx * pb_xy

                     + 1.5 * pa_yz * fx * pb_xyy

                     + pa_yyyz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yyyz,
                                    double pa_yyz,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyy,
                                    double r_0_0)
    {
        return r_0_0 * (22.5 * pa_yz * fx * fx * fz * pb_x

                     - 1.5 * pa_yz * fx * pb_x * fz * fgb

                     - 1.5 * pa_yz * fz * fga * pb_x * fx

                     - 3.0 * fx * fz * fga * pa_z * pb_xy

                     - pa_yyyz * pb_x * fz * fgb

                     + 15.0 * fx * fx * fz * pa_z * pb_xy

                     + 6.0 * pa_yyyz * fz * pb_x * fx

                     + 36.0 * pa_yyz * fx * fz * pb_xy

                     - 3.0 * pa_yz * fz * fga * pb_xyy

                     + 18.0 * pa_yz * fx * fz * pb_xyy

                     + 14.0 * pa_yyyz * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xyz_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pa_yyyz,
                                    double pa_yyz,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyz,
                                    double pb_xz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_x

                     + 0.75 * pa_yy * fx * fx * pb_x

                     + 0.75 * pa_y * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_z * pb_xz

                     + 0.5 * pa_yyy * fx * pb_xy

                     + 1.5 * pa_yyz * fx * pb_xz

                     + 1.5 * pa_yz * fx * pb_xyz

                     + pa_yyyz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xyz_r_0(double fga,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pa_yyyz,
                                    double pa_yyz,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyz,
                                    double pb_xz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fga * pb_x

                     + 3.0 * fx * fx * fx * fz * pb_x

                     + 7.5 * pa_yy * fx * fx * fz * pb_x

                     - 1.5 * pa_y * fz * fga * fx * pb_xy

                     - 1.5 * fx * fz * fga * pa_z * pb_xz

                     + 7.5 * pa_y * fx * fx * fz * pb_xy

                     + 7.5 * fx * fx * fz * pa_z * pb_xz

                     + 6.0 * pa_yyy * fz * fx * pb_xy

                     + 18.0 * pa_yyz * fx * fz * pb_xz

                     - 3.0 * pa_yz * fz * fga * pb_xyz

                     + 18.0 * pa_yz * fx * fz * pb_xyz

                     + 14.0 * pa_yyyz * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xzz_s_0(double fx,
                                    double pa_y,
                                    double pa_yyy,
                                    double pa_yyyz,
                                    double pa_yz,
                                    double pb_x,
                                    double pb_xz,
                                    double pb_xzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.75 * pa_yz * fx * fx * pb_x

                     + 1.5 * pa_y * fx * fx * pb_xz

                     + 0.5 * pa_yyyz * pb_x * fx

                     + pa_yyy * fx * pb_xz

                     + 1.5 * pa_yz * fx * pb_xzz

                     + pa_yyyz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yyy,
                                    double pa_yyyz,
                                    double pa_yz,
                                    double pb_x,
                                    double pb_xz,
                                    double pb_xzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_yz * fx * pb_x * fz * fgb

                     - 1.5 * pa_yz * fz * fga * pb_x * fx

                     - 3.0 * pa_y * fz * fga * fx * pb_xz

                     - pa_yyyz * pb_x * fz * fgb

                     + 7.5 * pa_yz * fx * fx * fz * pb_x

                     + 15.0 * pa_y * fx * fx * fz * pb_xz

                     + 6.0 * pa_yyyz * fz * pb_x * fx

                     + 12.0 * pa_yyy * fz * fx * pb_xz

                     - 3.0 * pa_yz * fz * fga * pb_xzz

                     + 18.0 * pa_yz * fx * fz * pb_xzz

                     + 14.0 * pa_yyyz * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_yyy_s_0(double fx,
                                    double pa_yyyz,
                                    double pa_yyz,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pa_z

                     + 2.25 * pa_yyz * fx * fx

                     + 6.75 * pa_yz * fx * fx * pb_y

                     + 2.25 * fx * fx * pa_z * pb_yy

                     + 1.5 * pa_yyyz * pb_y * fx

                     + 4.5 * pa_yyz * fx * pb_yy

                     + 1.5 * pa_yz * fx * pb_yyy

                     + pa_yyyz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_yyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yyyz,
                                    double pa_yyz,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double r_0_0)
    {
        return r_0_0 * (15.0 * fx * fx * fx * fz * pa_z

                     - 2.25 * fx * fx * pa_z * fz * fgb

                     - 2.25 * fx * fx * fz * fga * pa_z

                     - 4.5 * pa_yyz * fx * fz * fgb

                     + 22.5 * pa_yyz * fx * fx * fz

                     + 67.5 * pa_yz * fx * fx * fz * pb_y

                     - 4.5 * pa_yz * fx * pb_y * fz * fgb

                     - 4.5 * pa_yz * fz * fga * pb_y * fx

                     - 4.5 * fx * fz * fga * pa_z * pb_yy

                     - 3.0 * pa_yyyz * pb_y * fz * fgb

                     + 22.5 * fx * fx * fz * pa_z * pb_yy

                     + 18.0 * pa_yyyz * fz * pb_y * fx

                     + 54.0 * pa_yyz * fx * fz * pb_yy

                     - 3.0 * pa_yz * fz * fga * pb_yyy

                     + 18.0 * pa_yz * fx * fz * pb_yyy

                     + 14.0 * pa_yyyz * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_yyz_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pa_yyyz,
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
        return s_0_0 * (1.125 * pa_y * fx * fx * fx

                     + 0.75 * fx * fx * fx * pb_y

                     + 0.25 * pa_yyy * fx * fx

                     + 1.5 * pa_yy * fx * fx * pb_y

                     + 2.25 * pa_yz * fx * fx * pb_z

                     + 0.75 * pa_y * fx * fx * pb_yy

                     + 1.5 * fx * fx * pa_z * pb_yz

                     + 0.5 * pa_yyyz * fx * pb_z

                     + 0.5 * pa_yyy * fx * pb_yy

                     + 3.0 * pa_yyz * fx * pb_yz

                     + 1.5 * pa_yz * fx * pb_yyz

                     + pa_yyyz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_yyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pa_yyyz,
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
        return r_0_0 * (9.0 * pa_y * fx * fx * fx * fz

                     - 0.75 * pa_y * fx * fx * fz * fgb

                     - 0.75 * pa_y * fz * fga * fx * fx

                     - 1.5 * fx * fx * fz * fga * pb_y

                     - 0.5 * pa_yyy * fx * fz * fgb

                     + 6.0 * fx * fx * fx * fz * pb_y

                     + 2.5 * pa_yyy * fz * fx * fx

                     + 15.0 * pa_yy * fx * fx * fz * pb_y

                     + 22.5 * pa_yz * fx * fx * fz * pb_z

                     - 1.5 * pa_yz * fx * fz * fgb * pb_z

                     - 1.5 * pa_yz * fz * fga * fx * pb_z

                     - 1.5 * pa_y * fz * fga * fx * pb_yy

                     - 3.0 * fx * fz * fga * pa_z * pb_yz

                     - pa_yyyz * fz * fgb * pb_z

                     + 7.5 * pa_y * fx * fx * fz * pb_yy

                     + 15.0 * fx * fx * fz * pa_z * pb_yz

                     + 6.0 * pa_yyyz * fz * fx * pb_z

                     + 6.0 * pa_yyy * fz * fx * pb_yy

                     + 36.0 * pa_yyz * fx * fz * pb_yz

                     - 3.0 * pa_yz * fz * fga * pb_yyz

                     + 18.0 * pa_yz * fx * fz * pb_yyz

                     + 14.0 * pa_yyyz * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_yzz_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pa_yyyz,
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
        return s_0_0 * (0.375 * fx * fx * fx * pa_z

                     + 0.75 * fx * fx * fx * pb_z

                     + 0.75 * pa_yyz * fx * fx

                     + 1.5 * pa_yy * fx * fx * pb_z

                     + 0.75 * pa_yz * fx * fx * pb_y

                     + 1.5 * pa_y * fx * fx * pb_yz

                     + 0.75 * fx * fx * pa_z * pb_zz

                     + 0.5 * pa_yyyz * pb_y * fx

                     + pa_yyy * fx * pb_yz

                     + 1.5 * pa_yyz * fx * pb_zz

                     + 1.5 * pa_yz * fx * pb_yzz

                     + pa_yyyz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_yzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pa_yyyz,
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
        return r_0_0 * (- 0.75 * fx * fx * pa_z * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pa_z

                     - 1.5 * fx * fx * fz * fga * pb_z

                     - 1.5 * pa_yyz * fx * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_z

                     + 6.0 * fx * fx * fx * fz * pb_z

                     + 7.5 * pa_yyz * fx * fx * fz

                     + 15.0 * pa_yy * fx * fx * fz * pb_z

                     - 1.5 * pa_yz * fx * pb_y * fz * fgb

                     - 1.5 * pa_yz * fz * fga * pb_y * fx

                     - 3.0 * pa_y * fz * fga * fx * pb_yz

                     - 1.5 * fx * fz * fga * pa_z * pb_zz

                     - pa_yyyz * pb_y * fz * fgb

                     + 7.5 * pa_yz * fx * fx * fz * pb_y

                     + 15.0 * pa_y * fx * fx * fz * pb_yz

                     + 7.5 * fx * fx * fz * pa_z * pb_zz

                     + 6.0 * pa_yyyz * fz * pb_y * fx

                     + 12.0 * pa_yyy * fz * fx * pb_yz

                     + 18.0 * pa_yyz * fx * fz * pb_zz

                     - 3.0 * pa_yz * fz * fga * pb_yzz

                     + 18.0 * pa_yz * fx * fz * pb_yzz

                     + 14.0 * pa_yyyz * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_zzz_s_0(double fx,
                                    double pa_y,
                                    double pa_yyy,
                                    double pa_yyyz,
                                    double pa_yz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * pa_y * fx * fx * fx

                     + 0.75 * pa_yyy * fx * fx

                     + 2.25 * pa_yz * fx * fx * pb_z

                     + 2.25 * pa_y * fx * fx * pb_zz

                     + 1.5 * pa_yyyz * pb_z * fx

                     + 1.5 * pa_yyy * fx * pb_zz

                     + 1.5 * pa_yz * fx * pb_zzz

                     + pa_yyyz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_zzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yyy,
                                    double pa_yyyz,
                                    double pa_yz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 2.25 * pa_y * fx * fx * fz * fgb

                     - 2.25 * pa_y * fz * fga * fx * fx

                     - 1.5 * pa_yyy * fx * fz * fgb

                     + 9.0 * pa_y * fx * fx * fx * fz

                     + 7.5 * pa_yyy * fz * fx * fx

                     - 4.5 * pa_yz * fx * pb_z * fz * fgb

                     - 4.5 * pa_yz * fz * fga * pb_z * fx

                     - 4.5 * pa_y * fz * fga * fx * pb_zz

                     - 3.0 * pa_yyyz * pb_z * fz * fgb

                     + 22.5 * pa_yz * fx * fx * fz * pb_z

                     + 22.5 * pa_y * fx * fx * fz * pb_zz

                     + 18.0 * pa_yyyz * fz * pb_z * fx

                     + 18.0 * pa_yyy * fz * fx * pb_zz

                     - 3.0 * pa_yz * fz * fga * pb_zzz

                     + 18.0 * pa_yz * fx * fz * pb_zzz

                     + 14.0 * pa_yyyz * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xxx_s_0(double fx,
                                    double pa_yy,
                                    double pa_yyzz,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xxx,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_x

                     + 0.75 * pa_yy * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_zz * pb_x

                     + 1.5 * pa_yyzz * pb_x * fx

                     + 0.25 * fx * fx * pb_xxx

                     + 0.5 * pa_yy * fx * pb_xxx

                     + 0.5 * fx * pa_zz * pb_xxx

                     + pa_yyzz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yy,
                                    double pa_yyzz,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pb_x * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_x

                     - 1.5 * pa_yy * fx * pb_x * fz * fgb

                     - 1.5 * pa_yy * fz * fga * pb_x * fx

                     - 1.5 * fx * pa_zz * pb_x * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_x

                     - 1.5 * fz * fga * pa_zz * pb_x * fx

                     - 3.0 * pa_yyzz * pb_x * fz * fgb

                     + 7.5 * pa_yy * fz * fx * fx * pb_x

                     + 7.5 * fx * fx * pa_zz * fz * pb_x

                     + 18.0 * pa_yyzz * fz * pb_x * fx

                     - fx * fz * fga * pb_xxx

                     - pa_yy * fz * fga * pb_xxx

                     + 2.5 * fx * fx * fz * pb_xxx

                     - fz * fga * pa_zz * pb_xxx

                     + 6.0 * pa_yy * fz * fx * pb_xxx

                     + 6.0 * fx * pa_zz * fz * pb_xxx

                     + 14.0 * pa_yyzz * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xxy_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyzz,
                                    double pa_yzz,
                                    double pa_zz,
                                    double pb_xx,
                                    double pb_xxy,
                                    double pb_y,
                                    double s_0_0)
    {
        return s_0_0 * (0.25 * pa_y * fx * fx * fx

                     + 0.5 * pa_yzz * fx * fx

                     + 0.125 * fx * fx * fx * pb_y

                     + 0.25 * pa_yy * fx * fx * pb_y

                     + 0.5 * pa_y * fx * fx * pb_xx

                     + 0.25 * fx * fx * pa_zz * pb_y

                     + 0.5 * pa_yyzz * fx * pb_y

                     + pa_yzz * fx * pb_xx

                     + 0.25 * fx * fx * pb_xxy

                     + 0.5 * pa_yy * fx * pb_xxy

                     + 0.5 * fx * pa_zz * pb_xxy

                     + pa_yyzz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyzz,
                                    double pa_yzz,
                                    double pa_zz,
                                    double pb_xx,
                                    double pb_xxy,
                                    double pb_y,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_y * fx * fx * fz * fgb

                     - 0.5 * pa_y * fx * fx * fz * fga

                     - pa_yzz * fx * fz * fgb

                     + 2.0 * pa_y * fx * fx * fx * fz

                     + 5.0 * pa_yzz * fx * fx * fz

                     - 0.25 * fx * fx * fz * fgb * pb_y

                     - 0.5 * fx * fx * fz * fga * pb_y

                     - 0.5 * pa_yy * fx * fz * fgb * pb_y

                     - 0.5 * pa_yy * fz * fga * fx * pb_y

                     - pa_y * fx * fz * fga * pb_xx

                     - 0.5 * fx * pa_zz * fz * fgb * pb_y

                     + fx * fx * fx * fz * pb_y

                     - 0.5 * fz * fga * pa_zz * fx * pb_y

                     - pa_yyzz * fz * fgb * pb_y

                     + 2.5 * pa_yy * fz * fx * fx * pb_y

                     + 5.0 * pa_y * fx * fx * fz * pb_xx

                     + 2.5 * fx * fx * pa_zz * fz * pb_y

                     + 6.0 * pa_yyzz * fz * fx * pb_y

                     + 12.0 * pa_yzz * fx * fz * pb_xx

                     - fx * fz * fga * pb_xxy

                     - pa_yy * fz * fga * pb_xxy

                     + 2.5 * fx * fx * fz * pb_xxy

                     - fz * fga * pa_zz * pb_xxy

                     + 6.0 * pa_yy * fz * fx * pb_xxy

                     + 6.0 * fx * pa_zz * fz * pb_xxy

                     + 14.0 * pa_yyzz * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xxz_s_0(double fx,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_yyzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pb_xx,
                                    double pb_xxz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * fx * pa_z

                     + 0.5 * pa_yyz * fx * fx

                     + 0.125 * fx * fx * fx * pb_z

                     + 0.25 * pa_yy * fx * fx * pb_z

                     + 0.25 * fx * fx * pa_zz * pb_z

                     + 0.5 * fx * fx * pa_z * pb_xx

                     + 0.5 * pa_yyzz * fx * pb_z

                     + pa_yyz * fx * pb_xx

                     + 0.25 * fx * fx * pb_xxz

                     + 0.5 * pa_yy * fx * pb_xxz

                     + 0.5 * fx * pa_zz * pb_xxz

                     + pa_yyzz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_yyzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pb_xx,
                                    double pb_xxz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * pa_z * fz * fgb

                     - 0.5 * fz * fga * pa_z * fx * fx

                     - pa_yyz * fx * fz * fgb

                     + 2.0 * fx * fx * fx * pa_z * fz

                     + 5.0 * pa_yyz * fz * fx * fx

                     - 0.25 * fx * fx * fz * fgb * pb_z

                     - 0.5 * fx * fx * fz * fga * pb_z

                     - 0.5 * pa_yy * fx * fz * fgb * pb_z

                     - 0.5 * pa_yy * fz * fga * fx * pb_z

                     - 0.5 * fx * pa_zz * fz * fgb * pb_z

                     + fx * fx * fx * fz * pb_z

                     - 0.5 * fz * fga * pa_zz * fx * pb_z

                     - fz * fga * pa_z * fx * pb_xx

                     - pa_yyzz * fz * fgb * pb_z

                     + 2.5 * pa_yy * fz * fx * fx * pb_z

                     + 2.5 * fx * fx * pa_zz * fz * pb_z

                     + 5.0 * fx * fx * pa_z * fz * pb_xx

                     + 6.0 * pa_yyzz * fz * fx * pb_z

                     + 12.0 * pa_yyz * fz * fx * pb_xx

                     - fx * fz * fga * pb_xxz

                     - pa_yy * fz * fga * pb_xxz

                     + 2.5 * fx * fx * fz * pb_xxz

                     - fz * fga * pa_zz * pb_xxz

                     + 6.0 * pa_yy * fz * fx * pb_xxz

                     + 6.0 * fx * pa_zz * fz * pb_xxz

                     + 14.0 * pa_yyzz * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xyy_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyzz,
                                    double pa_yzz,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyy,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_zz * pb_x

                     + 0.25 * pa_yy * fx * fx * pb_x

                     + pa_y * fx * fx * pb_xy

                     + 0.5 * pa_yyzz * pb_x * fx

                     + 2.0 * pa_yzz * fx * pb_xy

                     + 0.25 * fx * fx * pb_xyy

                     + 0.5 * pa_yy * fx * pb_xyy

                     + 0.5 * fx * pa_zz * pb_xyy

                     + pa_yyzz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyzz,
                                    double pa_yzz,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyy,
                                    double r_0_0)
    {
        return r_0_0 * (- fx * fx * fz * fga * pb_x

                     + 3.0 * fx * fx * fx * fz * pb_x

                     + 7.5 * fx * fx * pa_zz * fz * pb_x

                     - 0.25 * fx * fx * pb_x * fz * fgb

                     - 0.5 * pa_yy * fx * pb_x * fz * fgb

                     - 0.5 * pa_yy * fz * fga * pb_x * fx

                     - 2.0 * pa_y * fx * fz * fga * pb_xy

                     - 0.5 * fx * pa_zz * pb_x * fz * fgb

                     - 0.5 * fz * fga * pa_zz * pb_x * fx

                     - pa_yyzz * pb_x * fz * fgb

                     + 2.5 * pa_yy * fz * fx * fx * pb_x

                     + 10.0 * pa_y * fx * fx * fz * pb_xy

                     + 6.0 * pa_yyzz * fz * pb_x * fx

                     + 24.0 * pa_yzz * fx * fz * pb_xy

                     - fx * fz * fga * pb_xyy

                     - pa_yy * fz * fga * pb_xyy

                     + 2.5 * fx * fx * fz * pb_xyy

                     - fz * fga * pa_zz * pb_xyy

                     + 6.0 * pa_yy * fz * fx * pb_xyy

                     + 6.0 * fx * pa_zz * fz * pb_xyy

                     + 14.0 * pa_yyzz * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xyz_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_yyzz,
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
        return s_0_0 * (pa_yz * fx * fx * pb_x

                     + 0.5 * pa_y * fx * fx * pb_xz

                     + 0.5 * fx * fx * pa_z * pb_xy

                     + pa_yyz * fx * pb_xy

                     + pa_yzz * fx * pb_xz

                     + 0.25 * fx * fx * pb_xyz

                     + 0.5 * pa_yy * fx * pb_xyz

                     + 0.5 * fx * pa_zz * pb_xyz

                     + pa_yyzz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xyz_r_0(double fga,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_yyzz,
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
        return r_0_0 * (10.0 * pa_yz * fx * fx * fz * pb_x

                     - pa_y * fx * fz * fga * pb_xz

                     - fz * fga * pa_z * fx * pb_xy

                     + 5.0 * pa_y * fx * fx * fz * pb_xz

                     + 5.0 * fx * fx * pa_z * fz * pb_xy

                     + 12.0 * pa_yyz * fz * fx * pb_xy

                     + 12.0 * pa_yzz * fx * fz * pb_xz

                     - fx * fz * fga * pb_xyz

                     - pa_yy * fz * fga * pb_xyz

                     + 2.5 * fx * fx * fz * pb_xyz

                     - fz * fga * pa_zz * pb_xyz

                     + 6.0 * pa_yy * fz * fx * pb_xyz

                     + 6.0 * fx * pa_zz * fz * pb_xyz

                     + 14.0 * pa_yyzz * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xzz_s_0(double fx,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_yyzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xz,
                                    double pb_xzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_x

                     + 0.75 * pa_yy * fx * fx * pb_x

                     + 0.25 * fx * fx * pa_zz * pb_x

                     + fx * fx * pa_z * pb_xz

                     + 0.5 * pa_yyzz * pb_x * fx

                     + 2.0 * pa_yyz * fx * pb_xz

                     + 0.25 * fx * fx * pb_xzz

                     + 0.5 * pa_yy * fx * pb_xzz

                     + 0.5 * fx * pa_zz * pb_xzz

                     + pa_yyzz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_yyzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xz,
                                    double pb_xzz,
                                    double r_0_0)
    {
        return r_0_0 * (- fz * fga * fx * fx * pb_x

                     + 3.0 * fx * fx * fx * fz * pb_x

                     + 7.5 * pa_yy * fx * fx * fz * pb_x

                     - 0.25 * fx * fx * pb_x * fz * fgb

                     - 0.5 * pa_yy * fx * pb_x * fz * fgb

                     - 0.5 * pa_yy * fz * fga * pb_x * fx

                     - 0.5 * fx * pa_zz * pb_x * fz * fgb

                     - 0.5 * fz * fga * pa_zz * pb_x * fx

                     - 2.0 * fz * fga * pa_z * fx * pb_xz

                     - pa_yyzz * pb_x * fz * fgb

                     + 2.5 * fx * fx * pa_zz * fz * pb_x

                     + 10.0 * fx * fx * pa_z * fz * pb_xz

                     + 6.0 * pa_yyzz * fz * pb_x * fx

                     + 24.0 * pa_yyz * fz * fx * pb_xz

                     - fx * fz * fga * pb_xzz

                     - pa_yy * fz * fga * pb_xzz

                     + 2.5 * fx * fx * fz * pb_xzz

                     - fz * fga * pa_zz * pb_xzz

                     + 6.0 * pa_yy * fz * fx * pb_xzz

                     + 6.0 * fx * pa_zz * fz * pb_xzz

                     + 14.0 * pa_yyzz * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_yyy_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyzz,
                                    double pa_yzz,
                                    double pa_zz,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx * fx

                     + 1.125 * fx * fx * fx * pb_y

                     + 1.5 * pa_yzz * fx * fx

                     + 2.25 * fx * fx * pa_zz * pb_y

                     + 0.75 * pa_yy * fx * fx * pb_y

                     + 1.5 * pa_y * fx * fx * pb_yy

                     + 1.5 * pa_yyzz * pb_y * fx

                     + 3.0 * pa_yzz * fx * pb_yy

                     + 0.25 * fx * fx * pb_yyy

                     + 0.5 * pa_yy * fx * pb_yyy

                     + 0.5 * fx * pa_zz * pb_yyy

                     + pa_yyzz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_yyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyzz,
                                    double pa_yzz,
                                    double pa_zz,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_y * fx * fx * fz * fgb

                     - 1.5 * pa_y * fx * fx * fz * fga

                     - 3.0 * fx * fx * fz * fga * pb_y

                     - 3.0 * pa_yzz * fx * fz * fgb

                     + 6.0 * pa_y * fx * fx * fx * fz

                     + 9.0 * fx * fx * fx * fz * pb_y

                     + 15.0 * pa_yzz * fx * fx * fz

                     + 22.5 * fx * fx * pa_zz * fz * pb_y

                     - 0.75 * fx * fx * pb_y * fz * fgb

                     - 1.5 * pa_yy * fx * pb_y * fz * fgb

                     - 1.5 * pa_yy * fz * fga * pb_y * fx

                     - 3.0 * pa_y * fx * fz * fga * pb_yy

                     - 1.5 * fx * pa_zz * pb_y * fz * fgb

                     - 1.5 * fz * fga * pa_zz * pb_y * fx

                     - 3.0 * pa_yyzz * pb_y * fz * fgb

                     + 7.5 * pa_yy * fz * fx * fx * pb_y

                     + 15.0 * pa_y * fx * fx * fz * pb_yy

                     + 18.0 * pa_yyzz * fz * pb_y * fx

                     + 36.0 * pa_yzz * fx * fz * pb_yy

                     - fx * fz * fga * pb_yyy

                     - pa_yy * fz * fga * pb_yyy

                     + 2.5 * fx * fx * fz * pb_yyy

                     - fz * fga * pa_zz * pb_yyy

                     + 6.0 * pa_yy * fz * fx * pb_yyy

                     + 6.0 * fx * pa_zz * fz * pb_yyy

                     + 14.0 * pa_yyzz * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_yyz_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_yyzz,
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
        return s_0_0 * (0.75 * fx * fx * fx * pa_z

                     + 0.375 * fx * fx * fx * pb_z

                     + 0.5 * pa_yyz * fx * fx

                     + 2.0 * pa_yz * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_zz * pb_z

                     + 0.25 * pa_yy * fx * fx * pb_z

                     + pa_y * fx * fx * pb_yz

                     + 0.5 * fx * fx * pa_z * pb_yy

                     + 0.5 * pa_yyzz * fx * pb_z

                     + pa_yyz * fx * pb_yy

                     + 2.0 * pa_yzz * fx * pb_yz

                     + 0.25 * fx * fx * pb_yyz

                     + 0.5 * pa_yy * fx * pb_yyz

                     + 0.5 * fx * pa_zz * pb_yyz

                     + pa_yyzz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_yyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_yyzz,
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
        return r_0_0 * (6.0 * fx * fx * fx * pa_z * fz

                     - 0.5 * fx * fx * pa_z * fz * fgb

                     - fx * fx * fz * fga * pb_z

                     - 0.5 * fz * fga * pa_z * fx * fx

                     - pa_yyz * fx * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_z

                     + 5.0 * pa_yyz * fz * fx * fx

                     + 20.0 * pa_yz * fx * fx * fz * pb_y

                     + 7.5 * fx * fx * pa_zz * fz * pb_z

                     - 0.25 * fx * fx * fz * fgb * pb_z

                     - 0.5 * pa_yy * fx * fz * fgb * pb_z

                     - 0.5 * pa_yy * fz * fga * fx * pb_z

                     - 2.0 * pa_y * fx * fz * fga * pb_yz

                     - 0.5 * fx * pa_zz * fz * fgb * pb_z

                     - 0.5 * fz * fga * pa_zz * fx * pb_z

                     - fz * fga * pa_z * fx * pb_yy

                     - pa_yyzz * fz * fgb * pb_z

                     + 2.5 * pa_yy * fz * fx * fx * pb_z

                     + 10.0 * pa_y * fx * fx * fz * pb_yz

                     + 5.0 * fx * fx * pa_z * fz * pb_yy

                     + 6.0 * pa_yyzz * fz * fx * pb_z

                     + 12.0 * pa_yyz * fz * fx * pb_yy

                     + 24.0 * pa_yzz * fx * fz * pb_yz

                     - fx * fz * fga * pb_yyz

                     - pa_yy * fz * fga * pb_yyz

                     + 2.5 * fx * fx * fz * pb_yyz

                     - fz * fga * pa_zz * pb_yyz

                     + 6.0 * pa_yy * fz * fx * pb_yyz

                     + 6.0 * fx * pa_zz * fz * pb_yyz

                     + 14.0 * pa_yyzz * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_yzz_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_yyzz,
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
        return s_0_0 * (0.75 * pa_y * fx * fx * fx

                     + 0.375 * fx * fx * fx * pb_y

                     + 0.75 * pa_yy * fx * fx * pb_y

                     + 0.5 * pa_yzz * fx * fx

                     + 2.0 * pa_yz * fx * fx * pb_z

                     + 0.5 * pa_y * fx * fx * pb_zz

                     + 0.25 * fx * fx * pa_zz * pb_y

                     + fx * fx * pa_z * pb_yz

                     + 0.5 * pa_yyzz * pb_y * fx

                     + 2.0 * pa_yyz * fx * pb_yz

                     + pa_yzz * fx * pb_zz

                     + 0.25 * fx * fx * pb_yzz

                     + 0.5 * pa_yy * fx * pb_yzz

                     + 0.5 * fx * pa_zz * pb_yzz

                     + pa_yyzz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_yzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_yyzz,
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
        return r_0_0 * (6.0 * pa_y * fx * fx * fx * fz

                     - 0.5 * pa_y * fx * fx * fz * fgb

                     - 0.5 * pa_y * fx * fx * fz * fga

                     - fz * fga * fx * fx * pb_y

                     - pa_yzz * fx * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_y

                     + 7.5 * pa_yy * fx * fx * fz * pb_y

                     + 5.0 * pa_yzz * fx * fx * fz

                     + 20.0 * pa_yz * fx * fx * fz * pb_z

                     - 0.25 * fx * fx * pb_y * fz * fgb

                     - 0.5 * pa_yy * fx * pb_y * fz * fgb

                     - 0.5 * pa_yy * fz * fga * pb_y * fx

                     - pa_y * fx * fz * fga * pb_zz

                     - 0.5 * fx * pa_zz * pb_y * fz * fgb

                     - 0.5 * fz * fga * pa_zz * pb_y * fx

                     - 2.0 * fz * fga * pa_z * fx * pb_yz

                     - pa_yyzz * pb_y * fz * fgb

                     + 5.0 * pa_y * fx * fx * fz * pb_zz

                     + 2.5 * fx * fx * pa_zz * fz * pb_y

                     + 10.0 * fx * fx * pa_z * fz * pb_yz

                     + 6.0 * pa_yyzz * fz * pb_y * fx

                     + 24.0 * pa_yyz * fz * fx * pb_yz

                     + 12.0 * pa_yzz * fx * fz * pb_zz

                     - fx * fz * fga * pb_yzz

                     - pa_yy * fz * fga * pb_yzz

                     + 2.5 * fx * fx * fz * pb_yzz

                     - fz * fga * pa_zz * pb_yzz

                     + 6.0 * pa_yy * fz * fx * pb_yzz

                     + 6.0 * fx * pa_zz * fz * pb_yzz

                     + 14.0 * pa_yyzz * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_zzz_s_0(double fx,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_yyzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * fx * pa_z

                     + 1.125 * fx * fx * fx * pb_z

                     + 1.5 * pa_yyz * fx * fx

                     + 2.25 * pa_yy * fx * fx * pb_z

                     + 0.75 * fx * fx * pa_zz * pb_z

                     + 1.5 * fx * fx * pa_z * pb_zz

                     + 1.5 * pa_yyzz * pb_z * fx

                     + 3.0 * pa_yyz * fx * pb_zz

                     + 0.25 * fx * fx * pb_zzz

                     + 0.5 * pa_yy * fx * pb_zzz

                     + 0.5 * fx * pa_zz * pb_zzz

                     + pa_yyzz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_zzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_yyzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_z * fz * fgb

                     - 1.5 * fz * fga * pa_z * fx * fx

                     - 3.0 * fz * fga * fx * fx * pb_z

                     - 3.0 * pa_yyz * fx * fz * fgb

                     + 6.0 * fx * fx * fx * pa_z * fz

                     + 9.0 * fx * fx * fx * fz * pb_z

                     + 15.0 * pa_yyz * fz * fx * fx

                     + 22.5 * pa_yy * fx * fx * fz * pb_z

                     - 0.75 * fx * fx * pb_z * fz * fgb

                     - 1.5 * pa_yy * fx * pb_z * fz * fgb

                     - 1.5 * pa_yy * fz * fga * pb_z * fx

                     - 1.5 * fx * pa_zz * pb_z * fz * fgb

                     - 1.5 * fz * fga * pa_zz * pb_z * fx

                     - 3.0 * fz * fga * pa_z * fx * pb_zz

                     - 3.0 * pa_yyzz * pb_z * fz * fgb

                     + 7.5 * fx * fx * pa_zz * fz * pb_z

                     + 15.0 * fx * fx * pa_z * fz * pb_zz

                     + 18.0 * pa_yyzz * fz * pb_z * fx

                     + 36.0 * pa_yyz * fz * fx * pb_zz

                     - fx * fz * fga * pb_zzz

                     - pa_yy * fz * fga * pb_zzz

                     + 2.5 * fx * fx * fz * pb_zzz

                     - fz * fga * pa_zz * pb_zzz

                     + 6.0 * pa_yy * fz * fx * pb_zzz

                     + 6.0 * fx * pa_zz * fz * pb_zzz

                     + 14.0 * pa_yyzz * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xxx_s_0(double fx,
                                    double pa_yz,
                                    double pa_yzzz,
                                    double pb_x,
                                    double pb_xxx,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_yz * fx * fx * pb_x

                     + 1.5 * pa_yzzz * pb_x * fx

                     + 1.5 * pa_yz * fx * pb_xxx

                     + pa_yzzz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yz,
                                    double pa_yzzz,
                                    double pb_x,
                                    double pb_xxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_yz * fx * pb_x * fz * fgb

                     - 4.5 * pa_yz * fz * fga * pb_x * fx

                     - 3.0 * pa_yzzz * pb_x * fz * fgb

                     + 22.5 * pa_yz * fz * fx * fx * pb_x

                     + 18.0 * pa_yzzz * fz * pb_x * fx

                     - 3.0 * pa_yz * fz * fga * pb_xxx

                     + 18.0 * pa_yz * fz * fx * pb_xxx

                     + 14.0 * pa_yzzz * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xxy_s_0(double fx,
                                    double pa_yz,
                                    double pa_yzzz,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_xx,
                                    double pb_xxy,
                                    double pb_y,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_z

                     + 0.25 * fx * fx * pa_zzz

                     + 0.75 * pa_yz * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_z * pb_xx

                     + 0.5 * pa_yzzz * fx * pb_y

                     + 0.5 * fx * pa_zzz * pb_xx

                     + 1.5 * pa_yz * fx * pb_xxy

                     + pa_yzzz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yz,
                                    double pa_yzzz,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_xx,
                                    double pb_xxy,
                                    double pb_y,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pa_z * fz * fgb

                     - 0.75 * fx * fx * pa_z * fz * fga

                     - 0.5 * fx * pa_zzz * fz * fgb

                     + 3.0 * fx * fx * fx * pa_z * fz

                     + 2.5 * fx * fx * pa_zzz * fz

                     - 1.5 * pa_yz * fx * fz * fgb * pb_y

                     - 1.5 * pa_yz * fz * fga * fx * pb_y

                     - 1.5 * fx * pa_z * fz * fga * pb_xx

                     - pa_yzzz * fz * fgb * pb_y

                     + 7.5 * pa_yz * fz * fx * fx * pb_y

                     + 7.5 * fx * fx * pa_z * fz * pb_xx

                     + 6.0 * pa_yzzz * fz * fx * pb_y

                     + 6.0 * fx * pa_zzz * fz * pb_xx

                     - 3.0 * pa_yz * fz * fga * pb_xxy

                     + 18.0 * pa_yz * fz * fx * pb_xxy

                     + 14.0 * pa_yzzz * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xxz_s_0(double fx,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pa_yzzz,
                                    double pb_xx,
                                    double pb_xxz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_y * fx * fx * fx

                     + 0.75 * pa_yzz * fx * fx

                     + 0.75 * pa_yz * fx * fx * pb_z

                     + 0.75 * pa_y * fx * fx * pb_xx

                     + 0.5 * pa_yzzz * fx * pb_z

                     + 1.5 * pa_yzz * fx * pb_xx

                     + 1.5 * pa_yz * fx * pb_xxz

                     + pa_yzzz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pa_yzzz,
                                    double pb_xx,
                                    double pb_xxz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * pa_y * fx * fx * fz * fgb

                     - 0.75 * pa_y * fx * fx * fz * fga

                     - 1.5 * pa_yzz * fx * fz * fgb

                     + 3.0 * pa_y * fx * fx * fx * fz

                     + 7.5 * pa_yzz * fz * fx * fx

                     - 1.5 * pa_yz * fx * fz * fgb * pb_z

                     - 1.5 * pa_yz * fz * fga * fx * pb_z

                     - 1.5 * pa_y * fx * fz * fga * pb_xx

                     - pa_yzzz * fz * fgb * pb_z

                     + 7.5 * pa_yz * fz * fx * fx * pb_z

                     + 7.5 * pa_y * fx * fx * fz * pb_xx

                     + 6.0 * pa_yzzz * fz * fx * pb_z

                     + 18.0 * pa_yzz * fz * fx * pb_xx

                     - 3.0 * pa_yz * fz * fga * pb_xxz

                     + 18.0 * pa_yz * fz * fx * pb_xxz

                     + 14.0 * pa_yzzz * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xyy_s_0(double fx,
                                    double pa_yz,
                                    double pa_yzzz,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyy,
                                    double s_0_0)
    {
        return s_0_0 * (0.75 * pa_yz * fx * fx * pb_x

                     + 1.5 * fx * fx * pa_z * pb_xy

                     + 0.5 * pa_yzzz * pb_x * fx

                     + fx * pa_zzz * pb_xy

                     + 1.5 * pa_yz * fx * pb_xyy

                     + pa_yzzz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yz,
                                    double pa_yzzz,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_yz * fx * pb_x * fz * fgb

                     - 1.5 * pa_yz * fz * fga * pb_x * fx

                     - 3.0 * fx * pa_z * fz * fga * pb_xy

                     - pa_yzzz * pb_x * fz * fgb

                     + 7.5 * pa_yz * fz * fx * fx * pb_x

                     + 15.0 * fx * fx * pa_z * fz * pb_xy

                     + 6.0 * pa_yzzz * fz * pb_x * fx

                     + 12.0 * fx * pa_zzz * fz * pb_xy

                     - 3.0 * pa_yz * fz * fga * pb_xyy

                     + 18.0 * pa_yz * fz * fx * pb_xyy

                     + 14.0 * pa_yzzz * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xyz_s_0(double fx,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pa_yzzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyz,
                                    double pb_xz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_zz * pb_x

                     + 0.75 * pa_y * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_z * pb_xz

                     + 1.5 * pa_yzz * fx * pb_xy

                     + 0.5 * fx * pa_zzz * pb_xz

                     + 1.5 * pa_yz * fx * pb_xyz

                     + pa_yzzz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xyz_r_0(double fga,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pa_yzzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyz,
                                    double pb_xz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fga * pb_x

                     + 3.0 * fx * fx * fx * fz * pb_x

                     + 7.5 * fx * fx * pa_zz * fz * pb_x

                     - 1.5 * pa_y * fx * fz * fga * pb_xy

                     - 1.5 * fx * pa_z * fz * fga * pb_xz

                     + 7.5 * pa_y * fx * fx * fz * pb_xy

                     + 7.5 * fx * fx * pa_z * fz * pb_xz

                     + 18.0 * pa_yzz * fz * fx * pb_xy

                     + 6.0 * fx * pa_zzz * fz * pb_xz

                     - 3.0 * pa_yz * fz * fga * pb_xyz

                     + 18.0 * pa_yz * fz * fx * pb_xyz

                     + 14.0 * pa_yzzz * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xzz_s_0(double fx,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pa_yzzz,
                                    double pb_x,
                                    double pb_xz,
                                    double pb_xzz,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_yz * fx * fx * pb_x

                     + 1.5 * pa_y * fx * fx * pb_xz

                     + 0.5 * pa_yzzz * pb_x * fx

                     + 3.0 * pa_yzz * fx * pb_xz

                     + 1.5 * pa_yz * fx * pb_xzz

                     + pa_yzzz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pa_yzzz,
                                    double pb_x,
                                    double pb_xz,
                                    double pb_xzz,
                                    double r_0_0)
    {
        return r_0_0 * (22.5 * pa_yz * fx * fx * fz * pb_x

                     - 1.5 * pa_yz * fx * pb_x * fz * fgb

                     - 1.5 * pa_yz * fz * fga * pb_x * fx

                     - 3.0 * pa_y * fx * fz * fga * pb_xz

                     - pa_yzzz * pb_x * fz * fgb

                     + 15.0 * pa_y * fx * fx * fz * pb_xz

                     + 6.0 * pa_yzzz * fz * pb_x * fx

                     + 36.0 * pa_yzz * fz * fx * pb_xz

                     - 3.0 * pa_yz * fz * fga * pb_xzz

                     + 18.0 * pa_yz * fz * fx * pb_xzz

                     + 14.0 * pa_yzzz * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_yyy_s_0(double fx,
                                    double pa_yz,
                                    double pa_yzzz,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pa_z

                     + 0.75 * fx * fx * pa_zzz

                     + 2.25 * pa_yz * fx * fx * pb_y

                     + 2.25 * fx * fx * pa_z * pb_yy

                     + 1.5 * pa_yzzz * pb_y * fx

                     + 1.5 * fx * pa_zzz * pb_yy

                     + 1.5 * pa_yz * fx * pb_yyy

                     + pa_yzzz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_yyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yz,
                                    double pa_yzzz,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pa_z * fz * fgb

                     - 2.25 * fx * fx * pa_z * fz * fga

                     - 1.5 * fx * pa_zzz * fz * fgb

                     + 9.0 * fx * fx * fx * pa_z * fz

                     + 7.5 * fx * fx * pa_zzz * fz

                     - 4.5 * pa_yz * fx * pb_y * fz * fgb

                     - 4.5 * pa_yz * fz * fga * pb_y * fx

                     - 4.5 * fx * pa_z * fz * fga * pb_yy

                     - 3.0 * pa_yzzz * pb_y * fz * fgb

                     + 22.5 * pa_yz * fz * fx * fx * pb_y

                     + 22.5 * fx * fx * pa_z * fz * pb_yy

                     + 18.0 * pa_yzzz * fz * pb_y * fx

                     + 18.0 * fx * pa_zzz * fz * pb_yy

                     - 3.0 * pa_yz * fz * fga * pb_yyy

                     + 18.0 * pa_yz * fz * fx * pb_yyy

                     + 14.0 * pa_yzzz * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_yyz_s_0(double fx,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pa_yzzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_y * fx * fx * fx

                     + 0.75 * fx * fx * fx * pb_y

                     + 0.75 * pa_yzz * fx * fx

                     + 1.5 * fx * fx * pa_zz * pb_y

                     + 0.75 * pa_yz * fx * fx * pb_z

                     + 0.75 * pa_y * fx * fx * pb_yy

                     + 1.5 * fx * fx * pa_z * pb_yz

                     + 0.5 * pa_yzzz * fx * pb_z

                     + 1.5 * pa_yzz * fx * pb_yy

                     + fx * pa_zzz * pb_yz

                     + 1.5 * pa_yz * fx * pb_yyz

                     + pa_yzzz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_yyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pa_yzzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * pa_y * fx * fx * fz * fgb

                     - 0.75 * pa_y * fx * fx * fz * fga

                     - 1.5 * fx * fx * fz * fga * pb_y

                     - 1.5 * pa_yzz * fx * fz * fgb

                     + 3.0 * pa_y * fx * fx * fx * fz

                     + 6.0 * fx * fx * fx * fz * pb_y

                     + 7.5 * pa_yzz * fz * fx * fx

                     + 15.0 * fx * fx * pa_zz * fz * pb_y

                     - 1.5 * pa_yz * fx * fz * fgb * pb_z

                     - 1.5 * pa_yz * fz * fga * fx * pb_z

                     - 1.5 * pa_y * fx * fz * fga * pb_yy

                     - 3.0 * fx * pa_z * fz * fga * pb_yz

                     - pa_yzzz * fz * fgb * pb_z

                     + 7.5 * pa_yz * fz * fx * fx * pb_z

                     + 7.5 * pa_y * fx * fx * fz * pb_yy

                     + 15.0 * fx * fx * pa_z * fz * pb_yz

                     + 6.0 * pa_yzzz * fz * fx * pb_z

                     + 18.0 * pa_yzz * fz * fx * pb_yy

                     + 12.0 * fx * pa_zzz * fz * pb_yz

                     - 3.0 * pa_yz * fz * fga * pb_yyz

                     + 18.0 * pa_yz * fz * fx * pb_yyz

                     + 14.0 * pa_yzzz * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_yzz_s_0(double fx,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pa_yzzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double pb_z,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pa_z

                     + 0.75 * fx * fx * fx * pb_z

                     + 2.25 * pa_yz * fx * fx * pb_y

                     + 0.25 * fx * fx * pa_zzz

                     + 1.5 * fx * fx * pa_zz * pb_z

                     + 1.5 * pa_y * fx * fx * pb_yz

                     + 0.75 * fx * fx * pa_z * pb_zz

                     + 0.5 * pa_yzzz * pb_y * fx

                     + 3.0 * pa_yzz * fx * pb_yz

                     + 0.5 * fx * pa_zzz * pb_zz

                     + 1.5 * pa_yz * fx * pb_yzz

                     + pa_yzzz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_yzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pa_yzzz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double pb_z,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (9.0 * fx * fx * fx * pa_z * fz

                     - 0.75 * fx * fx * pa_z * fz * fgb

                     - 0.75 * fx * fx * pa_z * fz * fga

                     - 1.5 * fx * fx * fz * fga * pb_z

                     - 0.5 * fx * pa_zzz * fz * fgb

                     + 6.0 * fx * fx * fx * fz * pb_z

                     + 22.5 * pa_yz * fx * fx * fz * pb_y

                     + 2.5 * fx * fx * pa_zzz * fz

                     + 15.0 * fx * fx * pa_zz * fz * pb_z

                     - 1.5 * pa_yz * fx * pb_y * fz * fgb

                     - 1.5 * pa_yz * fz * fga * pb_y * fx

                     - 3.0 * pa_y * fx * fz * fga * pb_yz

                     - 1.5 * fx * pa_z * fz * fga * pb_zz

                     - pa_yzzz * pb_y * fz * fgb

                     + 15.0 * pa_y * fx * fx * fz * pb_yz

                     + 7.5 * fx * fx * pa_z * fz * pb_zz

                     + 6.0 * pa_yzzz * fz * pb_y * fx

                     + 36.0 * pa_yzz * fz * fx * pb_yz

                     + 6.0 * fx * pa_zzz * fz * pb_zz

                     - 3.0 * pa_yz * fz * fga * pb_yzz

                     + 18.0 * pa_yz * fz * fx * pb_yzz

                     + 14.0 * pa_yzzz * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_zzz_s_0(double fx,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pa_yzzz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * pa_y * fx * fx * fx

                     + 2.25 * pa_yzz * fx * fx

                     + 6.75 * pa_yz * fx * fx * pb_z

                     + 2.25 * pa_y * fx * fx * pb_zz

                     + 1.5 * pa_yzzz * pb_z * fx

                     + 4.5 * pa_yzz * fx * pb_zz

                     + 1.5 * pa_yz * fx * pb_zzz

                     + pa_yzzz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_zzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pa_yzzz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (15.0 * pa_y * fx * fx * fx * fz

                     - 2.25 * pa_y * fx * fx * fz * fgb

                     - 2.25 * pa_y * fx * fx * fz * fga

                     - 4.5 * pa_yzz * fx * fz * fgb

                     + 22.5 * pa_yzz * fz * fx * fx

                     + 67.5 * pa_yz * fx * fx * fz * pb_z

                     - 4.5 * pa_yz * fx * pb_z * fz * fgb

                     - 4.5 * pa_yz * fz * fga * pb_z * fx

                     - 4.5 * pa_y * fx * fz * fga * pb_zz

                     - 3.0 * pa_yzzz * pb_z * fz * fgb

                     + 22.5 * pa_y * fx * fx * fz * pb_zz

                     + 18.0 * pa_yzzz * fz * pb_z * fx

                     + 54.0 * pa_yzz * fz * fx * pb_zz

                     - 3.0 * pa_yz * fz * fga * pb_zzz

                     + 18.0 * pa_yz * fz * fx * pb_zzz

                     + 14.0 * pa_yzzz * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xxx_s_0(double fx,
                                    double pa_zz,
                                    double pa_zzzz,
                                    double pb_x,
                                    double pb_xxx,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pb_x

                     + 4.5 * pa_zz * fx * fx * pb_x

                     + 1.5 * pa_zzzz * pb_x * fx

                     + 0.75 * fx * fx * pb_xxx

                     + 3.0 * pa_zz * fx * pb_xxx

                     + pa_zzzz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_zz,
                                    double pa_zzzz,
                                    double pb_x,
                                    double pb_xxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pb_x * fz * fgb

                     - 4.5 * fx * fx * fz * fga * pb_x

                     - 9.0 * pa_zz * fx * pb_x * fz * fgb

                     - 9.0 * pa_zz * fz * fga * pb_x * fx

                     + 9.0 * fx * fx * fx * fz * pb_x

                     - 3.0 * pa_zzzz * pb_x * fz * fgb

                     + 45.0 * pa_zz * fz * fx * fx * pb_x

                     + 18.0 * pa_zzzz * fz * pb_x * fx

                     - 3.0 * fx * fz * fga * pb_xxx

                     - 6.0 * pa_zz * fz * fga * pb_xxx

                     + 7.5 * fx * fx * fz * pb_xxx

                     + 36.0 * pa_zz * fz * fx * pb_xxx

                     + 14.0 * pa_zzzz * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xxy_s_0(double fx,
                                    double pa_zz,
                                    double pa_zzzz,
                                    double pb_xxy,
                                    double pb_y,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_y

                     + 1.5 * pa_zz * fx * fx * pb_y

                     + 0.5 * pa_zzzz * fx * pb_y

                     + 0.75 * fx * fx * pb_xxy

                     + 3.0 * pa_zz * fx * pb_xxy

                     + pa_zzzz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_zz,
                                    double pa_zzzz,
                                    double pb_xxy,
                                    double pb_y,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb * pb_y

                     - 1.5 * fx * fx * fz * fga * pb_y

                     - 3.0 * pa_zz * fx * fz * fgb * pb_y

                     - 3.0 * pa_zz * fz * fga * fx * pb_y

                     + 3.0 * fx * fx * fx * fz * pb_y

                     - pa_zzzz * fz * fgb * pb_y

                     + 15.0 * pa_zz * fz * fx * fx * pb_y

                     + 6.0 * pa_zzzz * fz * fx * pb_y

                     - 3.0 * fx * fz * fga * pb_xxy

                     - 6.0 * pa_zz * fz * fga * pb_xxy

                     + 7.5 * fx * fx * fz * pb_xxy

                     + 36.0 * pa_zz * fz * fx * pb_xxy

                     + 14.0 * pa_zzzz * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xxz_s_0(double fx,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pa_zzzz,
                                    double pb_xx,
                                    double pb_xxz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_z * fx * fx * fx

                     + pa_zzz * fx * fx

                     + 0.375 * fx * fx * fx * pb_z

                     + 1.5 * pa_zz * fx * fx * pb_z

                     + 3.0 * pa_z * fx * fx * pb_xx

                     + 0.5 * pa_zzzz * fx * pb_z

                     + 2.0 * pa_zzz * fx * pb_xx

                     + 0.75 * fx * fx * pb_xxz

                     + 3.0 * pa_zz * fx * pb_xxz

                     + pa_zzzz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pa_zzzz,
                                    double pb_xx,
                                    double pb_xxz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_z * fx * fx * fz * fgb

                     - 3.0 * pa_z * fx * fx * fz * fga

                     - 2.0 * pa_zzz * fx * fz * fgb

                     + 12.0 * pa_z * fx * fx * fx * fz

                     + 10.0 * pa_zzz * fz * fx * fx

                     - 0.75 * fx * fx * fz * fgb * pb_z

                     - 1.5 * fx * fx * fz * fga * pb_z

                     - 3.0 * pa_zz * fx * fz * fgb * pb_z

                     - 3.0 * pa_zz * fz * fga * fx * pb_z

                     - 6.0 * pa_z * fx * fz * fga * pb_xx

                     + 3.0 * fx * fx * fx * fz * pb_z

                     - pa_zzzz * fz * fgb * pb_z

                     + 15.0 * pa_zz * fz * fx * fx * pb_z

                     + 30.0 * pa_z * fx * fx * fz * pb_xx

                     + 6.0 * pa_zzzz * fz * fx * pb_z

                     + 24.0 * pa_zzz * fz * fx * pb_xx

                     - 3.0 * fx * fz * fga * pb_xxz

                     - 6.0 * pa_zz * fz * fga * pb_xxz

                     + 7.5 * fx * fx * fz * pb_xxz

                     + 36.0 * pa_zz * fz * fx * pb_xxz

                     + 14.0 * pa_zzzz * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xyy_s_0(double fx,
                                    double pa_zz,
                                    double pa_zzzz,
                                    double pb_x,
                                    double pb_xyy,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_x

                     + 1.5 * pa_zz * fx * fx * pb_x

                     + 0.5 * pa_zzzz * pb_x * fx

                     + 0.75 * fx * fx * pb_xyy

                     + 3.0 * pa_zz * fx * pb_xyy

                     + pa_zzzz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_zz,
                                    double pa_zzzz,
                                    double pb_x,
                                    double pb_xyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pb_x * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_x

                     - 3.0 * pa_zz * fx * pb_x * fz * fgb

                     - 3.0 * pa_zz * fz * fga * pb_x * fx

                     + 3.0 * fx * fx * fx * fz * pb_x

                     - pa_zzzz * pb_x * fz * fgb

                     + 15.0 * pa_zz * fz * fx * fx * pb_x

                     + 6.0 * pa_zzzz * fz * pb_x * fx

                     - 3.0 * fx * fz * fga * pb_xyy

                     - 6.0 * pa_zz * fz * fga * pb_xyy

                     + 7.5 * fx * fx * fz * pb_xyy

                     + 36.0 * pa_zz * fz * fx * pb_xyy

                     + 14.0 * pa_zzzz * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xyz_s_0(double fx,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pa_zzzz,
                                    double pb_xy,
                                    double pb_xyz,
                                    double s_0_0)
    {
        return s_0_0 * (3.0 * pa_z * fx * fx * pb_xy

                     + 2.0 * pa_zzz * fx * pb_xy

                     + 0.75 * fx * fx * pb_xyz

                     + 3.0 * pa_zz * fx * pb_xyz

                     + pa_zzzz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xyz_r_0(double fga,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pa_zzzz,
                                    double pb_xy,
                                    double pb_xyz,
                                    double r_0_0)
    {
        return r_0_0 * (- 6.0 * pa_z * fx * fz * fga * pb_xy

                     + 30.0 * pa_z * fx * fx * fz * pb_xy

                     + 24.0 * pa_zzz * fz * fx * pb_xy

                     - 3.0 * fx * fz * fga * pb_xyz

                     - 6.0 * pa_zz * fz * fga * pb_xyz

                     + 7.5 * fx * fx * fz * pb_xyz

                     + 36.0 * pa_zz * fz * fx * pb_xyz

                     + 14.0 * pa_zzzz * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xzz_s_0(double fx,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pa_zzzz,
                                    double pb_x,
                                    double pb_xz,
                                    double pb_xzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pb_x

                     + 4.5 * pa_zz * fx * fx * pb_x

                     + 6.0 * pa_z * fx * fx * pb_xz

                     + 0.5 * pa_zzzz * pb_x * fx

                     + 4.0 * pa_zzz * fx * pb_xz

                     + 0.75 * fx * fx * pb_xzz

                     + 3.0 * pa_zz * fx * pb_xzz

                     + pa_zzzz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pa_zzzz,
                                    double pb_x,
                                    double pb_xz,
                                    double pb_xzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * fx * fx * fz * fga * pb_x

                     + 15.0 * fx * fx * fx * fz * pb_x

                     + 45.0 * pa_zz * fx * fx * fz * pb_x

                     - 0.75 * fx * fx * pb_x * fz * fgb

                     - 3.0 * pa_zz * fx * pb_x * fz * fgb

                     - 3.0 * pa_zz * fz * fga * pb_x * fx

                     - 12.0 * pa_z * fx * fz * fga * pb_xz

                     - pa_zzzz * pb_x * fz * fgb

                     + 60.0 * pa_z * fx * fx * fz * pb_xz

                     + 6.0 * pa_zzzz * fz * pb_x * fx

                     + 48.0 * pa_zzz * fz * fx * pb_xz

                     - 3.0 * fx * fz * fga * pb_xzz

                     - 6.0 * pa_zz * fz * fga * pb_xzz

                     + 7.5 * fx * fx * fz * pb_xzz

                     + 36.0 * pa_zz * fz * fx * pb_xzz

                     + 14.0 * pa_zzzz * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_yyy_s_0(double fx,
                                    double pa_zz,
                                    double pa_zzzz,
                                    double pb_y,
                                    double pb_yyy,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pb_y

                     + 4.5 * pa_zz * fx * fx * pb_y

                     + 1.5 * pa_zzzz * pb_y * fx

                     + 0.75 * fx * fx * pb_yyy

                     + 3.0 * pa_zz * fx * pb_yyy

                     + pa_zzzz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_yyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_zz,
                                    double pa_zzzz,
                                    double pb_y,
                                    double pb_yyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pb_y * fz * fgb

                     - 4.5 * fx * fx * fz * fga * pb_y

                     - 9.0 * pa_zz * fx * pb_y * fz * fgb

                     - 9.0 * pa_zz * fz * fga * pb_y * fx

                     + 9.0 * fx * fx * fx * fz * pb_y

                     - 3.0 * pa_zzzz * pb_y * fz * fgb

                     + 45.0 * pa_zz * fz * fx * fx * pb_y

                     + 18.0 * pa_zzzz * fz * pb_y * fx

                     - 3.0 * fx * fz * fga * pb_yyy

                     - 6.0 * pa_zz * fz * fga * pb_yyy

                     + 7.5 * fx * fx * fz * pb_yyy

                     + 36.0 * pa_zz * fz * fx * pb_yyy

                     + 14.0 * pa_zzzz * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_yyz_s_0(double fx,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pa_zzzz,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_z * fx * fx * fx

                     + pa_zzz * fx * fx

                     + 0.375 * fx * fx * fx * pb_z

                     + 1.5 * pa_zz * fx * fx * pb_z

                     + 3.0 * pa_z * fx * fx * pb_yy

                     + 0.5 * pa_zzzz * fx * pb_z

                     + 2.0 * pa_zzz * fx * pb_yy

                     + 0.75 * fx * fx * pb_yyz

                     + 3.0 * pa_zz * fx * pb_yyz

                     + pa_zzzz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_yyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pa_zzzz,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_z * fx * fx * fz * fgb

                     - 3.0 * pa_z * fx * fx * fz * fga

                     - 2.0 * pa_zzz * fx * fz * fgb

                     + 12.0 * pa_z * fx * fx * fx * fz

                     + 10.0 * pa_zzz * fz * fx * fx

                     - 0.75 * fx * fx * fz * fgb * pb_z

                     - 1.5 * fx * fx * fz * fga * pb_z

                     - 3.0 * pa_zz * fx * fz * fgb * pb_z

                     - 3.0 * pa_zz * fz * fga * fx * pb_z

                     - 6.0 * pa_z * fx * fz * fga * pb_yy

                     + 3.0 * fx * fx * fx * fz * pb_z

                     - pa_zzzz * fz * fgb * pb_z

                     + 15.0 * pa_zz * fz * fx * fx * pb_z

                     + 30.0 * pa_z * fx * fx * fz * pb_yy

                     + 6.0 * pa_zzzz * fz * fx * pb_z

                     + 24.0 * pa_zzz * fz * fx * pb_yy

                     - 3.0 * fx * fz * fga * pb_yyz

                     - 6.0 * pa_zz * fz * fga * pb_yyz

                     + 7.5 * fx * fx * fz * pb_yyz

                     + 36.0 * pa_zz * fz * fx * pb_yyz

                     + 14.0 * pa_zzzz * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_yzz_s_0(double fx,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pa_zzzz,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pb_y

                     + 4.5 * pa_zz * fx * fx * pb_y

                     + 6.0 * pa_z * fx * fx * pb_yz

                     + 0.5 * pa_zzzz * pb_y * fx

                     + 4.0 * pa_zzz * fx * pb_yz

                     + 0.75 * fx * fx * pb_yzz

                     + 3.0 * pa_zz * fx * pb_yzz

                     + pa_zzzz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_yzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pa_zzzz,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * fx * fx * fz * fga * pb_y

                     + 15.0 * fx * fx * fx * fz * pb_y

                     + 45.0 * pa_zz * fx * fx * fz * pb_y

                     - 0.75 * fx * fx * pb_y * fz * fgb

                     - 3.0 * pa_zz * fx * pb_y * fz * fgb

                     - 3.0 * pa_zz * fz * fga * pb_y * fx

                     - 12.0 * pa_z * fx * fz * fga * pb_yz

                     - pa_zzzz * pb_y * fz * fgb

                     + 60.0 * pa_z * fx * fx * fz * pb_yz

                     + 6.0 * pa_zzzz * fz * pb_y * fx

                     + 48.0 * pa_zzz * fz * fx * pb_yz

                     - 3.0 * fx * fz * fga * pb_yzz

                     - 6.0 * pa_zz * fz * fga * pb_yzz

                     + 7.5 * fx * fx * fz * pb_yzz

                     + 36.0 * pa_zz * fz * fx * pb_yzz

                     + 14.0 * pa_zzzz * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_zzz_s_0(double fx,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pa_zzzz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (7.5 * pa_z * fx * fx * fx

                     + 5.625 * fx * fx * fx * pb_z

                     + 3.0 * pa_zzz * fx * fx

                     + 13.5 * pa_zz * fx * fx * pb_z

                     + 9.0 * pa_z * fx * fx * pb_zz

                     + 1.5 * pa_zzzz * pb_z * fx

                     + 6.0 * pa_zzz * fx * pb_zz

                     + 0.75 * fx * fx * pb_zzz

                     + 3.0 * pa_zz * fx * pb_zzz

                     + pa_zzzz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_zzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pa_zzzz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (60.0 * pa_z * fx * fx * fx * fz

                     - 9.0 * pa_z * fx * fx * fz * fgb

                     - 9.0 * pa_z * fx * fx * fz * fga

                     - 13.5 * fx * fx * fz * fga * pb_z

                     - 6.0 * pa_zzz * fx * fz * fgb

                     + 45.0 * fx * fx * fx * fz * pb_z

                     + 30.0 * pa_zzz * fz * fx * fx

                     + 135.0 * pa_zz * fx * fx * fz * pb_z

                     - 2.25 * fx * fx * pb_z * fz * fgb

                     - 9.0 * pa_zz * fx * pb_z * fz * fgb

                     - 9.0 * pa_zz * fz * fga * pb_z * fx

                     - 18.0 * pa_z * fx * fz * fga * pb_zz

                     - 3.0 * pa_zzzz * pb_z * fz * fgb

                     + 90.0 * pa_z * fx * fx * fz * pb_zz

                     + 18.0 * pa_zzzz * fz * pb_z * fx

                     + 72.0 * pa_zzz * fz * fx * pb_zz

                     - 3.0 * fx * fz * fga * pb_zzz

                     - 6.0 * pa_zz * fz * fga * pb_zzz

                     + 7.5 * fx * fx * fz * pb_zzz

                     + 36.0 * pa_zz * fz * fx * pb_zzz

                     + 14.0 * pa_zzzz * fz * pb_zzz);

    }


} // kinrecfunc namespace

