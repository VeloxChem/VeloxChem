//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

namespace kinvecfunc { // kinvecfunc namespace

    // SIMD elementary functions for (G|T|G) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xxxx_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxx,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxx,
                                     double s_0_0)
    {
        return s_0_0 * (6.5625 * fx * fx * fx * fx

                     + 11.25 * pa_xx * fx * fx * fx

                     + 30.0 * pa_x * fx * fx * fx * pb_x

                     + 11.25 * fx * fx * fx * pb_xx

                     + 0.75 * pa_xxxx * fx * fx

                     + 12.0 * pa_xxx * fx * fx * pb_x

                     + 27.0 * pa_xx * fx * fx * pb_xx

                     + 12.0 * pa_x * fx * fx * pb_xxx

                     + 3.0 * pa_xxxx * pb_xx * fx

                     + 8.0 * pa_xxx * fx * pb_xxx

                     + 0.75 * fx * fx * pb_xxxx

                     + 3.0 * pa_xx * fx * pb_xxxx

                     + pa_xxxx * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xxxx_r_0(double fga,
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
                                     double pb_xxxx,
                                     double r_0_0)
    {
        return r_0_0 * (52.5 * fx * fx * fx * fx * fz

                     - 11.25 * fx * fx * fx * fz * fgb

                     - 11.25 * fx * fx * fx * fz * fga

                     - 27.0 * pa_xx * fx * fx * fz * fgb

                     + 112.5 * pa_xx * fx * fx * fx * fz

                     + 300.0 * pa_x * fx * fx * fx * fz * pb_x

                     - 4.5 * pa_xx * fz * fga * fx * fx

                     - 36.0 * pa_x * fx * fx * pb_x * fz * fgb

                     - 36.0 * pa_x * fx * fx * fz * fga * pb_x

                     - 27.0 * fx * fx * fz * fga * pb_xx

                     - 3.0 * pa_xxxx * fx * fz * fgb

                     - 24.0 * pa_xxx * fx * pb_x * fz * fgb

                     + 112.5 * fx * fx * fx * fz * pb_xx

                     + 9.0 * pa_xxxx * fz * fx * fx

                     + 144.0 * pa_xxx * fz * fx * fx * pb_x

                     + 324.0 * pa_xx * fx * fx * fz * pb_xx

                     - 4.5 * fx * fx * pb_xx * fz * fgb

                     - 18.0 * pa_xx * fx * pb_xx * fz * fgb

                     - 18.0 * pa_xx * fz * fga * pb_xx * fx

                     - 24.0 * pa_x * fx * fz * fga * pb_xxx

                     - 6.0 * pa_xxxx * pb_xx * fz * fgb

                     + 144.0 * pa_x * fx * fx * fz * pb_xxx

                     + 42.0 * pa_xxxx * fz * pb_xx * fx

                     + 112.0 * pa_xxx * fz * fx * pb_xxx

                     - 3.0 * fx * fz * fga * pb_xxxx

                     - 6.0 * pa_xx * fz * fga * pb_xxxx

                     + 9.0 * fx * fx * fz * pb_xxxx

                     + 42.0 * pa_xx * fz * fx * pb_xxxx

                     + 16.0 * pa_xxxx * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xxxy_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxx,
                                     double pb_xxxy,
                                     double pb_xxy,
                                     double pb_xy,
                                     double pb_y,
                                     double s_0_0)
    {
        return s_0_0 * (7.5 * pa_x * fx * fx * fx * pb_y

                     + 5.625 * fx * fx * fx * pb_xy

                     + 3.0 * pa_xxx * fx * fx * pb_y

                     + 13.5 * pa_xx * fx * fx * pb_xy

                     + 9.0 * pa_x * fx * fx * pb_xxy

                     + 1.5 * pa_xxxx * pb_xy * fx

                     + 6.0 * pa_xxx * fx * pb_xxy

                     + 0.75 * fx * fx * pb_xxxy

                     + 3.0 * pa_xx * fx * pb_xxxy

                     + pa_xxxx * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xxxy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxx,
                                     double pb_xxxy,
                                     double pb_xxy,
                                     double pb_xy,
                                     double pb_y,
                                     double r_0_0)
    {
        return r_0_0 * (75.0 * pa_x * fx * fx * fx * fz * pb_y

                     - 9.0 * pa_x * fx * fx * fz * fgb * pb_y

                     - 9.0 * pa_x * fx * fx * fz * fga * pb_y

                     - 13.5 * fx * fx * fz * fga * pb_xy

                     - 6.0 * pa_xxx * fx * fz * fgb * pb_y

                     + 56.25 * fx * fx * fx * fz * pb_xy

                     + 36.0 * pa_xxx * fz * fx * fx * pb_y

                     + 162.0 * pa_xx * fx * fx * fz * pb_xy

                     - 2.25 * fx * fx * pb_xy * fz * fgb

                     - 9.0 * pa_xx * fx * pb_xy * fz * fgb

                     - 9.0 * pa_xx * fz * fga * pb_xy * fx

                     - 18.0 * pa_x * fx * fz * fga * pb_xxy

                     - 3.0 * pa_xxxx * pb_xy * fz * fgb

                     + 108.0 * pa_x * fx * fx * fz * pb_xxy

                     + 21.0 * pa_xxxx * fz * pb_xy * fx

                     + 84.0 * pa_xxx * fz * fx * pb_xxy

                     - 3.0 * fx * fz * fga * pb_xxxy

                     - 6.0 * pa_xx * fz * fga * pb_xxxy

                     + 9.0 * fx * fx * fz * pb_xxxy

                     + 42.0 * pa_xx * fz * fx * pb_xxxy

                     + 16.0 * pa_xxxx * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xxxz_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxx,
                                     double pb_xxxz,
                                     double pb_xxz,
                                     double pb_xz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (7.5 * pa_x * fx * fx * fx * pb_z

                     + 5.625 * fx * fx * fx * pb_xz

                     + 3.0 * pa_xxx * fx * fx * pb_z

                     + 13.5 * pa_xx * fx * fx * pb_xz

                     + 9.0 * pa_x * fx * fx * pb_xxz

                     + 1.5 * pa_xxxx * pb_xz * fx

                     + 6.0 * pa_xxx * fx * pb_xxz

                     + 0.75 * fx * fx * pb_xxxz

                     + 3.0 * pa_xx * fx * pb_xxxz

                     + pa_xxxx * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xxxz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxx,
                                     double pb_xxxz,
                                     double pb_xxz,
                                     double pb_xz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (75.0 * pa_x * fx * fx * fx * fz * pb_z

                     - 9.0 * pa_x * fx * fx * fz * fgb * pb_z

                     - 9.0 * pa_x * fx * fx * fz * fga * pb_z

                     - 13.5 * fx * fx * fz * fga * pb_xz

                     - 6.0 * pa_xxx * fx * fz * fgb * pb_z

                     + 56.25 * fx * fx * fx * fz * pb_xz

                     + 36.0 * pa_xxx * fz * fx * fx * pb_z

                     + 162.0 * pa_xx * fx * fx * fz * pb_xz

                     - 2.25 * fx * fx * pb_xz * fz * fgb

                     - 9.0 * pa_xx * fx * pb_xz * fz * fgb

                     - 9.0 * pa_xx * fz * fga * pb_xz * fx

                     - 18.0 * pa_x * fx * fz * fga * pb_xxz

                     - 3.0 * pa_xxxx * pb_xz * fz * fgb

                     + 108.0 * pa_x * fx * fx * fz * pb_xxz

                     + 21.0 * pa_xxxx * fz * pb_xz * fx

                     + 84.0 * pa_xxx * fz * fx * pb_xxz

                     - 3.0 * fx * fz * fga * pb_xxxz

                     - 6.0 * pa_xx * fz * fga * pb_xxxz

                     + 9.0 * fx * fx * fz * pb_xxxz

                     + 42.0 * pa_xx * fz * fx * pb_xxxz

                     + 16.0 * pa_xxxx * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xxyy_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxx,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxyy,
                                     double pb_xyy,
                                     double pb_yy,
                                     double s_0_0)
    {
        return s_0_0 * (0.9375 * fx * fx * fx * fx

                     + 2.25 * pa_xx * fx * fx * fx

                     + 3.0 * pa_x * fx * fx * fx * pb_x

                     + 1.875 * fx * fx * fx * pb_yy

                     + 0.25 * pa_xxxx * fx * fx

                     + 2.0 * pa_xxx * fx * fx * pb_x

                     + 4.5 * pa_xx * fx * fx * pb_yy

                     + 0.375 * fx * fx * fx * pb_xx

                     + 1.5 * pa_xx * fx * fx * pb_xx

                     + 6.0 * pa_x * fx * fx * pb_xyy

                     + 0.5 * pa_xxxx * pb_xx * fx

                     + 0.5 * pa_xxxx * fx * pb_yy

                     + 4.0 * pa_xxx * fx * pb_xyy

                     + 0.75 * fx * fx * pb_xxyy

                     + 3.0 * pa_xx * fx * pb_xxyy

                     + pa_xxxx * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xxyy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxx,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxyy,
                                     double pb_xyy,
                                     double pb_yy,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * fx * fz * fgb

                     - 2.25 * fx * fx * fx * fz * fga

                     - 6.0 * pa_xx * fx * fx * fz * fgb

                     + 7.5 * fx * fx * fx * fx * fz

                     + 22.5 * pa_xx * fx * fx * fx * fz

                     - 1.5 * pa_xx * fz * fga * fx * fx

                     - 6.0 * pa_x * fx * fx * pb_x * fz * fgb

                     - 6.0 * pa_x * fx * fx * fz * fga * pb_x

                     - 4.5 * fx * fx * fz * fga * pb_yy

                     - pa_xxxx * fx * fz * fgb

                     - 4.0 * pa_xxx * fx * pb_x * fz * fgb

                     + 30.0 * pa_x * fx * fx * fx * fz * pb_x

                     + 18.75 * fx * fx * fx * fz * pb_yy

                     + 3.0 * pa_xxxx * fz * fx * fx

                     + 24.0 * pa_xxx * fz * fx * fx * pb_x

                     + 54.0 * pa_xx * fx * fx * fz * pb_yy

                     - 0.75 * fx * fx * pb_xx * fz * fgb

                     - 0.75 * fx * fx * fz * fgb * pb_yy

                     - 1.5 * fx * fx * fz * fga * pb_xx

                     - 3.0 * pa_xx * fx * pb_xx * fz * fgb

                     - 3.0 * pa_xx * fx * fz * fgb * pb_yy

                     - 3.0 * pa_xx * fz * fga * pb_xx * fx

                     - 3.0 * pa_xx * fz * fga * fx * pb_yy

                     - 12.0 * pa_x * fx * fz * fga * pb_xyy

                     + 3.75 * fx * fx * fx * fz * pb_xx

                     - pa_xxxx * pb_xx * fz * fgb

                     - pa_xxxx * fz * fgb * pb_yy

                     + 18.0 * pa_xx * fz * fx * fx * pb_xx

                     + 72.0 * pa_x * fx * fx * fz * pb_xyy

                     + 7.0 * pa_xxxx * fz * pb_xx * fx

                     + 7.0 * pa_xxxx * fz * fx * pb_yy

                     + 56.0 * pa_xxx * fz * fx * pb_xyy

                     - 3.0 * fx * fz * fga * pb_xxyy

                     - 6.0 * pa_xx * fz * fga * pb_xxyy

                     + 9.0 * fx * fx * fz * pb_xxyy

                     + 42.0 * pa_xx * fz * fx * pb_xxyy

                     + 16.0 * pa_xxxx * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xxyz_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxx,
                                     double pb_xxyz,
                                     double pb_xyz,
                                     double pb_yz,
                                     double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pb_yz

                     + 4.5 * pa_xx * fx * fx * pb_yz

                     + 6.0 * pa_x * fx * fx * pb_xyz

                     + 0.5 * pa_xxxx * fx * pb_yz

                     + 4.0 * pa_xxx * fx * pb_xyz

                     + 0.75 * fx * fx * pb_xxyz

                     + 3.0 * pa_xx * fx * pb_xxyz

                     + pa_xxxx * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xxyz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxx,
                                     double pb_xxyz,
                                     double pb_xyz,
                                     double pb_yz,
                                     double r_0_0)
    {
        return r_0_0 * (- 4.5 * fx * fx * fz * fga * pb_yz

                     + 18.75 * fx * fx * fx * fz * pb_yz

                     + 54.0 * pa_xx * fx * fx * fz * pb_yz

                     - 0.75 * fx * fx * fz * fgb * pb_yz

                     - 3.0 * pa_xx * fx * fz * fgb * pb_yz

                     - 3.0 * pa_xx * fz * fga * fx * pb_yz

                     - 12.0 * pa_x * fx * fz * fga * pb_xyz

                     - pa_xxxx * fz * fgb * pb_yz

                     + 72.0 * pa_x * fx * fx * fz * pb_xyz

                     + 7.0 * pa_xxxx * fz * fx * pb_yz

                     + 56.0 * pa_xxx * fz * fx * pb_xyz

                     - 3.0 * fx * fz * fga * pb_xxyz

                     - 6.0 * pa_xx * fz * fga * pb_xxyz

                     + 9.0 * fx * fx * fz * pb_xxyz

                     + 42.0 * pa_xx * fz * fx * pb_xxyz

                     + 16.0 * pa_xxxx * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xxzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxx,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxzz,
                                     double pb_xzz,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.9375 * fx * fx * fx * fx

                     + 2.25 * pa_xx * fx * fx * fx

                     + 3.0 * pa_x * fx * fx * fx * pb_x

                     + 1.875 * fx * fx * fx * pb_zz

                     + 0.25 * pa_xxxx * fx * fx

                     + 2.0 * pa_xxx * fx * fx * pb_x

                     + 4.5 * pa_xx * fx * fx * pb_zz

                     + 0.375 * fx * fx * fx * pb_xx

                     + 1.5 * pa_xx * fx * fx * pb_xx

                     + 6.0 * pa_x * fx * fx * pb_xzz

                     + 0.5 * pa_xxxx * pb_xx * fx

                     + 0.5 * pa_xxxx * fx * pb_zz

                     + 4.0 * pa_xxx * fx * pb_xzz

                     + 0.75 * fx * fx * pb_xxzz

                     + 3.0 * pa_xx * fx * pb_xxzz

                     + pa_xxxx * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xxzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxx,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxzz,
                                     double pb_xzz,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * fx * fz * fgb

                     - 2.25 * fx * fx * fx * fz * fga

                     - 6.0 * pa_xx * fx * fx * fz * fgb

                     + 7.5 * fx * fx * fx * fx * fz

                     + 22.5 * pa_xx * fx * fx * fx * fz

                     - 1.5 * pa_xx * fz * fga * fx * fx

                     - 6.0 * pa_x * fx * fx * pb_x * fz * fgb

                     - 6.0 * pa_x * fx * fx * fz * fga * pb_x

                     - 4.5 * fx * fx * fz * fga * pb_zz

                     - pa_xxxx * fx * fz * fgb

                     - 4.0 * pa_xxx * fx * pb_x * fz * fgb

                     + 30.0 * pa_x * fx * fx * fx * fz * pb_x

                     + 18.75 * fx * fx * fx * fz * pb_zz

                     + 3.0 * pa_xxxx * fz * fx * fx

                     + 24.0 * pa_xxx * fz * fx * fx * pb_x

                     + 54.0 * pa_xx * fx * fx * fz * pb_zz

                     - 0.75 * fx * fx * pb_xx * fz * fgb

                     - 0.75 * fx * fx * fz * fgb * pb_zz

                     - 1.5 * fx * fx * fz * fga * pb_xx

                     - 3.0 * pa_xx * fx * pb_xx * fz * fgb

                     - 3.0 * pa_xx * fx * fz * fgb * pb_zz

                     - 3.0 * pa_xx * fz * fga * pb_xx * fx

                     - 3.0 * pa_xx * fz * fga * fx * pb_zz

                     - 12.0 * pa_x * fx * fz * fga * pb_xzz

                     + 3.75 * fx * fx * fx * fz * pb_xx

                     - pa_xxxx * pb_xx * fz * fgb

                     - pa_xxxx * fz * fgb * pb_zz

                     + 18.0 * pa_xx * fz * fx * fx * pb_xx

                     + 72.0 * pa_x * fx * fx * fz * pb_xzz

                     + 7.0 * pa_xxxx * fz * pb_xx * fx

                     + 7.0 * pa_xxxx * fz * fx * pb_zz

                     + 56.0 * pa_xxx * fz * fx * pb_xzz

                     - 3.0 * fx * fz * fga * pb_xxzz

                     - 6.0 * pa_xx * fz * fga * pb_xxzz

                     + 9.0 * fx * fx * fz * pb_xxzz

                     + 42.0 * pa_xx * fz * fx * pb_xxzz

                     + 16.0 * pa_xxxx * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xyyy_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxx,
                                     double pb_xy,
                                     double pb_xyyy,
                                     double pb_y,
                                     double pb_yyy,
                                     double s_0_0)
    {
        return s_0_0 * (4.5 * pa_x * fx * fx * fx * pb_y

                     + 3.0 * pa_xxx * fx * fx * pb_y

                     + 1.125 * fx * fx * fx * pb_xy

                     + 4.5 * pa_xx * fx * fx * pb_xy

                     + 3.0 * pa_x * fx * fx * pb_yyy

                     + 1.5 * pa_xxxx * pb_xy * fx

                     + 2.0 * pa_xxx * fx * pb_yyy

                     + 0.75 * fx * fx * pb_xyyy

                     + 3.0 * pa_xx * fx * pb_xyyy

                     + pa_xxxx * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xyyy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxx,
                                     double pb_xy,
                                     double pb_xyyy,
                                     double pb_y,
                                     double pb_yyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 9.0 * pa_x * fx * fx * pb_y * fz * fgb

                     - 9.0 * pa_x * fx * fx * fz * fga * pb_y

                     - 6.0 * pa_xxx * fx * pb_y * fz * fgb

                     + 45.0 * pa_x * fx * fx * fx * fz * pb_y

                     + 36.0 * pa_xxx * fz * fx * fx * pb_y

                     - 2.25 * fx * fx * pb_xy * fz * fgb

                     - 4.5 * fx * fx * fz * fga * pb_xy

                     - 9.0 * pa_xx * fx * pb_xy * fz * fgb

                     - 9.0 * pa_xx * fz * fga * pb_xy * fx

                     - 6.0 * pa_x * fx * fz * fga * pb_yyy

                     + 11.25 * fx * fx * fx * fz * pb_xy

                     - 3.0 * pa_xxxx * pb_xy * fz * fgb

                     + 54.0 * pa_xx * fz * fx * fx * pb_xy

                     + 36.0 * pa_x * fx * fx * fz * pb_yyy

                     + 21.0 * pa_xxxx * fz * pb_xy * fx

                     + 28.0 * pa_xxx * fz * fx * pb_yyy

                     - 3.0 * fx * fz * fga * pb_xyyy

                     - 6.0 * pa_xx * fz * fga * pb_xyyy

                     + 9.0 * fx * fx * fz * pb_xyyy

                     + 42.0 * pa_xx * fz * fx * pb_xyyy

                     + 16.0 * pa_xxxx * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xyyz_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxx,
                                     double pb_xyyz,
                                     double pb_xz,
                                     double pb_yyz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (1.5 * pa_x * fx * fx * fx * pb_z

                     + pa_xxx * fx * fx * pb_z

                     + 0.375 * fx * fx * fx * pb_xz

                     + 1.5 * pa_xx * fx * fx * pb_xz

                     + 3.0 * pa_x * fx * fx * pb_yyz

                     + 0.5 * pa_xxxx * pb_xz * fx

                     + 2.0 * pa_xxx * fx * pb_yyz

                     + 0.75 * fx * fx * pb_xyyz

                     + 3.0 * pa_xx * fx * pb_xyyz

                     + pa_xxxx * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xyyz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxx,
                                     double pb_xyyz,
                                     double pb_xz,
                                     double pb_yyz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_x * fx * fx * fz * fgb * pb_z

                     - 3.0 * pa_x * fx * fx * fz * fga * pb_z

                     - 2.0 * pa_xxx * fx * fz * fgb * pb_z

                     + 15.0 * pa_x * fx * fx * fx * fz * pb_z

                     + 12.0 * pa_xxx * fz * fx * fx * pb_z

                     - 0.75 * fx * fx * pb_xz * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_xz

                     - 3.0 * pa_xx * fx * pb_xz * fz * fgb

                     - 3.0 * pa_xx * fz * fga * pb_xz * fx

                     - 6.0 * pa_x * fx * fz * fga * pb_yyz

                     + 3.75 * fx * fx * fx * fz * pb_xz

                     - pa_xxxx * pb_xz * fz * fgb

                     + 18.0 * pa_xx * fz * fx * fx * pb_xz

                     + 36.0 * pa_x * fx * fx * fz * pb_yyz

                     + 7.0 * pa_xxxx * fz * pb_xz * fx

                     + 28.0 * pa_xxx * fz * fx * pb_yyz

                     - 3.0 * fx * fz * fga * pb_xyyz

                     - 6.0 * pa_xx * fz * fga * pb_xyyz

                     + 9.0 * fx * fx * fz * pb_xyyz

                     + 42.0 * pa_xx * fz * fx * pb_xyyz

                     + 16.0 * pa_xxxx * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xyzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxx,
                                     double pb_xy,
                                     double pb_xyzz,
                                     double pb_y,
                                     double pb_yzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.5 * pa_x * fx * fx * fx * pb_y

                     + pa_xxx * fx * fx * pb_y

                     + 0.375 * fx * fx * fx * pb_xy

                     + 1.5 * pa_xx * fx * fx * pb_xy

                     + 3.0 * pa_x * fx * fx * pb_yzz

                     + 0.5 * pa_xxxx * pb_xy * fx

                     + 2.0 * pa_xxx * fx * pb_yzz

                     + 0.75 * fx * fx * pb_xyzz

                     + 3.0 * pa_xx * fx * pb_xyzz

                     + pa_xxxx * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xyzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxx,
                                     double pb_xy,
                                     double pb_xyzz,
                                     double pb_y,
                                     double pb_yzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_x * fx * fx * pb_y * fz * fgb

                     - 3.0 * pa_x * fx * fx * fz * fga * pb_y

                     - 2.0 * pa_xxx * fx * pb_y * fz * fgb

                     + 15.0 * pa_x * fx * fx * fx * fz * pb_y

                     + 12.0 * pa_xxx * fz * fx * fx * pb_y

                     - 0.75 * fx * fx * pb_xy * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_xy

                     - 3.0 * pa_xx * fx * pb_xy * fz * fgb

                     - 3.0 * pa_xx * fz * fga * pb_xy * fx

                     - 6.0 * pa_x * fx * fz * fga * pb_yzz

                     + 3.75 * fx * fx * fx * fz * pb_xy

                     - pa_xxxx * pb_xy * fz * fgb

                     + 18.0 * pa_xx * fz * fx * fx * pb_xy

                     + 36.0 * pa_x * fx * fx * fz * pb_yzz

                     + 7.0 * pa_xxxx * fz * pb_xy * fx

                     + 28.0 * pa_xxx * fz * fx * pb_yzz

                     - 3.0 * fx * fz * fga * pb_xyzz

                     - 6.0 * pa_xx * fz * fga * pb_xyzz

                     + 9.0 * fx * fx * fz * pb_xyzz

                     + 42.0 * pa_xx * fz * fx * pb_xyzz

                     + 16.0 * pa_xxxx * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xzzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxx,
                                     double pb_xz,
                                     double pb_xzzz,
                                     double pb_z,
                                     double pb_zzz,
                                     double s_0_0)
    {
        return s_0_0 * (4.5 * pa_x * fx * fx * fx * pb_z

                     + 3.0 * pa_xxx * fx * fx * pb_z

                     + 1.125 * fx * fx * fx * pb_xz

                     + 4.5 * pa_xx * fx * fx * pb_xz

                     + 3.0 * pa_x * fx * fx * pb_zzz

                     + 1.5 * pa_xxxx * pb_xz * fx

                     + 2.0 * pa_xxx * fx * pb_zzz

                     + 0.75 * fx * fx * pb_xzzz

                     + 3.0 * pa_xx * fx * pb_xzzz

                     + pa_xxxx * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xzzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxx,
                                     double pb_xz,
                                     double pb_xzzz,
                                     double pb_z,
                                     double pb_zzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 9.0 * pa_x * fx * fx * pb_z * fz * fgb

                     - 9.0 * pa_x * fx * fx * fz * fga * pb_z

                     - 6.0 * pa_xxx * fx * pb_z * fz * fgb

                     + 45.0 * pa_x * fx * fx * fx * fz * pb_z

                     + 36.0 * pa_xxx * fz * fx * fx * pb_z

                     - 2.25 * fx * fx * pb_xz * fz * fgb

                     - 4.5 * fx * fx * fz * fga * pb_xz

                     - 9.0 * pa_xx * fx * pb_xz * fz * fgb

                     - 9.0 * pa_xx * fz * fga * pb_xz * fx

                     - 6.0 * pa_x * fx * fz * fga * pb_zzz

                     + 11.25 * fx * fx * fx * fz * pb_xz

                     - 3.0 * pa_xxxx * pb_xz * fz * fgb

                     + 54.0 * pa_xx * fz * fx * fx * pb_xz

                     + 36.0 * pa_x * fx * fx * fz * pb_zzz

                     + 21.0 * pa_xxxx * fz * pb_xz * fx

                     + 28.0 * pa_xxx * fz * fx * pb_zzz

                     - 3.0 * fx * fz * fga * pb_xzzz

                     - 6.0 * pa_xx * fz * fga * pb_xzzz

                     + 9.0 * fx * fx * fz * pb_xzzz

                     + 42.0 * pa_xx * fz * fx * pb_xzzz

                     + 16.0 * pa_xxxx * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_yyyy_s_0(double fx,
                                     double pa_xx,
                                     double pa_xxxx,
                                     double pb_yy,
                                     double pb_yyyy,
                                     double s_0_0)
    {
        return s_0_0 * (0.5625 * fx * fx * fx * fx

                     + 2.25 * pa_xx * fx * fx * fx

                     + 0.75 * pa_xxxx * fx * fx

                     + 2.25 * fx * fx * fx * pb_yy

                     + 9.0 * pa_xx * fx * fx * pb_yy

                     + 3.0 * pa_xxxx * pb_yy * fx

                     + 0.75 * fx * fx * pb_yyyy

                     + 3.0 * pa_xx * fx * pb_yyyy

                     + pa_xxxx * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_yyyy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xx,
                                     double pa_xxxx,
                                     double pb_yy,
                                     double pb_yyyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * fx * fz * fgb

                     - 2.25 * fx * fx * fx * fz * fga

                     - 9.0 * pa_xx * fx * fx * fz * fgb

                     - 4.5 * pa_xx * fz * fga * fx * fx

                     + 4.5 * fx * fx * fx * fx * fz

                     - 3.0 * pa_xxxx * fx * fz * fgb

                     + 22.5 * pa_xx * fz * fx * fx * fx

                     + 9.0 * pa_xxxx * fz * fx * fx

                     - 4.5 * fx * fx * pb_yy * fz * fgb

                     - 9.0 * fx * fx * fz * fga * pb_yy

                     - 18.0 * pa_xx * fx * pb_yy * fz * fgb

                     - 18.0 * pa_xx * fz * fga * pb_yy * fx

                     + 22.5 * fx * fx * fx * fz * pb_yy

                     - 6.0 * pa_xxxx * pb_yy * fz * fgb

                     + 108.0 * pa_xx * fz * fx * fx * pb_yy

                     + 42.0 * pa_xxxx * fz * pb_yy * fx

                     - 3.0 * fx * fz * fga * pb_yyyy

                     - 6.0 * pa_xx * fz * fga * pb_yyyy

                     + 9.0 * fx * fx * fz * pb_yyyy

                     + 42.0 * pa_xx * fz * fx * pb_yyyy

                     + 16.0 * pa_xxxx * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_yyyz_s_0(double fx,
                                     double pa_xx,
                                     double pa_xxxx,
                                     double pb_yyyz,
                                     double pb_yz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pb_yz

                     + 4.5 * pa_xx * fx * fx * pb_yz

                     + 1.5 * pa_xxxx * pb_yz * fx

                     + 0.75 * fx * fx * pb_yyyz

                     + 3.0 * pa_xx * fx * pb_yyyz

                     + pa_xxxx * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_yyyz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xx,
                                     double pa_xxxx,
                                     double pb_yyyz,
                                     double pb_yz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pb_yz * fz * fgb

                     - 4.5 * fx * fx * fz * fga * pb_yz

                     - 9.0 * pa_xx * fx * pb_yz * fz * fgb

                     - 9.0 * pa_xx * fz * fga * pb_yz * fx

                     + 11.25 * fx * fx * fx * fz * pb_yz

                     - 3.0 * pa_xxxx * pb_yz * fz * fgb

                     + 54.0 * pa_xx * fz * fx * fx * pb_yz

                     + 21.0 * pa_xxxx * fz * pb_yz * fx

                     - 3.0 * fx * fz * fga * pb_yyyz

                     - 6.0 * pa_xx * fz * fga * pb_yyyz

                     + 9.0 * fx * fx * fz * pb_yyyz

                     + 42.0 * pa_xx * fz * fx * pb_yyyz

                     + 16.0 * pa_xxxx * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_yyzz_s_0(double fx,
                                     double pa_xx,
                                     double pa_xxxx,
                                     double pb_yy,
                                     double pb_yyzz,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.75 * pa_xx * fx * fx * fx

                     + 0.25 * pa_xxxx * fx * fx

                     + 0.375 * fx * fx * fx * pb_yy

                     + 0.375 * fx * fx * fx * pb_zz

                     + 1.5 * pa_xx * fx * fx * pb_yy

                     + 1.5 * pa_xx * fx * fx * pb_zz

                     + 0.5 * pa_xxxx * pb_yy * fx

                     + 0.5 * pa_xxxx * fx * pb_zz

                     + 0.75 * fx * fx * pb_yyzz

                     + 3.0 * pa_xx * fx * pb_yyzz

                     + pa_xxxx * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_yyzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xx,
                                     double pa_xxxx,
                                     double pb_yy,
                                     double pb_yyzz,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fx * fz * fgb

                     - 0.75 * fx * fx * fx * fz * fga

                     - 3.0 * pa_xx * fx * fx * fz * fgb

                     - 1.5 * pa_xx * fz * fga * fx * fx

                     + 1.5 * fx * fx * fx * fx * fz

                     - pa_xxxx * fx * fz * fgb

                     + 7.5 * pa_xx * fz * fx * fx * fx

                     + 3.0 * pa_xxxx * fz * fx * fx

                     - 0.75 * fx * fx * pb_yy * fz * fgb

                     - 0.75 * fx * fx * fz * fgb * pb_zz

                     - 1.5 * fx * fx * fz * fga * pb_yy

                     - 1.5 * fx * fx * fz * fga * pb_zz

                     - 3.0 * pa_xx * fx * pb_yy * fz * fgb

                     - 3.0 * pa_xx * fx * fz * fgb * pb_zz

                     - 3.0 * pa_xx * fz * fga * pb_yy * fx

                     - 3.0 * pa_xx * fz * fga * fx * pb_zz

                     + 3.75 * fx * fx * fx * fz * pb_yy

                     + 3.75 * fx * fx * fx * fz * pb_zz

                     - pa_xxxx * pb_yy * fz * fgb

                     - pa_xxxx * fz * fgb * pb_zz

                     + 18.0 * pa_xx * fz * fx * fx * pb_yy

                     + 18.0 * pa_xx * fz * fx * fx * pb_zz

                     + 7.0 * pa_xxxx * fz * pb_yy * fx

                     + 7.0 * pa_xxxx * fz * fx * pb_zz

                     - 3.0 * fx * fz * fga * pb_yyzz

                     - 6.0 * pa_xx * fz * fga * pb_yyzz

                     + 9.0 * fx * fx * fz * pb_yyzz

                     + 42.0 * pa_xx * fz * fx * pb_yyzz

                     + 16.0 * pa_xxxx * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_yzzz_s_0(double fx,
                                     double pa_xx,
                                     double pa_xxxx,
                                     double pb_yz,
                                     double pb_yzzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pb_yz

                     + 4.5 * pa_xx * fx * fx * pb_yz

                     + 1.5 * pa_xxxx * pb_yz * fx

                     + 0.75 * fx * fx * pb_yzzz

                     + 3.0 * pa_xx * fx * pb_yzzz

                     + pa_xxxx * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_yzzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xx,
                                     double pa_xxxx,
                                     double pb_yz,
                                     double pb_yzzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pb_yz * fz * fgb

                     - 4.5 * fx * fx * fz * fga * pb_yz

                     - 9.0 * pa_xx * fx * pb_yz * fz * fgb

                     - 9.0 * pa_xx * fz * fga * pb_yz * fx

                     + 11.25 * fx * fx * fx * fz * pb_yz

                     - 3.0 * pa_xxxx * pb_yz * fz * fgb

                     + 54.0 * pa_xx * fz * fx * fx * pb_yz

                     + 21.0 * pa_xxxx * fz * pb_yz * fx

                     - 3.0 * fx * fz * fga * pb_yzzz

                     - 6.0 * pa_xx * fz * fga * pb_yzzz

                     + 9.0 * fx * fx * fz * pb_yzzz

                     + 42.0 * pa_xx * fz * fx * pb_yzzz

                     + 16.0 * pa_xxxx * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_zzzz_s_0(double fx,
                                     double pa_xx,
                                     double pa_xxxx,
                                     double pb_zz,
                                     double pb_zzzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.5625 * fx * fx * fx * fx

                     + 2.25 * pa_xx * fx * fx * fx

                     + 0.75 * pa_xxxx * fx * fx

                     + 2.25 * fx * fx * fx * pb_zz

                     + 9.0 * pa_xx * fx * fx * pb_zz

                     + 3.0 * pa_xxxx * pb_zz * fx

                     + 0.75 * fx * fx * pb_zzzz

                     + 3.0 * pa_xx * fx * pb_zzzz

                     + pa_xxxx * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_zzzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xx,
                                     double pa_xxxx,
                                     double pb_zz,
                                     double pb_zzzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * fx * fz * fgb

                     - 2.25 * fx * fx * fx * fz * fga

                     - 9.0 * pa_xx * fx * fx * fz * fgb

                     - 4.5 * pa_xx * fz * fga * fx * fx

                     + 4.5 * fx * fx * fx * fx * fz

                     - 3.0 * pa_xxxx * fx * fz * fgb

                     + 22.5 * pa_xx * fz * fx * fx * fx

                     + 9.0 * pa_xxxx * fz * fx * fx

                     - 4.5 * fx * fx * pb_zz * fz * fgb

                     - 9.0 * fx * fx * fz * fga * pb_zz

                     - 18.0 * pa_xx * fx * pb_zz * fz * fgb

                     - 18.0 * pa_xx * fz * fga * pb_zz * fx

                     + 22.5 * fx * fx * fx * fz * pb_zz

                     - 6.0 * pa_xxxx * pb_zz * fz * fgb

                     + 108.0 * pa_xx * fz * fx * fx * pb_zz

                     + 42.0 * pa_xxxx * fz * pb_zz * fx

                     - 3.0 * fx * fz * fga * pb_zzzz

                     - 6.0 * pa_xx * fz * fga * pb_zzzz

                     + 9.0 * fx * fx * fz * pb_zzzz

                     + 42.0 * pa_xx * fz * fx * pb_zzzz

                     + 16.0 * pa_xxxx * fz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xxxx_s_0(double fx,
                                     double pa_xxxy,
                                     double pa_xxy,
                                     double pa_xy,
                                     double pa_y,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxx,
                                     double s_0_0)
    {
        return s_0_0 * (5.625 * pa_xy * fx * fx * fx

                     + 7.5 * fx * fx * fx * pa_y * pb_x

                     + 0.75 * pa_xxxy * fx * fx

                     + 9.0 * pa_xxy * fx * fx * pb_x

                     + 13.5 * pa_xy * fx * fx * pb_xx

                     + 3.0 * fx * fx * pa_y * pb_xxx

                     + 3.0 * pa_xxxy * pb_xx * fx

                     + 6.0 * pa_xxy * fx * pb_xxx

                     + 1.5 * pa_xy * fx * pb_xxxx

                     + pa_xxxy * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xxxx_r_0(double fga,
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
                                     double pb_xxxx,
                                     double r_0_0)
    {
        return r_0_0 * (- 13.5 * pa_xy * fx * fx * fz * fgb

                     + 56.25 * pa_xy * fx * fx * fx * fz

                     + 75.0 * fx * fx * fx * fz * pa_y * pb_x

                     - 2.25 * pa_xy * fz * fga * fx * fx

                     - 9.0 * fx * fx * pa_y * pb_x * fz * fgb

                     - 9.0 * fx * fx * fz * fga * pa_y * pb_x

                     - 3.0 * pa_xxxy * fx * fz * fgb

                     - 18.0 * pa_xxy * fx * pb_x * fz * fgb

                     + 9.0 * pa_xxxy * fz * fx * fx

                     + 108.0 * pa_xxy * fx * fx * fz * pb_x

                     + 162.0 * pa_xy * fx * fx * fz * pb_xx

                     - 9.0 * pa_xy * fx * pb_xx * fz * fgb

                     - 9.0 * pa_xy * fz * fga * pb_xx * fx

                     - 6.0 * fx * fz * fga * pa_y * pb_xxx

                     - 6.0 * pa_xxxy * pb_xx * fz * fgb

                     + 36.0 * fx * fx * fz * pa_y * pb_xxx

                     + 42.0 * pa_xxxy * fz * pb_xx * fx

                     + 84.0 * pa_xxy * fx * fz * pb_xxx

                     - 3.0 * pa_xy * fz * fga * pb_xxxx

                     + 21.0 * pa_xy * fx * fz * pb_xxxx

                     + 16.0 * pa_xxxy * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xxxy_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxy,
                                     double pa_xxy,
                                     double pa_xy,
                                     double pa_y,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxy,
                                     double pb_xxy,
                                     double pb_xy,
                                     double pb_y,
                                     double s_0_0)
    {
        return s_0_0 * (0.9375 * fx * fx * fx * fx

                     + 1.125 * pa_xx * fx * fx * fx

                     + 3.375 * pa_x * fx * fx * fx * pb_x

                     + 1.875 * fx * fx * fx * pa_y * pb_y

                     + 1.125 * fx * fx * fx * pb_xx

                     + 0.75 * pa_xxx * fx * fx * pb_x

                     + 2.25 * pa_xxy * fx * fx * pb_y

                     + 2.25 * pa_xx * fx * fx * pb_xx

                     + 6.75 * pa_xy * fx * fx * pb_xy

                     + 0.75 * pa_x * fx * fx * pb_xxx

                     + 2.25 * fx * fx * pa_y * pb_xxy

                     + 1.5 * pa_xxxy * pb_xy * fx

                     + 0.5 * pa_xxx * fx * pb_xxx

                     + 4.5 * pa_xxy * fx * pb_xxy

                     + 1.5 * pa_xy * fx * pb_xxxy

                     + pa_xxxy * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xxxy_r_0(double fga,
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
                                     double pb_xxx,
                                     double pb_xxxy,
                                     double pb_xxy,
                                     double pb_xy,
                                     double pb_y,
                                     double r_0_0)
    {
        return r_0_0 * (7.5 * fx * fx * fx * fx * fz

                     - 1.125 * fx * fx * fx * fz * fgb

                     - 1.125 * fx * fx * fx * fz * fga

                     - 2.25 * pa_xx * fx * fx * fz * fgb

                     + 11.25 * pa_xx * fx * fx * fx * fz

                     + 33.75 * pa_x * fx * fx * fx * fz * pb_x

                     + 18.75 * fx * fx * fx * fz * pa_y * pb_y

                     - 2.25 * pa_x * fx * fx * pb_x * fz * fgb

                     - 2.25 * pa_x * fz * fga * fx * fx * pb_x

                     - 2.25 * fx * fx * pa_y * fz * fgb * pb_y

                     - 2.25 * fx * fx * fz * fga * pa_y * pb_y

                     - 2.25 * fx * fx * fz * fga * pb_xx

                     - 1.5 * pa_xxx * fx * pb_x * fz * fgb

                     - 4.5 * pa_xxy * fx * fz * fgb * pb_y

                     + 11.25 * fx * fx * fx * fz * pb_xx

                     + 9.0 * pa_xxx * fz * fx * fx * pb_x

                     + 27.0 * pa_xxy * fx * fx * fz * pb_y

                     + 27.0 * pa_xx * fx * fx * fz * pb_xx

                     + 81.0 * pa_xy * fx * fx * fz * pb_xy

                     - 4.5 * pa_xy * fx * pb_xy * fz * fgb

                     - 4.5 * pa_xy * fz * fga * pb_xy * fx

                     - 1.5 * pa_x * fz * fga * fx * pb_xxx

                     - 4.5 * fx * fz * fga * pa_y * pb_xxy

                     - 3.0 * pa_xxxy * pb_xy * fz * fgb

                     + 9.0 * pa_x * fx * fx * fz * pb_xxx

                     + 27.0 * fx * fx * fz * pa_y * pb_xxy

                     + 21.0 * pa_xxxy * fz * pb_xy * fx

                     + 7.0 * pa_xxx * fz * fx * pb_xxx

                     + 63.0 * pa_xxy * fx * fz * pb_xxy

                     - 3.0 * pa_xy * fz * fga * pb_xxxy

                     + 21.0 * pa_xy * fx * fz * pb_xxxy

                     + 16.0 * pa_xxxy * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xxxz_s_0(double fx,
                                     double pa_xxxy,
                                     double pa_xxy,
                                     double pa_xy,
                                     double pa_y,
                                     double pb_xxxz,
                                     double pb_xxz,
                                     double pb_xz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pa_y * pb_z

                     + 2.25 * pa_xxy * fx * fx * pb_z

                     + 6.75 * pa_xy * fx * fx * pb_xz

                     + 2.25 * fx * fx * pa_y * pb_xxz

                     + 1.5 * pa_xxxy * pb_xz * fx

                     + 4.5 * pa_xxy * fx * pb_xxz

                     + 1.5 * pa_xy * fx * pb_xxxz

                     + pa_xxxy * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xxxz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xxxy,
                                     double pa_xxy,
                                     double pa_xy,
                                     double pa_y,
                                     double pb_xxxz,
                                     double pb_xxz,
                                     double pb_xz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (18.75 * fx * fx * fx * fz * pa_y * pb_z

                     - 2.25 * fx * fx * pa_y * fz * fgb * pb_z

                     - 2.25 * fx * fx * fz * fga * pa_y * pb_z

                     - 4.5 * pa_xxy * fx * fz * fgb * pb_z

                     + 27.0 * pa_xxy * fx * fx * fz * pb_z

                     + 81.0 * pa_xy * fx * fx * fz * pb_xz

                     - 4.5 * pa_xy * fx * pb_xz * fz * fgb

                     - 4.5 * pa_xy * fz * fga * pb_xz * fx

                     - 4.5 * fx * fz * fga * pa_y * pb_xxz

                     - 3.0 * pa_xxxy * pb_xz * fz * fgb

                     + 27.0 * fx * fx * fz * pa_y * pb_xxz

                     + 21.0 * pa_xxxy * fz * pb_xz * fx

                     + 63.0 * pa_xxy * fx * fz * pb_xxz

                     - 3.0 * pa_xy * fz * fga * pb_xxxz

                     + 21.0 * pa_xy * fx * fz * pb_xxxz

                     + 16.0 * pa_xxxy * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xxyy_s_0(double fx,
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
                                     double pb_xxyy,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_y,
                                     double pb_yy,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_xy * fx * fx * fx

                     + 2.25 * pa_x * fx * fx * fx * pb_y

                     + 0.75 * fx * fx * fx * pa_y * pb_x

                     + 1.5 * fx * fx * fx * pb_xy

                     + 0.25 * pa_xxxy * fx * fx

                     + 0.5 * pa_xxx * fx * fx * pb_y

                     + 1.5 * pa_xxy * fx * fx * pb_x

                     + 3.0 * pa_xx * fx * fx * pb_xy

                     + 2.25 * pa_xy * fx * fx * pb_yy

                     + 0.75 * pa_xy * fx * fx * pb_xx

                     + 1.5 * pa_x * fx * fx * pb_xxy

                     + 1.5 * fx * fx * pa_y * pb_xyy

                     + 0.5 * pa_xxxy * pb_xx * fx

                     + 0.5 * pa_xxxy * fx * pb_yy

                     + pa_xxx * fx * pb_xxy

                     + 3.0 * pa_xxy * fx * pb_xyy

                     + 1.5 * pa_xy * fx * pb_xxyy

                     + pa_xxxy * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xxyy_r_0(double fga,
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
                                     double pb_xxyy,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_y,
                                     double pb_yy,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xy * fx * fx * fz * fgb

                     + 11.25 * pa_xy * fx * fx * fx * fz

                     + 22.5 * pa_x * fx * fx * fx * fz * pb_y

                     - 1.5 * pa_x * fx * fx * fz * fgb * pb_y

                     - 0.75 * pa_xy * fz * fga * fx * fx

                     - 1.5 * pa_x * fz * fga * fx * fx * pb_y

                     - 1.5 * fx * fx * pa_y * pb_x * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pa_y * pb_x

                     - 3.0 * fx * fx * fz * fga * pb_xy

                     - pa_xxxy * fx * fz * fgb

                     - pa_xxx * fx * fz * fgb * pb_y

                     - 3.0 * pa_xxy * fx * pb_x * fz * fgb

                     + 7.5 * fx * fx * fx * fz * pa_y * pb_x

                     + 15.0 * fx * fx * fx * fz * pb_xy

                     + 3.0 * pa_xxxy * fz * fx * fx

                     + 6.0 * pa_xxx * fz * fx * fx * pb_y

                     + 18.0 * pa_xxy * fx * fx * fz * pb_x

                     + 36.0 * pa_xx * fx * fx * fz * pb_xy

                     + 27.0 * pa_xy * fx * fx * fz * pb_yy

                     - 1.5 * pa_xy * fx * pb_xx * fz * fgb

                     - 1.5 * pa_xy * fx * fz * fgb * pb_yy

                     - 1.5 * pa_xy * fz * fga * pb_xx * fx

                     - 1.5 * pa_xy * fz * fga * fx * pb_yy

                     - 3.0 * pa_x * fz * fga * fx * pb_xxy

                     - 3.0 * fx * fz * fga * pa_y * pb_xyy

                     - pa_xxxy * pb_xx * fz * fgb

                     - pa_xxxy * fz * fgb * pb_yy

                     + 9.0 * pa_xy * fx * fx * fz * pb_xx

                     + 18.0 * pa_x * fx * fx * fz * pb_xxy

                     + 18.0 * fx * fx * fz * pa_y * pb_xyy

                     + 7.0 * pa_xxxy * fz * pb_xx * fx

                     + 7.0 * pa_xxxy * fz * fx * pb_yy

                     + 14.0 * pa_xxx * fz * fx * pb_xxy

                     + 42.0 * pa_xxy * fx * fz * pb_xyy

                     - 3.0 * pa_xy * fz * fga * pb_xxyy

                     + 21.0 * pa_xy * fx * fz * pb_xxyy

                     + 16.0 * pa_xxxy * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xxyz_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxy,
                                     double pa_xxy,
                                     double pa_xy,
                                     double pa_y,
                                     double pb_xxyz,
                                     double pb_xxz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_x * fx * fx * fx * pb_z

                     + 0.75 * fx * fx * fx * pb_xz

                     + 0.25 * pa_xxx * fx * fx * pb_z

                     + 1.5 * pa_xx * fx * fx * pb_xz

                     + 2.25 * pa_xy * fx * fx * pb_yz

                     + 0.75 * pa_x * fx * fx * pb_xxz

                     + 1.5 * fx * fx * pa_y * pb_xyz

                     + 0.5 * pa_xxxy * fx * pb_yz

                     + 0.5 * pa_xxx * fx * pb_xxz

                     + 3.0 * pa_xxy * fx * pb_xyz

                     + 1.5 * pa_xy * fx * pb_xxyz

                     + pa_xxxy * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xxyz_r_0(double fga,
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
                                     double pb_xxyz,
                                     double pb_xxz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (11.25 * pa_x * fx * fx * fx * fz * pb_z

                     - 0.75 * pa_x * fx * fx * fz * fgb * pb_z

                     - 0.75 * pa_x * fz * fga * fx * fx * pb_z

                     - 1.5 * fx * fx * fz * fga * pb_xz

                     - 0.5 * pa_xxx * fx * fz * fgb * pb_z

                     + 7.5 * fx * fx * fx * fz * pb_xz

                     + 3.0 * pa_xxx * fz * fx * fx * pb_z

                     + 18.0 * pa_xx * fx * fx * fz * pb_xz

                     + 27.0 * pa_xy * fx * fx * fz * pb_yz

                     - 1.5 * pa_xy * fx * fz * fgb * pb_yz

                     - 1.5 * pa_xy * fz * fga * fx * pb_yz

                     - 1.5 * pa_x * fz * fga * fx * pb_xxz

                     - 3.0 * fx * fz * fga * pa_y * pb_xyz

                     - pa_xxxy * fz * fgb * pb_yz

                     + 9.0 * pa_x * fx * fx * fz * pb_xxz

                     + 18.0 * fx * fx * fz * pa_y * pb_xyz

                     + 7.0 * pa_xxxy * fz * fx * pb_yz

                     + 7.0 * pa_xxx * fz * fx * pb_xxz

                     + 42.0 * pa_xxy * fx * fz * pb_xyz

                     - 3.0 * pa_xy * fz * fga * pb_xxyz

                     + 21.0 * pa_xy * fx * fz * pb_xxyz

                     + 16.0 * pa_xxxy * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xxzz_s_0(double fx,
                                     double pa_xxxy,
                                     double pa_xxy,
                                     double pa_xy,
                                     double pa_y,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxzz,
                                     double pb_xzz,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_xy * fx * fx * fx

                     + 0.75 * fx * fx * fx * pa_y * pb_x

                     + 0.25 * pa_xxxy * fx * fx

                     + 1.5 * pa_xxy * fx * fx * pb_x

                     + 2.25 * pa_xy * fx * fx * pb_zz

                     + 0.75 * pa_xy * fx * fx * pb_xx

                     + 1.5 * fx * fx * pa_y * pb_xzz

                     + 0.5 * pa_xxxy * pb_xx * fx

                     + 0.5 * pa_xxxy * fx * pb_zz

                     + 3.0 * pa_xxy * fx * pb_xzz

                     + 1.5 * pa_xy * fx * pb_xxzz

                     + pa_xxxy * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xxzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xxxy,
                                     double pa_xxy,
                                     double pa_xy,
                                     double pa_y,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxzz,
                                     double pb_xzz,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xy * fx * fx * fz * fgb

                     + 11.25 * pa_xy * fx * fx * fx * fz

                     - 0.75 * pa_xy * fz * fga * fx * fx

                     - 1.5 * fx * fx * pa_y * pb_x * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pa_y * pb_x

                     - pa_xxxy * fx * fz * fgb

                     - 3.0 * pa_xxy * fx * pb_x * fz * fgb

                     + 7.5 * fx * fx * fx * fz * pa_y * pb_x

                     + 3.0 * pa_xxxy * fz * fx * fx

                     + 18.0 * pa_xxy * fx * fx * fz * pb_x

                     + 27.0 * pa_xy * fx * fx * fz * pb_zz

                     - 1.5 * pa_xy * fx * pb_xx * fz * fgb

                     - 1.5 * pa_xy * fx * fz * fgb * pb_zz

                     - 1.5 * pa_xy * fz * fga * pb_xx * fx

                     - 1.5 * pa_xy * fz * fga * fx * pb_zz

                     - 3.0 * fx * fz * fga * pa_y * pb_xzz

                     - pa_xxxy * pb_xx * fz * fgb

                     - pa_xxxy * fz * fgb * pb_zz

                     + 9.0 * pa_xy * fx * fx * fz * pb_xx

                     + 18.0 * fx * fx * fz * pa_y * pb_xzz

                     + 7.0 * pa_xxxy * fz * pb_xx * fx

                     + 7.0 * pa_xxxy * fz * fx * pb_zz

                     + 42.0 * pa_xxy * fx * fz * pb_xzz

                     - 3.0 * pa_xy * fz * fga * pb_xxzz

                     + 21.0 * pa_xy * fx * fz * pb_xxzz

                     + 16.0 * pa_xxxy * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xyyy_s_0(double fx,
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
                                     double pb_xyyy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double s_0_0)
    {
        return s_0_0 * (0.5625 * fx * fx * fx * fx

                     + 1.125 * pa_xx * fx * fx * fx

                     + 1.125 * pa_x * fx * fx * fx * pb_x

                     + 1.125 * fx * fx * fx * pa_y * pb_y

                     + 1.125 * fx * fx * fx * pb_yy

                     + 0.75 * pa_xxx * fx * fx * pb_x

                     + 2.25 * pa_xxy * fx * fx * pb_y

                     + 2.25 * pa_xx * fx * fx * pb_yy

                     + 2.25 * pa_xy * fx * fx * pb_xy

                     + 2.25 * pa_x * fx * fx * pb_xyy

                     + 0.75 * fx * fx * pa_y * pb_yyy

                     + 1.5 * pa_xxxy * pb_xy * fx

                     + 1.5 * pa_xxx * fx * pb_xyy

                     + 1.5 * pa_xxy * fx * pb_yyy

                     + 1.5 * pa_xy * fx * pb_xyyy

                     + pa_xxxy * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xyyy_r_0(double fga,
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
                                     double pb_xyyy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.125 * fx * fx * fx * fz * fgb

                     - 1.125 * fx * fx * fx * fz * fga

                     - 2.25 * pa_xx * fx * fx * fz * fgb

                     + 4.5 * fx * fx * fx * fx * fz

                     + 11.25 * pa_xx * fx * fx * fx * fz

                     - 2.25 * pa_x * fx * fx * pb_x * fz * fgb

                     - 2.25 * pa_x * fz * fga * fx * fx * pb_x

                     - 2.25 * fx * fx * pa_y * pb_y * fz * fgb

                     - 2.25 * fx * fx * fz * fga * pa_y * pb_y

                     - 2.25 * fx * fx * fz * fga * pb_yy

                     - 1.5 * pa_xxx * fx * pb_x * fz * fgb

                     - 4.5 * pa_xxy * fx * pb_y * fz * fgb

                     + 11.25 * pa_x * fx * fx * fx * fz * pb_x

                     + 11.25 * fx * fx * fx * fz * pa_y * pb_y

                     + 11.25 * fx * fx * fx * fz * pb_yy

                     + 9.0 * pa_xxx * fz * fx * fx * pb_x

                     + 27.0 * pa_xxy * fx * fx * fz * pb_y

                     + 27.0 * pa_xx * fx * fx * fz * pb_yy

                     - 4.5 * pa_xy * fx * pb_xy * fz * fgb

                     - 4.5 * pa_xy * fz * fga * pb_xy * fx

                     - 4.5 * pa_x * fz * fga * fx * pb_xyy

                     - 1.5 * fx * fz * fga * pa_y * pb_yyy

                     - 3.0 * pa_xxxy * pb_xy * fz * fgb

                     + 27.0 * pa_xy * fx * fx * fz * pb_xy

                     + 27.0 * pa_x * fx * fx * fz * pb_xyy

                     + 9.0 * fx * fx * fz * pa_y * pb_yyy

                     + 21.0 * pa_xxxy * fz * pb_xy * fx

                     + 21.0 * pa_xxx * fz * fx * pb_xyy

                     + 21.0 * pa_xxy * fx * fz * pb_yyy

                     - 3.0 * pa_xy * fz * fga * pb_xyyy

                     + 21.0 * pa_xy * fx * fz * pb_xyyy

                     + 16.0 * pa_xxxy * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xyyz_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxy,
                                     double pa_xxy,
                                     double pa_xy,
                                     double pa_y,
                                     double pb_xyyz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_y * pb_z

                     + 0.75 * fx * fx * fx * pb_yz

                     + 0.75 * pa_xxy * fx * fx * pb_z

                     + 1.5 * pa_xx * fx * fx * pb_yz

                     + 0.75 * pa_xy * fx * fx * pb_xz

                     + 1.5 * pa_x * fx * fx * pb_xyz

                     + 0.75 * fx * fx * pa_y * pb_yyz

                     + 0.5 * pa_xxxy * pb_xz * fx

                     + pa_xxx * fx * pb_xyz

                     + 1.5 * pa_xxy * fx * pb_yyz

                     + 1.5 * pa_xy * fx * pb_xyyz

                     + pa_xxxy * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xyyz_r_0(double fga,
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
                                     double pb_xyyz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pa_y * fz * fgb * pb_z

                     - 0.75 * fx * fx * fz * fga * pa_y * pb_z

                     - 1.5 * fx * fx * fz * fga * pb_yz

                     - 1.5 * pa_xxy * fx * fz * fgb * pb_z

                     + 3.75 * fx * fx * fx * fz * pa_y * pb_z

                     + 7.5 * fx * fx * fx * fz * pb_yz

                     + 9.0 * pa_xxy * fx * fx * fz * pb_z

                     + 18.0 * pa_xx * fx * fx * fz * pb_yz

                     - 1.5 * pa_xy * fx * pb_xz * fz * fgb

                     - 1.5 * pa_xy * fz * fga * pb_xz * fx

                     - 3.0 * pa_x * fz * fga * fx * pb_xyz

                     - 1.5 * fx * fz * fga * pa_y * pb_yyz

                     - pa_xxxy * pb_xz * fz * fgb

                     + 9.0 * pa_xy * fx * fx * fz * pb_xz

                     + 18.0 * pa_x * fx * fx * fz * pb_xyz

                     + 9.0 * fx * fx * fz * pa_y * pb_yyz

                     + 7.0 * pa_xxxy * fz * pb_xz * fx

                     + 14.0 * pa_xxx * fz * fx * pb_xyz

                     + 21.0 * pa_xxy * fx * fz * pb_yyz

                     - 3.0 * pa_xy * fz * fga * pb_xyyz

                     + 21.0 * pa_xy * fx * fz * pb_xyyz

                     + 16.0 * pa_xxxy * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xyzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxy,
                                     double pa_xxy,
                                     double pa_xy,
                                     double pa_y,
                                     double pb_x,
                                     double pb_xy,
                                     double pb_xyzz,
                                     double pb_xzz,
                                     double pb_y,
                                     double pb_yzz,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * pa_xx * fx * fx * fx

                     + 0.375 * pa_x * fx * fx * fx * pb_x

                     + 0.375 * fx * fx * fx * pa_y * pb_y

                     + 0.375 * fx * fx * fx * pb_zz

                     + 0.25 * pa_xxx * fx * fx * pb_x

                     + 0.75 * pa_xxy * fx * fx * pb_y

                     + 0.75 * pa_xx * fx * fx * pb_zz

                     + 0.75 * pa_xy * fx * fx * pb_xy

                     + 0.75 * pa_x * fx * fx * pb_xzz

                     + 0.75 * fx * fx * pa_y * pb_yzz

                     + 0.5 * pa_xxxy * pb_xy * fx

                     + 0.5 * pa_xxx * fx * pb_xzz

                     + 1.5 * pa_xxy * fx * pb_yzz

                     + 1.5 * pa_xy * fx * pb_xyzz

                     + pa_xxxy * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xyzz_r_0(double fga,
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
                                     double pb_xyzz,
                                     double pb_xzz,
                                     double pb_y,
                                     double pb_yzz,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.375 * fx * fx * fx * fz * fgb

                     - 0.375 * fx * fx * fx * fz * fga

                     - 0.75 * pa_xx * fx * fx * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     + 3.75 * pa_xx * fx * fx * fx * fz

                     - 0.75 * pa_x * fx * fx * pb_x * fz * fgb

                     - 0.75 * pa_x * fz * fga * fx * fx * pb_x

                     - 0.75 * fx * fx * pa_y * pb_y * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pa_y * pb_y

                     - 0.75 * fx * fx * fz * fga * pb_zz

                     - 0.5 * pa_xxx * fx * pb_x * fz * fgb

                     - 1.5 * pa_xxy * fx * pb_y * fz * fgb

                     + 3.75 * pa_x * fx * fx * fx * fz * pb_x

                     + 3.75 * fx * fx * fx * fz * pa_y * pb_y

                     + 3.75 * fx * fx * fx * fz * pb_zz

                     + 3.0 * pa_xxx * fz * fx * fx * pb_x

                     + 9.0 * pa_xxy * fx * fx * fz * pb_y

                     + 9.0 * pa_xx * fx * fx * fz * pb_zz

                     - 1.5 * pa_xy * fx * pb_xy * fz * fgb

                     - 1.5 * pa_xy * fz * fga * pb_xy * fx

                     - 1.5 * pa_x * fz * fga * fx * pb_xzz

                     - 1.5 * fx * fz * fga * pa_y * pb_yzz

                     - pa_xxxy * pb_xy * fz * fgb

                     + 9.0 * pa_xy * fx * fx * fz * pb_xy

                     + 9.0 * pa_x * fx * fx * fz * pb_xzz

                     + 9.0 * fx * fx * fz * pa_y * pb_yzz

                     + 7.0 * pa_xxxy * fz * pb_xy * fx

                     + 7.0 * pa_xxx * fz * fx * pb_xzz

                     + 21.0 * pa_xxy * fx * fz * pb_yzz

                     - 3.0 * pa_xy * fz * fga * pb_xyzz

                     + 21.0 * pa_xy * fx * fz * pb_xyzz

                     + 16.0 * pa_xxxy * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xzzz_s_0(double fx,
                                     double pa_xxxy,
                                     double pa_xxy,
                                     double pa_xy,
                                     double pa_y,
                                     double pb_xz,
                                     double pb_xzzz,
                                     double pb_z,
                                     double pb_zzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pa_y * pb_z

                     + 2.25 * pa_xxy * fx * fx * pb_z

                     + 2.25 * pa_xy * fx * fx * pb_xz

                     + 0.75 * fx * fx * pa_y * pb_zzz

                     + 1.5 * pa_xxxy * pb_xz * fx

                     + 1.5 * pa_xxy * fx * pb_zzz

                     + 1.5 * pa_xy * fx * pb_xzzz

                     + pa_xxxy * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xzzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xxxy,
                                     double pa_xxy,
                                     double pa_xy,
                                     double pa_y,
                                     double pb_xz,
                                     double pb_xzzz,
                                     double pb_z,
                                     double pb_zzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pa_y * pb_z * fz * fgb

                     - 2.25 * fx * fx * fz * fga * pa_y * pb_z

                     - 4.5 * pa_xxy * fx * pb_z * fz * fgb

                     + 11.25 * fx * fx * fx * fz * pa_y * pb_z

                     + 27.0 * pa_xxy * fx * fx * fz * pb_z

                     - 4.5 * pa_xy * fx * pb_xz * fz * fgb

                     - 4.5 * pa_xy * fz * fga * pb_xz * fx

                     - 1.5 * fx * fz * fga * pa_y * pb_zzz

                     - 3.0 * pa_xxxy * pb_xz * fz * fgb

                     + 27.0 * pa_xy * fx * fx * fz * pb_xz

                     + 9.0 * fx * fx * fz * pa_y * pb_zzz

                     + 21.0 * pa_xxxy * fz * pb_xz * fx

                     + 21.0 * pa_xxy * fx * fz * pb_zzz

                     - 3.0 * pa_xy * fz * fga * pb_xzzz

                     + 21.0 * pa_xy * fx * fz * pb_xzzz

                     + 16.0 * pa_xxxy * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_yyyy_s_0(double fx,
                                     double pa_x,
                                     double pa_xxx,
                                     double pa_xxxy,
                                     double pa_xy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double pb_yyyy,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_xy * fx * fx * fx

                     + 4.5 * pa_x * fx * fx * fx * pb_y

                     + 0.75 * pa_xxxy * fx * fx

                     + 3.0 * pa_xxx * fx * fx * pb_y

                     + 4.5 * pa_xy * fx * fx * pb_yy

                     + 3.0 * pa_x * fx * fx * pb_yyy

                     + 3.0 * pa_xxxy * pb_yy * fx

                     + 2.0 * pa_xxx * fx * pb_yyy

                     + 1.5 * pa_xy * fx * pb_yyyy

                     + pa_xxxy * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_yyyy_r_0(double fga,
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
                                     double pb_yyyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_xy * fx * fx * fz * fgb

                     - 9.0 * pa_x * fx * fx * pb_y * fz * fgb

                     - 2.25 * pa_xy * fz * fga * fx * fx

                     - 9.0 * pa_x * fz * fga * fx * fx * pb_y

                     - 3.0 * pa_xxxy * fx * fz * fgb

                     - 6.0 * pa_xxx * fx * pb_y * fz * fgb

                     + 11.25 * pa_xy * fx * fx * fx * fz

                     + 45.0 * pa_x * fx * fx * fx * fz * pb_y

                     + 9.0 * pa_xxxy * fz * fx * fx

                     + 36.0 * pa_xxx * fz * fx * fx * pb_y

                     - 9.0 * pa_xy * fx * pb_yy * fz * fgb

                     - 9.0 * pa_xy * fz * fga * pb_yy * fx

                     - 6.0 * pa_x * fz * fga * fx * pb_yyy

                     - 6.0 * pa_xxxy * pb_yy * fz * fgb

                     + 54.0 * pa_xy * fx * fx * fz * pb_yy

                     + 36.0 * pa_x * fx * fx * fz * pb_yyy

                     + 42.0 * pa_xxxy * fz * pb_yy * fx

                     + 28.0 * pa_xxx * fz * fx * pb_yyy

                     - 3.0 * pa_xy * fz * fga * pb_yyyy

                     + 21.0 * pa_xy * fx * fz * pb_yyyy

                     + 16.0 * pa_xxxy * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_yyyz_s_0(double fx,
                                     double pa_x,
                                     double pa_xxx,
                                     double pa_xxxy,
                                     double pa_xy,
                                     double pb_yyyz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_x * fx * fx * fx * pb_z

                     + 0.75 * pa_xxx * fx * fx * pb_z

                     + 2.25 * pa_xy * fx * fx * pb_yz

                     + 2.25 * pa_x * fx * fx * pb_yyz

                     + 1.5 * pa_xxxy * pb_yz * fx

                     + 1.5 * pa_xxx * fx * pb_yyz

                     + 1.5 * pa_xy * fx * pb_yyyz

                     + pa_xxxy * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_yyyz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xxx,
                                     double pa_xxxy,
                                     double pa_xy,
                                     double pb_yyyz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * pa_x * fx * fx * fz * fgb * pb_z

                     - 2.25 * pa_x * fz * fga * fx * fx * pb_z

                     - 1.5 * pa_xxx * fx * fz * fgb * pb_z

                     + 11.25 * pa_x * fx * fx * fx * fz * pb_z

                     + 9.0 * pa_xxx * fz * fx * fx * pb_z

                     - 4.5 * pa_xy * fx * pb_yz * fz * fgb

                     - 4.5 * pa_xy * fz * fga * pb_yz * fx

                     - 4.5 * pa_x * fz * fga * fx * pb_yyz

                     - 3.0 * pa_xxxy * pb_yz * fz * fgb

                     + 27.0 * pa_xy * fx * fx * fz * pb_yz

                     + 27.0 * pa_x * fx * fx * fz * pb_yyz

                     + 21.0 * pa_xxxy * fz * pb_yz * fx

                     + 21.0 * pa_xxx * fz * fx * pb_yyz

                     - 3.0 * pa_xy * fz * fga * pb_yyyz

                     + 21.0 * pa_xy * fx * fz * pb_yyyz

                     + 16.0 * pa_xxxy * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_yyzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xxx,
                                     double pa_xxxy,
                                     double pa_xy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyzz,
                                     double pb_yzz,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * pa_xy * fx * fx * fx

                     + 0.75 * pa_x * fx * fx * fx * pb_y

                     + 0.25 * pa_xxxy * fx * fx

                     + 0.5 * pa_xxx * fx * fx * pb_y

                     + 0.75 * pa_xy * fx * fx * pb_yy

                     + 0.75 * pa_xy * fx * fx * pb_zz

                     + 1.5 * pa_x * fx * fx * pb_yzz

                     + 0.5 * pa_xxxy * pb_yy * fx

                     + 0.5 * pa_xxxy * fx * pb_zz

                     + pa_xxx * fx * pb_yzz

                     + 1.5 * pa_xy * fx * pb_yyzz

                     + pa_xxxy * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_yyzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xxx,
                                     double pa_xxxy,
                                     double pa_xy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyzz,
                                     double pb_yzz,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xy * fx * fx * fz * fgb

                     - 1.5 * pa_x * fx * fx * pb_y * fz * fgb

                     - 0.75 * pa_xy * fz * fga * fx * fx

                     - 1.5 * pa_x * fz * fga * fx * fx * pb_y

                     - pa_xxxy * fx * fz * fgb

                     - pa_xxx * fx * pb_y * fz * fgb

                     + 3.75 * pa_xy * fx * fx * fx * fz

                     + 7.5 * pa_x * fx * fx * fx * fz * pb_y

                     + 3.0 * pa_xxxy * fz * fx * fx

                     + 6.0 * pa_xxx * fz * fx * fx * pb_y

                     - 1.5 * pa_xy * fx * pb_yy * fz * fgb

                     - 1.5 * pa_xy * fx * fz * fgb * pb_zz

                     - 1.5 * pa_xy * fz * fga * pb_yy * fx

                     - 1.5 * pa_xy * fz * fga * fx * pb_zz

                     - 3.0 * pa_x * fz * fga * fx * pb_yzz

                     - pa_xxxy * pb_yy * fz * fgb

                     - pa_xxxy * fz * fgb * pb_zz

                     + 9.0 * pa_xy * fx * fx * fz * pb_yy

                     + 9.0 * pa_xy * fx * fx * fz * pb_zz

                     + 18.0 * pa_x * fx * fx * fz * pb_yzz

                     + 7.0 * pa_xxxy * fz * pb_yy * fx

                     + 7.0 * pa_xxxy * fz * fx * pb_zz

                     + 14.0 * pa_xxx * fz * fx * pb_yzz

                     - 3.0 * pa_xy * fz * fga * pb_yyzz

                     + 21.0 * pa_xy * fx * fz * pb_yyzz

                     + 16.0 * pa_xxxy * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_yzzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xxx,
                                     double pa_xxxy,
                                     double pa_xy,
                                     double pb_yz,
                                     double pb_yzzz,
                                     double pb_z,
                                     double pb_zzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_x * fx * fx * fx * pb_z

                     + 0.75 * pa_xxx * fx * fx * pb_z

                     + 2.25 * pa_xy * fx * fx * pb_yz

                     + 0.75 * pa_x * fx * fx * pb_zzz

                     + 1.5 * pa_xxxy * pb_yz * fx

                     + 0.5 * pa_xxx * fx * pb_zzz

                     + 1.5 * pa_xy * fx * pb_yzzz

                     + pa_xxxy * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_yzzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xxx,
                                     double pa_xxxy,
                                     double pa_xy,
                                     double pb_yz,
                                     double pb_yzzz,
                                     double pb_z,
                                     double pb_zzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * pa_x * fx * fx * pb_z * fz * fgb

                     - 2.25 * pa_x * fz * fga * fx * fx * pb_z

                     - 1.5 * pa_xxx * fx * pb_z * fz * fgb

                     + 11.25 * pa_x * fx * fx * fx * fz * pb_z

                     + 9.0 * pa_xxx * fz * fx * fx * pb_z

                     - 4.5 * pa_xy * fx * pb_yz * fz * fgb

                     - 4.5 * pa_xy * fz * fga * pb_yz * fx

                     - 1.5 * pa_x * fz * fga * fx * pb_zzz

                     - 3.0 * pa_xxxy * pb_yz * fz * fgb

                     + 27.0 * pa_xy * fx * fx * fz * pb_yz

                     + 9.0 * pa_x * fx * fx * fz * pb_zzz

                     + 21.0 * pa_xxxy * fz * pb_yz * fx

                     + 7.0 * pa_xxx * fz * fx * pb_zzz

                     - 3.0 * pa_xy * fz * fga * pb_yzzz

                     + 21.0 * pa_xy * fx * fz * pb_yzzz

                     + 16.0 * pa_xxxy * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_zzzz_s_0(double fx,
                                     double pa_xxxy,
                                     double pa_xy,
                                     double pb_zz,
                                     double pb_zzzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_xy * fx * fx * fx

                     + 0.75 * pa_xxxy * fx * fx

                     + 4.5 * pa_xy * fx * fx * pb_zz

                     + 3.0 * pa_xxxy * pb_zz * fx

                     + 1.5 * pa_xy * fx * pb_zzzz

                     + pa_xxxy * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_zzzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xxxy,
                                     double pa_xy,
                                     double pb_zz,
                                     double pb_zzzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_xy * fx * fx * fz * fgb

                     - 2.25 * pa_xy * fz * fga * fx * fx

                     - 3.0 * pa_xxxy * fx * fz * fgb

                     + 11.25 * pa_xy * fx * fx * fx * fz

                     + 9.0 * pa_xxxy * fz * fx * fx

                     - 9.0 * pa_xy * fx * pb_zz * fz * fgb

                     - 9.0 * pa_xy * fz * fga * pb_zz * fx

                     - 6.0 * pa_xxxy * pb_zz * fz * fgb

                     + 54.0 * pa_xy * fx * fx * fz * pb_zz

                     + 42.0 * pa_xxxy * fz * pb_zz * fx

                     - 3.0 * pa_xy * fz * fga * pb_zzzz

                     + 21.0 * pa_xy * fx * fz * pb_zzzz

                     + 16.0 * pa_xxxy * fz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xxxx_s_0(double fx,
                                     double pa_xxxz,
                                     double pa_xxz,
                                     double pa_xz,
                                     double pa_z,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxx,
                                     double s_0_0)
    {
        return s_0_0 * (5.625 * pa_xz * fx * fx * fx

                     + 7.5 * fx * fx * fx * pa_z * pb_x

                     + 0.75 * pa_xxxz * fx * fx

                     + 9.0 * pa_xxz * fx * fx * pb_x

                     + 13.5 * pa_xz * fx * fx * pb_xx

                     + 3.0 * fx * fx * pa_z * pb_xxx

                     + 3.0 * pa_xxxz * pb_xx * fx

                     + 6.0 * pa_xxz * fx * pb_xxx

                     + 1.5 * pa_xz * fx * pb_xxxx

                     + pa_xxxz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xxxx_r_0(double fga,
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
                                     double pb_xxxx,
                                     double r_0_0)
    {
        return r_0_0 * (- 13.5 * pa_xz * fx * fx * fz * fgb

                     + 56.25 * pa_xz * fx * fx * fx * fz

                     + 75.0 * fx * fx * fx * fz * pa_z * pb_x

                     - 2.25 * pa_xz * fz * fga * fx * fx

                     - 9.0 * fx * fx * pa_z * pb_x * fz * fgb

                     - 9.0 * fx * fx * fz * fga * pa_z * pb_x

                     - 3.0 * pa_xxxz * fx * fz * fgb

                     - 18.0 * pa_xxz * fx * pb_x * fz * fgb

                     + 9.0 * pa_xxxz * fz * fx * fx

                     + 108.0 * pa_xxz * fx * fx * fz * pb_x

                     + 162.0 * pa_xz * fx * fx * fz * pb_xx

                     - 9.0 * pa_xz * fx * pb_xx * fz * fgb

                     - 9.0 * pa_xz * fz * fga * pb_xx * fx

                     - 6.0 * fx * fz * fga * pa_z * pb_xxx

                     - 6.0 * pa_xxxz * pb_xx * fz * fgb

                     + 36.0 * fx * fx * fz * pa_z * pb_xxx

                     + 42.0 * pa_xxxz * fz * pb_xx * fx

                     + 84.0 * pa_xxz * fx * fz * pb_xxx

                     - 3.0 * pa_xz * fz * fga * pb_xxxx

                     + 21.0 * pa_xz * fx * fz * pb_xxxx

                     + 16.0 * pa_xxxz * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xxxy_s_0(double fx,
                                     double pa_xxxz,
                                     double pa_xxz,
                                     double pa_xz,
                                     double pa_z,
                                     double pb_xxxy,
                                     double pb_xxy,
                                     double pb_xy,
                                     double pb_y,
                                     double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pa_z * pb_y

                     + 2.25 * pa_xxz * fx * fx * pb_y

                     + 6.75 * pa_xz * fx * fx * pb_xy

                     + 2.25 * fx * fx * pa_z * pb_xxy

                     + 1.5 * pa_xxxz * pb_xy * fx

                     + 4.5 * pa_xxz * fx * pb_xxy

                     + 1.5 * pa_xz * fx * pb_xxxy

                     + pa_xxxz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xxxy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xxxz,
                                     double pa_xxz,
                                     double pa_xz,
                                     double pa_z,
                                     double pb_xxxy,
                                     double pb_xxy,
                                     double pb_xy,
                                     double pb_y,
                                     double r_0_0)
    {
        return r_0_0 * (18.75 * fx * fx * fx * fz * pa_z * pb_y

                     - 2.25 * fx * fx * pa_z * fz * fgb * pb_y

                     - 2.25 * fx * fx * fz * fga * pa_z * pb_y

                     - 4.5 * pa_xxz * fx * fz * fgb * pb_y

                     + 27.0 * pa_xxz * fx * fx * fz * pb_y

                     + 81.0 * pa_xz * fx * fx * fz * pb_xy

                     - 4.5 * pa_xz * fx * pb_xy * fz * fgb

                     - 4.5 * pa_xz * fz * fga * pb_xy * fx

                     - 4.5 * fx * fz * fga * pa_z * pb_xxy

                     - 3.0 * pa_xxxz * pb_xy * fz * fgb

                     + 27.0 * fx * fx * fz * pa_z * pb_xxy

                     + 21.0 * pa_xxxz * fz * pb_xy * fx

                     + 63.0 * pa_xxz * fx * fz * pb_xxy

                     - 3.0 * pa_xz * fz * fga * pb_xxxy

                     + 21.0 * pa_xz * fx * fz * pb_xxxy

                     + 16.0 * pa_xxxz * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xxxz_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxz,
                                     double pa_xxz,
                                     double pa_xz,
                                     double pa_z,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxz,
                                     double pb_xxz,
                                     double pb_xz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.9375 * fx * fx * fx * fx

                     + 1.125 * pa_xx * fx * fx * fx

                     + 3.375 * pa_x * fx * fx * fx * pb_x

                     + 1.875 * fx * fx * fx * pa_z * pb_z

                     + 1.125 * fx * fx * fx * pb_xx

                     + 0.75 * pa_xxx * fx * fx * pb_x

                     + 2.25 * pa_xxz * fx * fx * pb_z

                     + 2.25 * pa_xx * fx * fx * pb_xx

                     + 6.75 * pa_xz * fx * fx * pb_xz

                     + 0.75 * pa_x * fx * fx * pb_xxx

                     + 2.25 * fx * fx * pa_z * pb_xxz

                     + 1.5 * pa_xxxz * pb_xz * fx

                     + 0.5 * pa_xxx * fx * pb_xxx

                     + 4.5 * pa_xxz * fx * pb_xxz

                     + 1.5 * pa_xz * fx * pb_xxxz

                     + pa_xxxz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xxxz_r_0(double fga,
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
                                     double pb_xxx,
                                     double pb_xxxz,
                                     double pb_xxz,
                                     double pb_xz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (7.5 * fx * fx * fx * fx * fz

                     - 1.125 * fx * fx * fx * fz * fgb

                     - 1.125 * fx * fx * fx * fz * fga

                     - 2.25 * pa_xx * fx * fx * fz * fgb

                     + 11.25 * pa_xx * fx * fx * fx * fz

                     + 33.75 * pa_x * fx * fx * fx * fz * pb_x

                     + 18.75 * fx * fx * fx * fz * pa_z * pb_z

                     - 2.25 * pa_x * fx * fx * pb_x * fz * fgb

                     - 2.25 * pa_x * fz * fga * fx * fx * pb_x

                     - 2.25 * fx * fx * pa_z * fz * fgb * pb_z

                     - 2.25 * fx * fx * fz * fga * pa_z * pb_z

                     - 2.25 * fx * fx * fz * fga * pb_xx

                     - 1.5 * pa_xxx * fx * pb_x * fz * fgb

                     - 4.5 * pa_xxz * fx * fz * fgb * pb_z

                     + 11.25 * fx * fx * fx * fz * pb_xx

                     + 9.0 * pa_xxx * fz * fx * fx * pb_x

                     + 27.0 * pa_xxz * fx * fx * fz * pb_z

                     + 27.0 * pa_xx * fx * fx * fz * pb_xx

                     + 81.0 * pa_xz * fx * fx * fz * pb_xz

                     - 4.5 * pa_xz * fx * pb_xz * fz * fgb

                     - 4.5 * pa_xz * fz * fga * pb_xz * fx

                     - 1.5 * pa_x * fz * fga * fx * pb_xxx

                     - 4.5 * fx * fz * fga * pa_z * pb_xxz

                     - 3.0 * pa_xxxz * pb_xz * fz * fgb

                     + 9.0 * pa_x * fx * fx * fz * pb_xxx

                     + 27.0 * fx * fx * fz * pa_z * pb_xxz

                     + 21.0 * pa_xxxz * fz * pb_xz * fx

                     + 7.0 * pa_xxx * fz * fx * pb_xxx

                     + 63.0 * pa_xxz * fx * fz * pb_xxz

                     - 3.0 * pa_xz * fz * fga * pb_xxxz

                     + 21.0 * pa_xz * fx * fz * pb_xxxz

                     + 16.0 * pa_xxxz * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xxyy_s_0(double fx,
                                     double pa_xxxz,
                                     double pa_xxz,
                                     double pa_xz,
                                     double pa_z,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxyy,
                                     double pb_xyy,
                                     double pb_yy,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_xz * fx * fx * fx

                     + 0.75 * fx * fx * fx * pa_z * pb_x

                     + 0.25 * pa_xxxz * fx * fx

                     + 1.5 * pa_xxz * fx * fx * pb_x

                     + 2.25 * pa_xz * fx * fx * pb_yy

                     + 0.75 * pa_xz * fx * fx * pb_xx

                     + 1.5 * fx * fx * pa_z * pb_xyy

                     + 0.5 * pa_xxxz * pb_xx * fx

                     + 0.5 * pa_xxxz * fx * pb_yy

                     + 3.0 * pa_xxz * fx * pb_xyy

                     + 1.5 * pa_xz * fx * pb_xxyy

                     + pa_xxxz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xxyy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xxxz,
                                     double pa_xxz,
                                     double pa_xz,
                                     double pa_z,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxyy,
                                     double pb_xyy,
                                     double pb_yy,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xz * fx * fx * fz * fgb

                     + 11.25 * pa_xz * fx * fx * fx * fz

                     - 0.75 * pa_xz * fz * fga * fx * fx

                     - 1.5 * fx * fx * pa_z * pb_x * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pa_z * pb_x

                     - pa_xxxz * fx * fz * fgb

                     - 3.0 * pa_xxz * fx * pb_x * fz * fgb

                     + 7.5 * fx * fx * fx * fz * pa_z * pb_x

                     + 3.0 * pa_xxxz * fz * fx * fx

                     + 18.0 * pa_xxz * fx * fx * fz * pb_x

                     + 27.0 * pa_xz * fx * fx * fz * pb_yy

                     - 1.5 * pa_xz * fx * pb_xx * fz * fgb

                     - 1.5 * pa_xz * fx * fz * fgb * pb_yy

                     - 1.5 * pa_xz * fz * fga * pb_xx * fx

                     - 1.5 * pa_xz * fz * fga * fx * pb_yy

                     - 3.0 * fx * fz * fga * pa_z * pb_xyy

                     - pa_xxxz * pb_xx * fz * fgb

                     - pa_xxxz * fz * fgb * pb_yy

                     + 9.0 * pa_xz * fx * fx * fz * pb_xx

                     + 18.0 * fx * fx * fz * pa_z * pb_xyy

                     + 7.0 * pa_xxxz * fz * pb_xx * fx

                     + 7.0 * pa_xxxz * fz * fx * pb_yy

                     + 42.0 * pa_xxz * fx * fz * pb_xyy

                     - 3.0 * pa_xz * fz * fga * pb_xxyy

                     + 21.0 * pa_xz * fx * fz * pb_xxyy

                     + 16.0 * pa_xxxz * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xxyz_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxz,
                                     double pa_xxz,
                                     double pa_xz,
                                     double pa_z,
                                     double pb_xxy,
                                     double pb_xxyz,
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_y,
                                     double pb_yz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_x * fx * fx * fx * pb_y

                     + 0.75 * fx * fx * fx * pb_xy

                     + 0.25 * pa_xxx * fx * fx * pb_y

                     + 1.5 * pa_xx * fx * fx * pb_xy

                     + 2.25 * pa_xz * fx * fx * pb_yz

                     + 0.75 * pa_x * fx * fx * pb_xxy

                     + 1.5 * fx * fx * pa_z * pb_xyz

                     + 0.5 * pa_xxxz * fx * pb_yz

                     + 0.5 * pa_xxx * fx * pb_xxy

                     + 3.0 * pa_xxz * fx * pb_xyz

                     + 1.5 * pa_xz * fx * pb_xxyz

                     + pa_xxxz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xxyz_r_0(double fga,
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
                                     double pb_xxy,
                                     double pb_xxyz,
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_y,
                                     double pb_yz,
                                     double r_0_0)
    {
        return r_0_0 * (11.25 * pa_x * fx * fx * fx * fz * pb_y

                     - 0.75 * pa_x * fx * fx * fz * fgb * pb_y

                     - 0.75 * pa_x * fz * fga * fx * fx * pb_y

                     - 1.5 * fx * fx * fz * fga * pb_xy

                     - 0.5 * pa_xxx * fx * fz * fgb * pb_y

                     + 7.5 * fx * fx * fx * fz * pb_xy

                     + 3.0 * pa_xxx * fz * fx * fx * pb_y

                     + 18.0 * pa_xx * fx * fx * fz * pb_xy

                     + 27.0 * pa_xz * fx * fx * fz * pb_yz

                     - 1.5 * pa_xz * fx * fz * fgb * pb_yz

                     - 1.5 * pa_xz * fz * fga * fx * pb_yz

                     - 1.5 * pa_x * fz * fga * fx * pb_xxy

                     - 3.0 * fx * fz * fga * pa_z * pb_xyz

                     - pa_xxxz * fz * fgb * pb_yz

                     + 9.0 * pa_x * fx * fx * fz * pb_xxy

                     + 18.0 * fx * fx * fz * pa_z * pb_xyz

                     + 7.0 * pa_xxxz * fz * fx * pb_yz

                     + 7.0 * pa_xxx * fz * fx * pb_xxy

                     + 42.0 * pa_xxz * fx * fz * pb_xyz

                     - 3.0 * pa_xz * fz * fga * pb_xxyz

                     + 21.0 * pa_xz * fx * fz * pb_xxyz

                     + 16.0 * pa_xxxz * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xxzz_s_0(double fx,
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
                                     double pb_xxzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_xz * fx * fx * fx

                     + 2.25 * pa_x * fx * fx * fx * pb_z

                     + 0.75 * fx * fx * fx * pa_z * pb_x

                     + 1.5 * fx * fx * fx * pb_xz

                     + 0.25 * pa_xxxz * fx * fx

                     + 0.5 * pa_xxx * fx * fx * pb_z

                     + 1.5 * pa_xxz * fx * fx * pb_x

                     + 3.0 * pa_xx * fx * fx * pb_xz

                     + 2.25 * pa_xz * fx * fx * pb_zz

                     + 0.75 * pa_xz * fx * fx * pb_xx

                     + 1.5 * pa_x * fx * fx * pb_xxz

                     + 1.5 * fx * fx * pa_z * pb_xzz

                     + 0.5 * pa_xxxz * pb_xx * fx

                     + 0.5 * pa_xxxz * fx * pb_zz

                     + pa_xxx * fx * pb_xxz

                     + 3.0 * pa_xxz * fx * pb_xzz

                     + 1.5 * pa_xz * fx * pb_xxzz

                     + pa_xxxz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xxzz_r_0(double fga,
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
                                     double pb_xxzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xz * fx * fx * fz * fgb

                     + 11.25 * pa_xz * fx * fx * fx * fz

                     + 22.5 * pa_x * fx * fx * fx * fz * pb_z

                     - 1.5 * pa_x * fx * fx * fz * fgb * pb_z

                     - 0.75 * pa_xz * fz * fga * fx * fx

                     - 1.5 * pa_x * fz * fga * fx * fx * pb_z

                     - 1.5 * fx * fx * pa_z * pb_x * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pa_z * pb_x

                     - 3.0 * fx * fx * fz * fga * pb_xz

                     - pa_xxxz * fx * fz * fgb

                     - pa_xxx * fx * fz * fgb * pb_z

                     - 3.0 * pa_xxz * fx * pb_x * fz * fgb

                     + 7.5 * fx * fx * fx * fz * pa_z * pb_x

                     + 15.0 * fx * fx * fx * fz * pb_xz

                     + 3.0 * pa_xxxz * fz * fx * fx

                     + 6.0 * pa_xxx * fz * fx * fx * pb_z

                     + 18.0 * pa_xxz * fx * fx * fz * pb_x

                     + 36.0 * pa_xx * fx * fx * fz * pb_xz

                     + 27.0 * pa_xz * fx * fx * fz * pb_zz

                     - 1.5 * pa_xz * fx * pb_xx * fz * fgb

                     - 1.5 * pa_xz * fx * fz * fgb * pb_zz

                     - 1.5 * pa_xz * fz * fga * pb_xx * fx

                     - 1.5 * pa_xz * fz * fga * fx * pb_zz

                     - 3.0 * pa_x * fz * fga * fx * pb_xxz

                     - 3.0 * fx * fz * fga * pa_z * pb_xzz

                     - pa_xxxz * pb_xx * fz * fgb

                     - pa_xxxz * fz * fgb * pb_zz

                     + 9.0 * pa_xz * fx * fx * fz * pb_xx

                     + 18.0 * pa_x * fx * fx * fz * pb_xxz

                     + 18.0 * fx * fx * fz * pa_z * pb_xzz

                     + 7.0 * pa_xxxz * fz * pb_xx * fx

                     + 7.0 * pa_xxxz * fz * fx * pb_zz

                     + 14.0 * pa_xxx * fz * fx * pb_xxz

                     + 42.0 * pa_xxz * fx * fz * pb_xzz

                     - 3.0 * pa_xz * fz * fga * pb_xxzz

                     + 21.0 * pa_xz * fx * fz * pb_xxzz

                     + 16.0 * pa_xxxz * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xyyy_s_0(double fx,
                                     double pa_xxxz,
                                     double pa_xxz,
                                     double pa_xz,
                                     double pa_z,
                                     double pb_xy,
                                     double pb_xyyy,
                                     double pb_y,
                                     double pb_yyy,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pa_z * pb_y

                     + 2.25 * pa_xxz * fx * fx * pb_y

                     + 2.25 * pa_xz * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_z * pb_yyy

                     + 1.5 * pa_xxxz * pb_xy * fx

                     + 1.5 * pa_xxz * fx * pb_yyy

                     + 1.5 * pa_xz * fx * pb_xyyy

                     + pa_xxxz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xyyy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xxxz,
                                     double pa_xxz,
                                     double pa_xz,
                                     double pa_z,
                                     double pb_xy,
                                     double pb_xyyy,
                                     double pb_y,
                                     double pb_yyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pa_z * pb_y * fz * fgb

                     - 2.25 * fx * fx * fz * fga * pa_z * pb_y

                     - 4.5 * pa_xxz * fx * pb_y * fz * fgb

                     + 11.25 * fx * fx * fx * fz * pa_z * pb_y

                     + 27.0 * pa_xxz * fx * fx * fz * pb_y

                     - 4.5 * pa_xz * fx * pb_xy * fz * fgb

                     - 4.5 * pa_xz * fz * fga * pb_xy * fx

                     - 1.5 * fx * fz * fga * pa_z * pb_yyy

                     - 3.0 * pa_xxxz * pb_xy * fz * fgb

                     + 27.0 * pa_xz * fx * fx * fz * pb_xy

                     + 9.0 * fx * fx * fz * pa_z * pb_yyy

                     + 21.0 * pa_xxxz * fz * pb_xy * fx

                     + 21.0 * pa_xxz * fx * fz * pb_yyy

                     - 3.0 * pa_xz * fz * fga * pb_xyyy

                     + 21.0 * pa_xz * fx * fz * pb_xyyy

                     + 16.0 * pa_xxxz * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xyyz_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxz,
                                     double pa_xxz,
                                     double pa_xz,
                                     double pa_z,
                                     double pb_x,
                                     double pb_xyy,
                                     double pb_xyyz,
                                     double pb_xz,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * pa_xx * fx * fx * fx

                     + 0.375 * pa_x * fx * fx * fx * pb_x

                     + 0.375 * fx * fx * fx * pa_z * pb_z

                     + 0.375 * fx * fx * fx * pb_yy

                     + 0.25 * pa_xxx * fx * fx * pb_x

                     + 0.75 * pa_xxz * fx * fx * pb_z

                     + 0.75 * pa_xx * fx * fx * pb_yy

                     + 0.75 * pa_xz * fx * fx * pb_xz

                     + 0.75 * pa_x * fx * fx * pb_xyy

                     + 0.75 * fx * fx * pa_z * pb_yyz

                     + 0.5 * pa_xxxz * pb_xz * fx

                     + 0.5 * pa_xxx * fx * pb_xyy

                     + 1.5 * pa_xxz * fx * pb_yyz

                     + 1.5 * pa_xz * fx * pb_xyyz

                     + pa_xxxz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xyyz_r_0(double fga,
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
                                     double pb_xyy,
                                     double pb_xyyz,
                                     double pb_xz,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.375 * fx * fx * fx * fz * fgb

                     - 0.375 * fx * fx * fx * fz * fga

                     - 0.75 * pa_xx * fx * fx * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     + 3.75 * pa_xx * fx * fx * fx * fz

                     - 0.75 * pa_x * fx * fx * pb_x * fz * fgb

                     - 0.75 * pa_x * fz * fga * fx * fx * pb_x

                     - 0.75 * fx * fx * pa_z * fz * fgb * pb_z

                     - 0.75 * fx * fx * fz * fga * pa_z * pb_z

                     - 0.75 * fx * fx * fz * fga * pb_yy

                     - 0.5 * pa_xxx * fx * pb_x * fz * fgb

                     - 1.5 * pa_xxz * fx * fz * fgb * pb_z

                     + 3.75 * pa_x * fx * fx * fx * fz * pb_x

                     + 3.75 * fx * fx * fx * fz * pa_z * pb_z

                     + 3.75 * fx * fx * fx * fz * pb_yy

                     + 3.0 * pa_xxx * fz * fx * fx * pb_x

                     + 9.0 * pa_xxz * fx * fx * fz * pb_z

                     + 9.0 * pa_xx * fx * fx * fz * pb_yy

                     - 1.5 * pa_xz * fx * pb_xz * fz * fgb

                     - 1.5 * pa_xz * fz * fga * pb_xz * fx

                     - 1.5 * pa_x * fz * fga * fx * pb_xyy

                     - 1.5 * fx * fz * fga * pa_z * pb_yyz

                     - pa_xxxz * pb_xz * fz * fgb

                     + 9.0 * pa_xz * fx * fx * fz * pb_xz

                     + 9.0 * pa_x * fx * fx * fz * pb_xyy

                     + 9.0 * fx * fx * fz * pa_z * pb_yyz

                     + 7.0 * pa_xxxz * fz * pb_xz * fx

                     + 7.0 * pa_xxx * fz * fx * pb_xyy

                     + 21.0 * pa_xxz * fx * fz * pb_yyz

                     - 3.0 * pa_xz * fz * fga * pb_xyyz

                     + 21.0 * pa_xz * fx * fz * pb_xyyz

                     + 16.0 * pa_xxxz * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xyzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxx,
                                     double pa_xxxz,
                                     double pa_xxz,
                                     double pa_xz,
                                     double pa_z,
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_xyzz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_yzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_z * pb_y

                     + 0.75 * fx * fx * fx * pb_yz

                     + 0.75 * pa_xxz * fx * fx * pb_y

                     + 1.5 * pa_xx * fx * fx * pb_yz

                     + 0.75 * pa_xz * fx * fx * pb_xy

                     + 1.5 * pa_x * fx * fx * pb_xyz

                     + 0.75 * fx * fx * pa_z * pb_yzz

                     + 0.5 * pa_xxxz * pb_xy * fx

                     + pa_xxx * fx * pb_xyz

                     + 1.5 * pa_xxz * fx * pb_yzz

                     + 1.5 * pa_xz * fx * pb_xyzz

                     + pa_xxxz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xyzz_r_0(double fga,
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
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_xyzz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_yzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pa_z * pb_y * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pa_z * pb_y

                     - 1.5 * fx * fx * fz * fga * pb_yz

                     - 1.5 * pa_xxz * fx * pb_y * fz * fgb

                     + 3.75 * fx * fx * fx * fz * pa_z * pb_y

                     + 7.5 * fx * fx * fx * fz * pb_yz

                     + 9.0 * pa_xxz * fx * fx * fz * pb_y

                     + 18.0 * pa_xx * fx * fx * fz * pb_yz

                     - 1.5 * pa_xz * fx * pb_xy * fz * fgb

                     - 1.5 * pa_xz * fz * fga * pb_xy * fx

                     - 3.0 * pa_x * fz * fga * fx * pb_xyz

                     - 1.5 * fx * fz * fga * pa_z * pb_yzz

                     - pa_xxxz * pb_xy * fz * fgb

                     + 9.0 * pa_xz * fx * fx * fz * pb_xy

                     + 18.0 * pa_x * fx * fx * fz * pb_xyz

                     + 9.0 * fx * fx * fz * pa_z * pb_yzz

                     + 7.0 * pa_xxxz * fz * pb_xy * fx

                     + 14.0 * pa_xxx * fz * fx * pb_xyz

                     + 21.0 * pa_xxz * fx * fz * pb_yzz

                     - 3.0 * pa_xz * fz * fga * pb_xyzz

                     + 21.0 * pa_xz * fx * fz * pb_xyzz

                     + 16.0 * pa_xxxz * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xzzz_s_0(double fx,
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
                                     double pb_xzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.5625 * fx * fx * fx * fx

                     + 1.125 * pa_xx * fx * fx * fx

                     + 1.125 * pa_x * fx * fx * fx * pb_x

                     + 1.125 * fx * fx * fx * pa_z * pb_z

                     + 1.125 * fx * fx * fx * pb_zz

                     + 0.75 * pa_xxx * fx * fx * pb_x

                     + 2.25 * pa_xxz * fx * fx * pb_z

                     + 2.25 * pa_xx * fx * fx * pb_zz

                     + 2.25 * pa_xz * fx * fx * pb_xz

                     + 2.25 * pa_x * fx * fx * pb_xzz

                     + 0.75 * fx * fx * pa_z * pb_zzz

                     + 1.5 * pa_xxxz * pb_xz * fx

                     + 1.5 * pa_xxx * fx * pb_xzz

                     + 1.5 * pa_xxz * fx * pb_zzz

                     + 1.5 * pa_xz * fx * pb_xzzz

                     + pa_xxxz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xzzz_r_0(double fga,
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
                                     double pb_xzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.125 * fx * fx * fx * fz * fgb

                     - 1.125 * fx * fx * fx * fz * fga

                     - 2.25 * pa_xx * fx * fx * fz * fgb

                     + 4.5 * fx * fx * fx * fx * fz

                     + 11.25 * pa_xx * fx * fx * fx * fz

                     - 2.25 * pa_x * fx * fx * pb_x * fz * fgb

                     - 2.25 * pa_x * fz * fga * fx * fx * pb_x

                     - 2.25 * fx * fx * pa_z * pb_z * fz * fgb

                     - 2.25 * fx * fx * fz * fga * pa_z * pb_z

                     - 2.25 * fx * fx * fz * fga * pb_zz

                     - 1.5 * pa_xxx * fx * pb_x * fz * fgb

                     - 4.5 * pa_xxz * fx * pb_z * fz * fgb

                     + 11.25 * pa_x * fx * fx * fx * fz * pb_x

                     + 11.25 * fx * fx * fx * fz * pa_z * pb_z

                     + 11.25 * fx * fx * fx * fz * pb_zz

                     + 9.0 * pa_xxx * fz * fx * fx * pb_x

                     + 27.0 * pa_xxz * fx * fx * fz * pb_z

                     + 27.0 * pa_xx * fx * fx * fz * pb_zz

                     - 4.5 * pa_xz * fx * pb_xz * fz * fgb

                     - 4.5 * pa_xz * fz * fga * pb_xz * fx

                     - 4.5 * pa_x * fz * fga * fx * pb_xzz

                     - 1.5 * fx * fz * fga * pa_z * pb_zzz

                     - 3.0 * pa_xxxz * pb_xz * fz * fgb

                     + 27.0 * pa_xz * fx * fx * fz * pb_xz

                     + 27.0 * pa_x * fx * fx * fz * pb_xzz

                     + 9.0 * fx * fx * fz * pa_z * pb_zzz

                     + 21.0 * pa_xxxz * fz * pb_xz * fx

                     + 21.0 * pa_xxx * fz * fx * pb_xzz

                     + 21.0 * pa_xxz * fx * fz * pb_zzz

                     - 3.0 * pa_xz * fz * fga * pb_xzzz

                     + 21.0 * pa_xz * fx * fz * pb_xzzz

                     + 16.0 * pa_xxxz * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_yyyy_s_0(double fx,
                                     double pa_xxxz,
                                     double pa_xz,
                                     double pb_yy,
                                     double pb_yyyy,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_xz * fx * fx * fx

                     + 0.75 * pa_xxxz * fx * fx

                     + 4.5 * pa_xz * fx * fx * pb_yy

                     + 3.0 * pa_xxxz * pb_yy * fx

                     + 1.5 * pa_xz * fx * pb_yyyy

                     + pa_xxxz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_yyyy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xxxz,
                                     double pa_xz,
                                     double pb_yy,
                                     double pb_yyyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_xz * fx * fx * fz * fgb

                     - 2.25 * pa_xz * fz * fga * fx * fx

                     - 3.0 * pa_xxxz * fx * fz * fgb

                     + 11.25 * pa_xz * fx * fx * fx * fz

                     + 9.0 * pa_xxxz * fz * fx * fx

                     - 9.0 * pa_xz * fx * pb_yy * fz * fgb

                     - 9.0 * pa_xz * fz * fga * pb_yy * fx

                     - 6.0 * pa_xxxz * pb_yy * fz * fgb

                     + 54.0 * pa_xz * fx * fx * fz * pb_yy

                     + 42.0 * pa_xxxz * fz * pb_yy * fx

                     - 3.0 * pa_xz * fz * fga * pb_yyyy

                     + 21.0 * pa_xz * fx * fz * pb_yyyy

                     + 16.0 * pa_xxxz * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_yyyz_s_0(double fx,
                                     double pa_x,
                                     double pa_xxx,
                                     double pa_xxxz,
                                     double pa_xz,
                                     double pb_y,
                                     double pb_yyy,
                                     double pb_yyyz,
                                     double pb_yz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_x * fx * fx * fx * pb_y

                     + 0.75 * pa_xxx * fx * fx * pb_y

                     + 2.25 * pa_xz * fx * fx * pb_yz

                     + 0.75 * pa_x * fx * fx * pb_yyy

                     + 1.5 * pa_xxxz * pb_yz * fx

                     + 0.5 * pa_xxx * fx * pb_yyy

                     + 1.5 * pa_xz * fx * pb_yyyz

                     + pa_xxxz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_yyyz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xxx,
                                     double pa_xxxz,
                                     double pa_xz,
                                     double pb_y,
                                     double pb_yyy,
                                     double pb_yyyz,
                                     double pb_yz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * pa_x * fx * fx * pb_y * fz * fgb

                     - 2.25 * pa_x * fz * fga * fx * fx * pb_y

                     - 1.5 * pa_xxx * fx * pb_y * fz * fgb

                     + 11.25 * pa_x * fx * fx * fx * fz * pb_y

                     + 9.0 * pa_xxx * fz * fx * fx * pb_y

                     - 4.5 * pa_xz * fx * pb_yz * fz * fgb

                     - 4.5 * pa_xz * fz * fga * pb_yz * fx

                     - 1.5 * pa_x * fz * fga * fx * pb_yyy

                     - 3.0 * pa_xxxz * pb_yz * fz * fgb

                     + 27.0 * pa_xz * fx * fx * fz * pb_yz

                     + 9.0 * pa_x * fx * fx * fz * pb_yyy

                     + 21.0 * pa_xxxz * fz * pb_yz * fx

                     + 7.0 * pa_xxx * fz * fx * pb_yyy

                     - 3.0 * pa_xz * fz * fga * pb_yyyz

                     + 21.0 * pa_xz * fx * fz * pb_yyyz

                     + 16.0 * pa_xxxz * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_yyzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xxx,
                                     double pa_xxxz,
                                     double pa_xz,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_yyzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * pa_xz * fx * fx * fx

                     + 0.75 * pa_x * fx * fx * fx * pb_z

                     + 0.25 * pa_xxxz * fx * fx

                     + 0.5 * pa_xxx * fx * fx * pb_z

                     + 0.75 * pa_xz * fx * fx * pb_yy

                     + 0.75 * pa_xz * fx * fx * pb_zz

                     + 1.5 * pa_x * fx * fx * pb_yyz

                     + 0.5 * pa_xxxz * pb_yy * fx

                     + 0.5 * pa_xxxz * fx * pb_zz

                     + pa_xxx * fx * pb_yyz

                     + 1.5 * pa_xz * fx * pb_yyzz

                     + pa_xxxz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_yyzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xxx,
                                     double pa_xxxz,
                                     double pa_xz,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_yyzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xz * fx * fx * fz * fgb

                     - 1.5 * pa_x * fx * fx * fz * fgb * pb_z

                     - 0.75 * pa_xz * fz * fga * fx * fx

                     - 1.5 * pa_x * fz * fga * fx * fx * pb_z

                     - pa_xxxz * fx * fz * fgb

                     - pa_xxx * fx * fz * fgb * pb_z

                     + 3.75 * pa_xz * fx * fx * fx * fz

                     + 7.5 * pa_x * fx * fx * fx * fz * pb_z

                     + 3.0 * pa_xxxz * fz * fx * fx

                     + 6.0 * pa_xxx * fz * fx * fx * pb_z

                     - 1.5 * pa_xz * fx * pb_yy * fz * fgb

                     - 1.5 * pa_xz * fx * fz * fgb * pb_zz

                     - 1.5 * pa_xz * fz * fga * pb_yy * fx

                     - 1.5 * pa_xz * fz * fga * fx * pb_zz

                     - 3.0 * pa_x * fz * fga * fx * pb_yyz

                     - pa_xxxz * pb_yy * fz * fgb

                     - pa_xxxz * fz * fgb * pb_zz

                     + 9.0 * pa_xz * fx * fx * fz * pb_yy

                     + 9.0 * pa_xz * fx * fx * fz * pb_zz

                     + 18.0 * pa_x * fx * fx * fz * pb_yyz

                     + 7.0 * pa_xxxz * fz * pb_yy * fx

                     + 7.0 * pa_xxxz * fz * fx * pb_zz

                     + 14.0 * pa_xxx * fz * fx * pb_yyz

                     - 3.0 * pa_xz * fz * fga * pb_yyzz

                     + 21.0 * pa_xz * fx * fz * pb_yyzz

                     + 16.0 * pa_xxxz * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_yzzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xxx,
                                     double pa_xxxz,
                                     double pa_xz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_yzzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_x * fx * fx * fx * pb_y

                     + 0.75 * pa_xxx * fx * fx * pb_y

                     + 2.25 * pa_xz * fx * fx * pb_yz

                     + 2.25 * pa_x * fx * fx * pb_yzz

                     + 1.5 * pa_xxxz * pb_yz * fx

                     + 1.5 * pa_xxx * fx * pb_yzz

                     + 1.5 * pa_xz * fx * pb_yzzz

                     + pa_xxxz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_yzzz_r_0(double fga,
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
                                     double pb_yzzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * pa_x * fx * fx * pb_y * fz * fgb

                     - 2.25 * pa_x * fz * fga * fx * fx * pb_y

                     - 1.5 * pa_xxx * fx * pb_y * fz * fgb

                     + 11.25 * pa_x * fx * fx * fx * fz * pb_y

                     + 9.0 * pa_xxx * fz * fx * fx * pb_y

                     - 4.5 * pa_xz * fx * pb_yz * fz * fgb

                     - 4.5 * pa_xz * fz * fga * pb_yz * fx

                     - 4.5 * pa_x * fz * fga * fx * pb_yzz

                     - 3.0 * pa_xxxz * pb_yz * fz * fgb

                     + 27.0 * pa_xz * fx * fx * fz * pb_yz

                     + 27.0 * pa_x * fx * fx * fz * pb_yzz

                     + 21.0 * pa_xxxz * fz * pb_yz * fx

                     + 21.0 * pa_xxx * fz * fx * pb_yzz

                     - 3.0 * pa_xz * fz * fga * pb_yzzz

                     + 21.0 * pa_xz * fx * fz * pb_yzzz

                     + 16.0 * pa_xxxz * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_zzzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xxx,
                                     double pa_xxxz,
                                     double pa_xz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double pb_zzzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_xz * fx * fx * fx

                     + 4.5 * pa_x * fx * fx * fx * pb_z

                     + 0.75 * pa_xxxz * fx * fx

                     + 3.0 * pa_xxx * fx * fx * pb_z

                     + 4.5 * pa_xz * fx * fx * pb_zz

                     + 3.0 * pa_x * fx * fx * pb_zzz

                     + 3.0 * pa_xxxz * pb_zz * fx

                     + 2.0 * pa_xxx * fx * pb_zzz

                     + 1.5 * pa_xz * fx * pb_zzzz

                     + pa_xxxz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_zzzz_r_0(double fga,
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
                                     double pb_zzzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_xz * fx * fx * fz * fgb

                     - 9.0 * pa_x * fx * fx * pb_z * fz * fgb

                     - 2.25 * pa_xz * fz * fga * fx * fx

                     - 9.0 * pa_x * fz * fga * fx * fx * pb_z

                     - 3.0 * pa_xxxz * fx * fz * fgb

                     - 6.0 * pa_xxx * fx * pb_z * fz * fgb

                     + 11.25 * pa_xz * fx * fx * fx * fz

                     + 45.0 * pa_x * fx * fx * fx * fz * pb_z

                     + 9.0 * pa_xxxz * fz * fx * fx

                     + 36.0 * pa_xxx * fz * fx * fx * pb_z

                     - 9.0 * pa_xz * fx * pb_zz * fz * fgb

                     - 9.0 * pa_xz * fz * fga * pb_zz * fx

                     - 6.0 * pa_x * fz * fga * fx * pb_zzz

                     - 6.0 * pa_xxxz * pb_zz * fz * fgb

                     + 54.0 * pa_xz * fx * fx * fz * pb_zz

                     + 36.0 * pa_x * fx * fx * fz * pb_zzz

                     + 42.0 * pa_xxxz * fz * pb_zz * fx

                     + 28.0 * pa_xxx * fz * fx * pb_zzz

                     - 3.0 * pa_xz * fz * fga * pb_zzzz

                     + 21.0 * pa_xz * fx * fz * pb_zzzz

                     + 16.0 * pa_xxxz * fz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xxxx_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxyy,
                                     double pa_xyy,
                                     double pa_yy,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxx,
                                     double s_0_0)
    {
        return s_0_0 * (0.9375 * fx * fx * fx * fx

                     + 1.875 * fx * fx * fx * pa_yy

                     + 0.375 * pa_xx * fx * fx * fx

                     + 3.0 * pa_x * fx * fx * fx * pb_x

                     + 2.25 * fx * fx * fx * pb_xx

                     + 0.75 * pa_xxyy * fx * fx

                     + 6.0 * pa_xyy * fx * fx * pb_x

                     + 4.5 * fx * fx * pa_yy * pb_xx

                     + 1.5 * pa_xx * fx * fx * pb_xx

                     + 2.0 * pa_x * fx * fx * pb_xxx

                     + 3.0 * pa_xxyy * pb_xx * fx

                     + 4.0 * pa_xyy * fx * pb_xxx

                     + 0.25 * fx * fx * pb_xxxx

                     + 0.5 * pa_xx * fx * pb_xxxx

                     + 0.5 * fx * pa_yy * pb_xxxx

                     + pa_xxyy * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xxxx_r_0(double fga,
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
                                     double pb_xxxx,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * fx * fz * fgb

                     - 2.25 * fx * fx * fx * fz * fga

                     - 4.5 * fx * fx * pa_yy * fz * fgb

                     + 7.5 * fx * fx * fx * fx * fz

                     + 18.75 * fx * fx * fx * pa_yy * fz

                     - 1.5 * pa_xx * fx * fx * fz * fgb

                     - 0.75 * pa_xx * fz * fga * fx * fx

                     - 6.0 * pa_x * fx * fx * pb_x * fz * fgb

                     - 6.0 * pa_x * fx * fx * fz * fga * pb_x

                     - 6.0 * fx * fx * fz * fga * pb_xx

                     - 0.75 * fz * fga * pa_yy * fx * fx

                     - 3.0 * pa_xxyy * fx * fz * fgb

                     + 3.75 * pa_xx * fz * fx * fx * fx

                     - 12.0 * pa_xyy * fx * pb_x * fz * fgb

                     + 30.0 * pa_x * fx * fx * fx * fz * pb_x

                     + 22.5 * fx * fx * fx * fz * pb_xx

                     + 9.0 * pa_xxyy * fz * fx * fx

                     + 72.0 * pa_xyy * fx * fx * fz * pb_x

                     + 54.0 * fx * fx * pa_yy * fz * pb_xx

                     - 1.5 * fx * fx * pb_xx * fz * fgb

                     - 3.0 * pa_xx * fx * pb_xx * fz * fgb

                     - 3.0 * pa_xx * fz * fga * pb_xx * fx

                     - 4.0 * pa_x * fx * fz * fga * pb_xxx

                     - 3.0 * fx * pa_yy * pb_xx * fz * fgb

                     - 3.0 * fz * fga * pa_yy * pb_xx * fx

                     - 6.0 * pa_xxyy * pb_xx * fz * fgb

                     + 18.0 * pa_xx * fz * fx * fx * pb_xx

                     + 24.0 * pa_x * fx * fx * fz * pb_xxx

                     + 42.0 * pa_xxyy * fz * pb_xx * fx

                     + 56.0 * pa_xyy * fx * fz * pb_xxx

                     - fx * fz * fga * pb_xxxx

                     - pa_xx * fz * fga * pb_xxxx

                     + 3.0 * fx * fx * fz * pb_xxxx

                     - fz * fga * pa_yy * pb_xxxx

                     + 7.0 * pa_xx * fz * fx * pb_xxxx

                     + 7.0 * fx * pa_yy * fz * pb_xxxx

                     + 16.0 * pa_xxyy * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xxxy_s_0(double fx,
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
                                     double pb_xxx,
                                     double pb_xxxy,
                                     double pb_xxy,
                                     double pb_xy,
                                     double pb_y,
                                     double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xy * fx * fx * fx

                     + 2.25 * fx * fx * fx * pa_y * pb_x

                     + 0.75 * pa_x * fx * fx * fx * pb_y

                     + 1.125 * fx * fx * fx * pb_xy

                     + 1.5 * pa_xxy * fx * fx * pb_x

                     + 1.5 * pa_xyy * fx * fx * pb_y

                     + 3.0 * pa_xy * fx * fx * pb_xx

                     + 2.25 * fx * fx * pa_yy * pb_xy

                     + 0.75 * pa_xx * fx * fx * pb_xy

                     + 1.5 * pa_x * fx * fx * pb_xxy

                     + 0.5 * fx * fx * pa_y * pb_xxx

                     + 1.5 * pa_xxyy * pb_xy * fx

                     + pa_xxy * fx * pb_xxx

                     + 3.0 * pa_xyy * fx * pb_xxy

                     + 0.25 * fx * fx * pb_xxxy

                     + 0.5 * pa_xx * fx * pb_xxxy

                     + 0.5 * fx * pa_yy * pb_xxxy

                     + pa_xxyy * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xxxy_r_0(double fga,
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
                                     double pb_xxx,
                                     double pb_xxxy,
                                     double pb_xxy,
                                     double pb_xy,
                                     double pb_y,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xy * fx * fx * fz * fgb

                     + 15.0 * pa_xy * fx * fx * fx * fz

                     + 22.5 * fx * fx * fx * pa_y * fz * pb_x

                     - 1.5 * pa_x * fx * fx * fz * fgb * pb_y

                     - 1.5 * pa_x * fx * fx * fz * fga * pb_y

                     - 1.5 * fx * fx * pa_y * pb_x * fz * fgb

                     - 3.0 * fx * fx * fz * fga * pb_xy

                     - 1.5 * fz * fga * pa_y * fx * fx * pb_x

                     - 3.0 * pa_xxy * fx * pb_x * fz * fgb

                     - 3.0 * pa_xyy * fx * fz * fgb * pb_y

                     + 7.5 * pa_x * fx * fx * fx * fz * pb_y

                     + 11.25 * fx * fx * fx * fz * pb_xy

                     + 18.0 * pa_xxy * fz * fx * fx * pb_x

                     + 18.0 * pa_xyy * fx * fx * fz * pb_y

                     + 36.0 * pa_xy * fx * fx * fz * pb_xx

                     + 27.0 * fx * fx * pa_yy * fz * pb_xy

                     - 0.75 * fx * fx * pb_xy * fz * fgb

                     - 1.5 * pa_xx * fx * pb_xy * fz * fgb

                     - 1.5 * pa_xx * fz * fga * pb_xy * fx

                     - 3.0 * pa_x * fx * fz * fga * pb_xxy

                     - 1.5 * fx * pa_yy * pb_xy * fz * fgb

                     - 1.5 * fz * fga * pa_yy * pb_xy * fx

                     - fz * fga * pa_y * fx * pb_xxx

                     - 3.0 * pa_xxyy * pb_xy * fz * fgb

                     + 9.0 * pa_xx * fz * fx * fx * pb_xy

                     + 18.0 * pa_x * fx * fx * fz * pb_xxy

                     + 6.0 * fx * fx * pa_y * fz * pb_xxx

                     + 21.0 * pa_xxyy * fz * pb_xy * fx

                     + 14.0 * pa_xxy * fz * fx * pb_xxx

                     + 42.0 * pa_xyy * fx * fz * pb_xxy

                     - fx * fz * fga * pb_xxxy

                     - pa_xx * fz * fga * pb_xxxy

                     + 3.0 * fx * fx * fz * pb_xxxy

                     - fz * fga * pa_yy * pb_xxxy

                     + 7.0 * pa_xx * fz * fx * pb_xxxy

                     + 7.0 * fx * pa_yy * fz * pb_xxxy

                     + 16.0 * pa_xxyy * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xxxz_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxyy,
                                     double pa_xyy,
                                     double pa_yy,
                                     double pb_xxxz,
                                     double pb_xxz,
                                     double pb_xz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * fx * pb_z

                     + 1.125 * fx * fx * fx * pb_xz

                     + 1.5 * pa_xyy * fx * fx * pb_z

                     + 2.25 * fx * fx * pa_yy * pb_xz

                     + 0.75 * pa_xx * fx * fx * pb_xz

                     + 1.5 * pa_x * fx * fx * pb_xxz

                     + 1.5 * pa_xxyy * pb_xz * fx

                     + 3.0 * pa_xyy * fx * pb_xxz

                     + 0.25 * fx * fx * pb_xxxz

                     + 0.5 * pa_xx * fx * pb_xxxz

                     + 0.5 * fx * pa_yy * pb_xxxz

                     + pa_xxyy * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xxxz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxyy,
                                     double pa_xyy,
                                     double pa_yy,
                                     double pb_xxxz,
                                     double pb_xxz,
                                     double pb_xz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fx * fz * fgb * pb_z

                     - 1.5 * pa_x * fx * fx * fz * fga * pb_z

                     - 3.0 * fx * fx * fz * fga * pb_xz

                     - 3.0 * pa_xyy * fx * fz * fgb * pb_z

                     + 7.5 * pa_x * fx * fx * fx * fz * pb_z

                     + 11.25 * fx * fx * fx * fz * pb_xz

                     + 18.0 * pa_xyy * fx * fx * fz * pb_z

                     + 27.0 * fx * fx * pa_yy * fz * pb_xz

                     - 0.75 * fx * fx * pb_xz * fz * fgb

                     - 1.5 * pa_xx * fx * pb_xz * fz * fgb

                     - 1.5 * pa_xx * fz * fga * pb_xz * fx

                     - 3.0 * pa_x * fx * fz * fga * pb_xxz

                     - 1.5 * fx * pa_yy * pb_xz * fz * fgb

                     - 1.5 * fz * fga * pa_yy * pb_xz * fx

                     - 3.0 * pa_xxyy * pb_xz * fz * fgb

                     + 9.0 * pa_xx * fz * fx * fx * pb_xz

                     + 18.0 * pa_x * fx * fx * fz * pb_xxz

                     + 21.0 * pa_xxyy * fz * pb_xz * fx

                     + 42.0 * pa_xyy * fx * fz * pb_xxz

                     - fx * fz * fga * pb_xxxz

                     - pa_xx * fz * fga * pb_xxxz

                     + 3.0 * fx * fx * fz * pb_xxxz

                     - fz * fga * pa_yy * pb_xxxz

                     + 7.0 * pa_xx * fz * fx * pb_xxxz

                     + 7.0 * fx * pa_yy * fz * pb_xxxz

                     + 16.0 * pa_xxyy * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xxyy_s_0(double fx,
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
                                     double pb_xxyy,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_y,
                                     double pb_yy,
                                     double s_0_0)
    {
        return s_0_0 * (0.5625 * fx * fx * fx * fx

                     + 0.375 * pa_xx * fx * fx * fx

                     + 1.5 * pa_x * fx * fx * fx * pb_x

                     + 0.375 * fx * fx * fx * pa_yy

                     + 1.5 * fx * fx * fx * pa_y * pb_y

                     + 0.375 * fx * fx * fx * pb_xx

                     + 0.375 * fx * fx * fx * pb_yy

                     + 0.25 * pa_xxyy * fx * fx

                     + pa_xxy * fx * fx * pb_y

                     + 0.75 * pa_xx * fx * fx * pb_xx

                     + pa_xyy * fx * fx * pb_x

                     + 4.0 * pa_xy * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_yy * pb_yy

                     + 0.25 * pa_xx * fx * fx * pb_yy

                     + pa_x * fx * fx * pb_xyy

                     + 0.25 * fx * fx * pa_yy * pb_xx

                     + fx * fx * pa_y * pb_xxy

                     + 0.5 * pa_xxyy * pb_xx * fx

                     + 0.5 * pa_xxyy * fx * pb_yy

                     + 2.0 * pa_xxy * fx * pb_xxy

                     + 2.0 * pa_xyy * fx * pb_xyy

                     + 0.25 * fx * fx * pb_xxyy

                     + 0.5 * pa_xx * fx * pb_xxyy

                     + 0.5 * fx * pa_yy * pb_xxyy

                     + pa_xxyy * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xxyy_r_0(double fga,
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
                                     double pb_xxyy,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_y,
                                     double pb_yy,
                                     double r_0_0)
    {
        return r_0_0 * (4.5 * fx * fx * fx * fx * fz

                     - 0.75 * fx * fx * fx * fz * fgb

                     - 0.75 * fx * fx * fx * fz * fga

                     - pa_xx * fx * fx * fz * fgb

                     - fx * fx * pa_yy * fz * fgb

                     + 3.75 * pa_xx * fx * fx * fx * fz

                     + 15.0 * pa_x * fx * fx * fx * fz * pb_x

                     + 3.75 * fx * fx * fx * pa_yy * fz

                     + 15.0 * fx * fx * fx * pa_y * fz * pb_y

                     - 0.25 * pa_xx * fz * fga * fx * fx

                     - pa_x * fx * fx * pb_x * fz * fgb

                     - pa_x * fx * fx * fz * fga * pb_x

                     - fx * fx * pa_y * fz * fgb * pb_y

                     - fx * fx * fz * fga * pb_yy

                     - 0.25 * fz * fga * pa_yy * fx * fx

                     - fz * fga * pa_y * fx * fx * pb_y

                     - fz * fga * fx * fx * pb_xx

                     - pa_xxyy * fx * fz * fgb

                     - 2.0 * pa_xxy * fx * fz * fgb * pb_y

                     - 2.0 * pa_xyy * fx * pb_x * fz * fgb

                     + 3.75 * fx * fx * fx * fz * pb_xx

                     + 3.75 * fx * fx * fx * fz * pb_yy

                     + 3.0 * pa_xxyy * fz * fx * fx

                     + 12.0 * pa_xxy * fz * fx * fx * pb_y

                     + 9.0 * pa_xx * fx * fx * fz * pb_xx

                     + 12.0 * pa_xyy * fx * fx * fz * pb_x

                     + 48.0 * pa_xy * fx * fx * fz * pb_xy

                     + 9.0 * fx * fx * pa_yy * fz * pb_yy

                     - 0.25 * fx * fx * pb_xx * fz * fgb

                     - 0.25 * fx * fx * fz * fgb * pb_yy

                     - 0.5 * pa_xx * fx * pb_xx * fz * fgb

                     - 0.5 * pa_xx * fx * fz * fgb * pb_yy

                     - 0.5 * pa_xx * fz * fga * pb_xx * fx

                     - 0.5 * pa_xx * fz * fga * fx * pb_yy

                     - 2.0 * pa_x * fx * fz * fga * pb_xyy

                     - 0.5 * fx * pa_yy * pb_xx * fz * fgb

                     - 0.5 * fx * pa_yy * fz * fgb * pb_yy

                     - 0.5 * fz * fga * pa_yy * pb_xx * fx

                     - 0.5 * fz * fga * pa_yy * fx * pb_yy

                     - 2.0 * fz * fga * pa_y * fx * pb_xxy

                     - pa_xxyy * pb_xx * fz * fgb

                     - pa_xxyy * fz * fgb * pb_yy

                     + 3.0 * pa_xx * fz * fx * fx * pb_yy

                     + 12.0 * pa_x * fx * fx * fz * pb_xyy

                     + 3.0 * fx * fx * pa_yy * fz * pb_xx

                     + 12.0 * fx * fx * pa_y * fz * pb_xxy

                     + 7.0 * pa_xxyy * fz * pb_xx * fx

                     + 7.0 * pa_xxyy * fz * fx * pb_yy

                     + 28.0 * pa_xxy * fz * fx * pb_xxy

                     + 28.0 * pa_xyy * fx * fz * pb_xyy

                     - fx * fz * fga * pb_xxyy

                     - pa_xx * fz * fga * pb_xxyy

                     + 3.0 * fx * fx * fz * pb_xxyy

                     - fz * fga * pa_yy * pb_xxyy

                     + 7.0 * pa_xx * fz * fx * pb_xxyy

                     + 7.0 * fx * pa_yy * fz * pb_xxyy

                     + 16.0 * pa_xxyy * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xxyz_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxy,
                                     double pa_xxyy,
                                     double pa_xy,
                                     double pa_xyy,
                                     double pa_y,
                                     double pa_yy,
                                     double pb_xxyz,
                                     double pb_xxz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * fx * pa_y * pb_z

                     + 0.375 * fx * fx * fx * pb_yz

                     + 0.5 * pa_xxy * fx * fx * pb_z

                     + 2.0 * pa_xy * fx * fx * pb_xz

                     + 0.75 * fx * fx * pa_yy * pb_yz

                     + 0.25 * pa_xx * fx * fx * pb_yz

                     + pa_x * fx * fx * pb_xyz

                     + 0.5 * fx * fx * pa_y * pb_xxz

                     + 0.5 * pa_xxyy * fx * pb_yz

                     + pa_xxy * fx * pb_xxz

                     + 2.0 * pa_xyy * fx * pb_xyz

                     + 0.25 * fx * fx * pb_xxyz

                     + 0.5 * pa_xx * fx * pb_xxyz

                     + 0.5 * fx * pa_yy * pb_xxyz

                     + pa_xxyy * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xxyz_r_0(double fga,
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
                                     double pb_xxyz,
                                     double pb_xxz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (7.5 * fx * fx * fx * pa_y * fz * pb_z

                     - 0.5 * fx * fx * pa_y * fz * fgb * pb_z

                     - fx * fx * fz * fga * pb_yz

                     - 0.5 * fz * fga * pa_y * fx * fx * pb_z

                     - pa_xxy * fx * fz * fgb * pb_z

                     + 3.75 * fx * fx * fx * fz * pb_yz

                     + 6.0 * pa_xxy * fz * fx * fx * pb_z

                     + 24.0 * pa_xy * fx * fx * fz * pb_xz

                     + 9.0 * fx * fx * pa_yy * fz * pb_yz

                     - 0.25 * fx * fx * fz * fgb * pb_yz

                     - 0.5 * pa_xx * fx * fz * fgb * pb_yz

                     - 0.5 * pa_xx * fz * fga * fx * pb_yz

                     - 2.0 * pa_x * fx * fz * fga * pb_xyz

                     - 0.5 * fx * pa_yy * fz * fgb * pb_yz

                     - 0.5 * fz * fga * pa_yy * fx * pb_yz

                     - fz * fga * pa_y * fx * pb_xxz

                     - pa_xxyy * fz * fgb * pb_yz

                     + 3.0 * pa_xx * fz * fx * fx * pb_yz

                     + 12.0 * pa_x * fx * fx * fz * pb_xyz

                     + 6.0 * fx * fx * pa_y * fz * pb_xxz

                     + 7.0 * pa_xxyy * fz * fx * pb_yz

                     + 14.0 * pa_xxy * fz * fx * pb_xxz

                     + 28.0 * pa_xyy * fx * fz * pb_xyz

                     - fx * fz * fga * pb_xxyz

                     - pa_xx * fz * fga * pb_xxyz

                     + 3.0 * fx * fx * fz * pb_xxyz

                     - fz * fga * pa_yy * pb_xxyz

                     + 7.0 * pa_xx * fz * fx * pb_xxyz

                     + 7.0 * fx * pa_yy * fz * pb_xxyz

                     + 16.0 * pa_xxyy * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xxzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxyy,
                                     double pa_xyy,
                                     double pa_yy,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxzz,
                                     double pb_xzz,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * fx * fx * fx * pa_yy

                     + 0.125 * pa_xx * fx * fx * fx

                     + 0.5 * pa_x * fx * fx * fx * pb_x

                     + 0.375 * fx * fx * fx * pb_zz

                     + 0.25 * pa_xxyy * fx * fx

                     + pa_xyy * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_yy * pb_zz

                     + 0.125 * fx * fx * fx * pb_xx

                     + 0.25 * pa_xx * fx * fx * pb_xx

                     + 0.25 * pa_xx * fx * fx * pb_zz

                     + pa_x * fx * fx * pb_xzz

                     + 0.25 * fx * fx * pa_yy * pb_xx

                     + 0.5 * pa_xxyy * pb_xx * fx

                     + 0.5 * pa_xxyy * fx * pb_zz

                     + 2.0 * pa_xyy * fx * pb_xzz

                     + 0.25 * fx * fx * pb_xxzz

                     + 0.5 * pa_xx * fx * pb_xxzz

                     + 0.5 * fx * pa_yy * pb_xxzz

                     + pa_xxyy * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xxzz_r_0(double fga,
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
                                     double pb_xxzz,
                                     double pb_xzz,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * fx * fz * fgb

                     - 0.5 * fx * fx * fx * fz * fga

                     - fx * fx * pa_yy * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     + 3.75 * fx * fx * fx * pa_yy * fz

                     - 0.5 * pa_xx * fx * fx * fz * fgb

                     - 0.25 * pa_xx * fz * fga * fx * fx

                     - pa_x * fx * fx * pb_x * fz * fgb

                     - pa_x * fx * fx * fz * fga * pb_x

                     - fx * fx * fz * fga * pb_zz

                     - 0.25 * fz * fga * pa_yy * fx * fx

                     - pa_xxyy * fx * fz * fgb

                     + 1.25 * pa_xx * fz * fx * fx * fx

                     - 2.0 * pa_xyy * fx * pb_x * fz * fgb

                     + 5.0 * pa_x * fx * fx * fx * fz * pb_x

                     + 3.75 * fx * fx * fx * fz * pb_zz

                     + 3.0 * pa_xxyy * fz * fx * fx

                     + 12.0 * pa_xyy * fx * fx * fz * pb_x

                     + 9.0 * fx * fx * pa_yy * fz * pb_zz

                     - 0.25 * fx * fx * pb_xx * fz * fgb

                     - 0.25 * fx * fx * fz * fgb * pb_zz

                     - 0.5 * fx * fx * fz * fga * pb_xx

                     - 0.5 * pa_xx * fx * pb_xx * fz * fgb

                     - 0.5 * pa_xx * fx * fz * fgb * pb_zz

                     - 0.5 * pa_xx * fz * fga * pb_xx * fx

                     - 0.5 * pa_xx * fz * fga * fx * pb_zz

                     - 2.0 * pa_x * fx * fz * fga * pb_xzz

                     - 0.5 * fx * pa_yy * pb_xx * fz * fgb

                     - 0.5 * fx * pa_yy * fz * fgb * pb_zz

                     + 1.25 * fx * fx * fx * fz * pb_xx

                     - 0.5 * fz * fga * pa_yy * pb_xx * fx

                     - 0.5 * fz * fga * pa_yy * fx * pb_zz

                     - pa_xxyy * pb_xx * fz * fgb

                     - pa_xxyy * fz * fgb * pb_zz

                     + 3.0 * pa_xx * fz * fx * fx * pb_xx

                     + 3.0 * pa_xx * fz * fx * fx * pb_zz

                     + 12.0 * pa_x * fx * fx * fz * pb_xzz

                     + 3.0 * fx * fx * pa_yy * fz * pb_xx

                     + 7.0 * pa_xxyy * fz * pb_xx * fx

                     + 7.0 * pa_xxyy * fz * fx * pb_zz

                     + 28.0 * pa_xyy * fx * fz * pb_xzz

                     - fx * fz * fga * pb_xxzz

                     - pa_xx * fz * fga * pb_xxzz

                     + 3.0 * fx * fx * fz * pb_xxzz

                     - fz * fga * pa_yy * pb_xxzz

                     + 7.0 * pa_xx * fz * fx * pb_xxzz

                     + 7.0 * fx * pa_yy * fz * pb_xxzz

                     + 16.0 * pa_xxyy * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xyyy_s_0(double fx,
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
                                     double pb_xyyy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xy * fx * fx * fx

                     + 2.25 * pa_x * fx * fx * fx * pb_y

                     + 0.75 * fx * fx * fx * pa_y * pb_x

                     + 1.125 * fx * fx * fx * pb_xy

                     + 1.5 * pa_xxy * fx * fx * pb_x

                     + 2.25 * pa_xx * fx * fx * pb_xy

                     + 1.5 * pa_xyy * fx * fx * pb_y

                     + 3.0 * pa_xy * fx * fx * pb_yy

                     + 0.5 * pa_x * fx * fx * pb_yyy

                     + 0.75 * fx * fx * pa_yy * pb_xy

                     + 1.5 * fx * fx * pa_y * pb_xyy

                     + 1.5 * pa_xxyy * pb_xy * fx

                     + 3.0 * pa_xxy * fx * pb_xyy

                     + pa_xyy * fx * pb_yyy

                     + 0.25 * fx * fx * pb_xyyy

                     + 0.5 * pa_xx * fx * pb_xyyy

                     + 0.5 * fx * pa_yy * pb_xyyy

                     + pa_xxyy * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xyyy_r_0(double fga,
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
                                     double pb_xyyy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xy * fx * fx * fz * fgb

                     + 15.0 * pa_xy * fx * fx * fx * fz

                     + 22.5 * pa_x * fx * fx * fx * fz * pb_y

                     - 1.5 * pa_x * fx * fx * pb_y * fz * fgb

                     - 1.5 * pa_x * fx * fx * fz * fga * pb_y

                     - 1.5 * fx * fx * pa_y * pb_x * fz * fgb

                     - 1.5 * fz * fga * pa_y * fx * fx * pb_x

                     - 3.0 * fz * fga * fx * fx * pb_xy

                     - 3.0 * pa_xxy * fx * pb_x * fz * fgb

                     - 3.0 * pa_xyy * fx * pb_y * fz * fgb

                     + 7.5 * fx * fx * fx * pa_y * fz * pb_x

                     + 11.25 * fx * fx * fx * fz * pb_xy

                     + 18.0 * pa_xxy * fz * fx * fx * pb_x

                     + 27.0 * pa_xx * fx * fx * fz * pb_xy

                     + 18.0 * pa_xyy * fx * fx * fz * pb_y

                     + 36.0 * pa_xy * fx * fx * fz * pb_yy

                     - 0.75 * fx * fx * pb_xy * fz * fgb

                     - 1.5 * pa_xx * fx * pb_xy * fz * fgb

                     - 1.5 * pa_xx * fz * fga * pb_xy * fx

                     - pa_x * fx * fz * fga * pb_yyy

                     - 1.5 * fx * pa_yy * pb_xy * fz * fgb

                     - 1.5 * fz * fga * pa_yy * pb_xy * fx

                     - 3.0 * fz * fga * pa_y * fx * pb_xyy

                     - 3.0 * pa_xxyy * pb_xy * fz * fgb

                     + 6.0 * pa_x * fx * fx * fz * pb_yyy

                     + 9.0 * fx * fx * pa_yy * fz * pb_xy

                     + 18.0 * fx * fx * pa_y * fz * pb_xyy

                     + 21.0 * pa_xxyy * fz * pb_xy * fx

                     + 42.0 * pa_xxy * fz * fx * pb_xyy

                     + 14.0 * pa_xyy * fx * fz * pb_yyy

                     - fx * fz * fga * pb_xyyy

                     - pa_xx * fz * fga * pb_xyyy

                     + 3.0 * fx * fx * fz * pb_xyyy

                     - fz * fga * pa_yy * pb_xyyy

                     + 7.0 * pa_xx * fz * fx * pb_xyyy

                     + 7.0 * fx * pa_yy * fz * pb_xyyy

                     + 16.0 * pa_xxyy * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xyyz_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxy,
                                     double pa_xxyy,
                                     double pa_xy,
                                     double pa_xyy,
                                     double pa_y,
                                     double pa_yy,
                                     double pb_xyyz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * fx * pb_z

                     + 0.375 * fx * fx * fx * pb_xz

                     + 0.75 * pa_xx * fx * fx * pb_xz

                     + 0.5 * pa_xyy * fx * fx * pb_z

                     + 2.0 * pa_xy * fx * fx * pb_yz

                     + 0.5 * pa_x * fx * fx * pb_yyz

                     + 0.25 * fx * fx * pa_yy * pb_xz

                     + fx * fx * pa_y * pb_xyz

                     + 0.5 * pa_xxyy * pb_xz * fx

                     + 2.0 * pa_xxy * fx * pb_xyz

                     + pa_xyy * fx * pb_yyz

                     + 0.25 * fx * fx * pb_xyyz

                     + 0.5 * pa_xx * fx * pb_xyyz

                     + 0.5 * fx * pa_yy * pb_xyyz

                     + pa_xxyy * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xyyz_r_0(double fga,
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
                                     double pb_xyyz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (7.5 * pa_x * fx * fx * fx * fz * pb_z

                     - 0.5 * pa_x * fx * fx * fz * fgb * pb_z

                     - 0.5 * pa_x * fx * fx * fz * fga * pb_z

                     - fz * fga * fx * fx * pb_xz

                     - pa_xyy * fx * fz * fgb * pb_z

                     + 3.75 * fx * fx * fx * fz * pb_xz

                     + 9.0 * pa_xx * fx * fx * fz * pb_xz

                     + 6.0 * pa_xyy * fx * fx * fz * pb_z

                     + 24.0 * pa_xy * fx * fx * fz * pb_yz

                     - 0.25 * fx * fx * pb_xz * fz * fgb

                     - 0.5 * pa_xx * fx * pb_xz * fz * fgb

                     - 0.5 * pa_xx * fz * fga * pb_xz * fx

                     - pa_x * fx * fz * fga * pb_yyz

                     - 0.5 * fx * pa_yy * pb_xz * fz * fgb

                     - 0.5 * fz * fga * pa_yy * pb_xz * fx

                     - 2.0 * fz * fga * pa_y * fx * pb_xyz

                     - pa_xxyy * pb_xz * fz * fgb

                     + 6.0 * pa_x * fx * fx * fz * pb_yyz

                     + 3.0 * fx * fx * pa_yy * fz * pb_xz

                     + 12.0 * fx * fx * pa_y * fz * pb_xyz

                     + 7.0 * pa_xxyy * fz * pb_xz * fx

                     + 28.0 * pa_xxy * fz * fx * pb_xyz

                     + 14.0 * pa_xyy * fx * fz * pb_yyz

                     - fx * fz * fga * pb_xyyz

                     - pa_xx * fz * fga * pb_xyyz

                     + 3.0 * fx * fx * fz * pb_xyyz

                     - fz * fga * pa_yy * pb_xyyz

                     + 7.0 * pa_xx * fz * fx * pb_xyyz

                     + 7.0 * fx * pa_yy * fz * pb_xyyz

                     + 16.0 * pa_xxyy * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xyzz_s_0(double fx,
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
                                     double pb_xyzz,
                                     double pb_xzz,
                                     double pb_y,
                                     double pb_yzz,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xy * fx * fx * fx

                     + 0.25 * pa_x * fx * fx * fx * pb_y

                     + 0.25 * fx * fx * fx * pa_y * pb_x

                     + 0.5 * pa_xxy * fx * fx * pb_x

                     + 0.5 * pa_xyy * fx * fx * pb_y

                     + pa_xy * fx * fx * pb_zz

                     + 0.125 * fx * fx * fx * pb_xy

                     + 0.25 * pa_xx * fx * fx * pb_xy

                     + 0.5 * pa_x * fx * fx * pb_yzz

                     + 0.25 * fx * fx * pa_yy * pb_xy

                     + 0.5 * fx * fx * pa_y * pb_xzz

                     + 0.5 * pa_xxyy * pb_xy * fx

                     + pa_xxy * fx * pb_xzz

                     + pa_xyy * fx * pb_yzz

                     + 0.25 * fx * fx * pb_xyzz

                     + 0.5 * pa_xx * fx * pb_xyzz

                     + 0.5 * fx * pa_yy * pb_xyzz

                     + pa_xxyy * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xyzz_r_0(double fga,
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
                                     double pb_xyzz,
                                     double pb_xzz,
                                     double pb_y,
                                     double pb_yzz,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- pa_xy * fx * fx * fz * fgb

                     + 5.0 * pa_xy * fx * fx * fx * fz

                     - 0.5 * pa_x * fx * fx * pb_y * fz * fgb

                     - 0.5 * pa_x * fx * fx * fz * fga * pb_y

                     - 0.5 * fx * fx * pa_y * pb_x * fz * fgb

                     - 0.5 * fz * fga * pa_y * fx * fx * pb_x

                     - pa_xxy * fx * pb_x * fz * fgb

                     - pa_xyy * fx * pb_y * fz * fgb

                     + 2.5 * pa_x * fx * fx * fx * fz * pb_y

                     + 2.5 * fx * fx * fx * pa_y * fz * pb_x

                     + 6.0 * pa_xxy * fz * fx * fx * pb_x

                     + 6.0 * pa_xyy * fx * fx * fz * pb_y

                     + 12.0 * pa_xy * fx * fx * fz * pb_zz

                     - 0.25 * fx * fx * pb_xy * fz * fgb

                     - 0.5 * fx * fx * fz * fga * pb_xy

                     - 0.5 * pa_xx * fx * pb_xy * fz * fgb

                     - 0.5 * pa_xx * fz * fga * pb_xy * fx

                     - pa_x * fx * fz * fga * pb_yzz

                     - 0.5 * fx * pa_yy * pb_xy * fz * fgb

                     + 1.25 * fx * fx * fx * fz * pb_xy

                     - 0.5 * fz * fga * pa_yy * pb_xy * fx

                     - fz * fga * pa_y * fx * pb_xzz

                     - pa_xxyy * pb_xy * fz * fgb

                     + 3.0 * pa_xx * fz * fx * fx * pb_xy

                     + 6.0 * pa_x * fx * fx * fz * pb_yzz

                     + 3.0 * fx * fx * pa_yy * fz * pb_xy

                     + 6.0 * fx * fx * pa_y * fz * pb_xzz

                     + 7.0 * pa_xxyy * fz * pb_xy * fx

                     + 14.0 * pa_xxy * fz * fx * pb_xzz

                     + 14.0 * pa_xyy * fx * fz * pb_yzz

                     - fx * fz * fga * pb_xyzz

                     - pa_xx * fz * fga * pb_xyzz

                     + 3.0 * fx * fx * fz * pb_xyzz

                     - fz * fga * pa_yy * pb_xyzz

                     + 7.0 * pa_xx * fz * fx * pb_xyzz

                     + 7.0 * fx * pa_yy * fz * pb_xyzz

                     + 16.0 * pa_xxyy * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xzzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxyy,
                                     double pa_xyy,
                                     double pa_yy,
                                     double pb_xz,
                                     double pb_xzzz,
                                     double pb_z,
                                     double pb_zzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * fx * pb_z

                     + 1.5 * pa_xyy * fx * fx * pb_z

                     + 0.375 * fx * fx * fx * pb_xz

                     + 0.75 * pa_xx * fx * fx * pb_xz

                     + 0.5 * pa_x * fx * fx * pb_zzz

                     + 0.75 * fx * fx * pa_yy * pb_xz

                     + 1.5 * pa_xxyy * pb_xz * fx

                     + pa_xyy * fx * pb_zzz

                     + 0.25 * fx * fx * pb_xzzz

                     + 0.5 * pa_xx * fx * pb_xzzz

                     + 0.5 * fx * pa_yy * pb_xzzz

                     + pa_xxyy * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xzzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxyy,
                                     double pa_xyy,
                                     double pa_yy,
                                     double pb_xz,
                                     double pb_xzzz,
                                     double pb_z,
                                     double pb_zzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fx * pb_z * fz * fgb

                     - 1.5 * pa_x * fx * fx * fz * fga * pb_z

                     - 3.0 * pa_xyy * fx * pb_z * fz * fgb

                     + 7.5 * pa_x * fx * fx * fx * fz * pb_z

                     + 18.0 * pa_xyy * fx * fx * fz * pb_z

                     - 0.75 * fx * fx * pb_xz * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_xz

                     - 1.5 * pa_xx * fx * pb_xz * fz * fgb

                     - 1.5 * pa_xx * fz * fga * pb_xz * fx

                     - pa_x * fx * fz * fga * pb_zzz

                     - 1.5 * fx * pa_yy * pb_xz * fz * fgb

                     + 3.75 * fx * fx * fx * fz * pb_xz

                     - 1.5 * fz * fga * pa_yy * pb_xz * fx

                     - 3.0 * pa_xxyy * pb_xz * fz * fgb

                     + 9.0 * pa_xx * fz * fx * fx * pb_xz

                     + 6.0 * pa_x * fx * fx * fz * pb_zzz

                     + 9.0 * fx * fx * pa_yy * fz * pb_xz

                     + 21.0 * pa_xxyy * fz * pb_xz * fx

                     + 14.0 * pa_xyy * fx * fz * pb_zzz

                     - fx * fz * fga * pb_xzzz

                     - pa_xx * fz * fga * pb_xzzz

                     + 3.0 * fx * fx * fz * pb_xzzz

                     - fz * fga * pa_yy * pb_xzzz

                     + 7.0 * pa_xx * fz * fx * pb_xzzz

                     + 7.0 * fx * pa_yy * fz * pb_xzzz

                     + 16.0 * pa_xxyy * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_yyyy_s_0(double fx,
                                     double pa_xx,
                                     double pa_xxy,
                                     double pa_xxyy,
                                     double pa_y,
                                     double pa_yy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double pb_yyyy,
                                     double s_0_0)
    {
        return s_0_0 * (0.9375 * fx * fx * fx * fx

                     + 1.875 * pa_xx * fx * fx * fx

                     + 0.375 * fx * fx * fx * pa_yy

                     + 3.0 * fx * fx * fx * pa_y * pb_y

                     + 2.25 * fx * fx * fx * pb_yy

                     + 0.75 * pa_xxyy * fx * fx

                     + 6.0 * pa_xxy * fx * fx * pb_y

                     + 4.5 * pa_xx * fx * fx * pb_yy

                     + 1.5 * fx * fx * pa_yy * pb_yy

                     + 2.0 * fx * fx * pa_y * pb_yyy

                     + 3.0 * pa_xxyy * pb_yy * fx

                     + 4.0 * pa_xxy * fx * pb_yyy

                     + 0.25 * fx * fx * pb_yyyy

                     + 0.5 * pa_xx * fx * pb_yyyy

                     + 0.5 * fx * pa_yy * pb_yyyy

                     + pa_xxyy * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_yyyy_r_0(double fga,
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
                                     double pb_yyyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * fx * fz * fgb

                     - 2.25 * fz * fga * fx * fx * fx

                     - 4.5 * pa_xx * fx * fx * fz * fgb

                     + 7.5 * fx * fx * fx * fx * fz

                     + 18.75 * pa_xx * fx * fx * fx * fz

                     - 0.75 * pa_xx * fz * fga * fx * fx

                     - 1.5 * fx * fx * pa_yy * fz * fgb

                     - 6.0 * fx * fx * pa_y * pb_y * fz * fgb

                     - 0.75 * fz * fga * pa_yy * fx * fx

                     - 6.0 * fz * fga * pa_y * fx * fx * pb_y

                     - 6.0 * fz * fga * fx * fx * pb_yy

                     - 3.0 * pa_xxyy * fx * fz * fgb

                     - 12.0 * pa_xxy * fx * pb_y * fz * fgb

                     + 3.75 * fx * fx * fx * pa_yy * fz

                     + 30.0 * fx * fx * fx * pa_y * fz * pb_y

                     + 22.5 * fx * fx * fx * fz * pb_yy

                     + 9.0 * pa_xxyy * fz * fx * fx

                     + 72.0 * pa_xxy * fz * fx * fx * pb_y

                     + 54.0 * pa_xx * fx * fx * fz * pb_yy

                     - 1.5 * fx * fx * pb_yy * fz * fgb

                     - 3.0 * pa_xx * fx * pb_yy * fz * fgb

                     - 3.0 * pa_xx * fz * fga * pb_yy * fx

                     - 3.0 * fx * pa_yy * pb_yy * fz * fgb

                     - 3.0 * fz * fga * pa_yy * pb_yy * fx

                     - 4.0 * fz * fga * pa_y * fx * pb_yyy

                     - 6.0 * pa_xxyy * pb_yy * fz * fgb

                     + 18.0 * fx * fx * pa_yy * fz * pb_yy

                     + 24.0 * fx * fx * pa_y * fz * pb_yyy

                     + 42.0 * pa_xxyy * fz * pb_yy * fx

                     + 56.0 * pa_xxy * fz * fx * pb_yyy

                     - fx * fz * fga * pb_yyyy

                     - pa_xx * fz * fga * pb_yyyy

                     + 3.0 * fx * fx * fz * pb_yyyy

                     - fz * fga * pa_yy * pb_yyyy

                     + 7.0 * pa_xx * fz * fx * pb_yyyy

                     + 7.0 * fx * pa_yy * fz * pb_yyyy

                     + 16.0 * pa_xxyy * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_yyyz_s_0(double fx,
                                     double pa_xx,
                                     double pa_xxy,
                                     double pa_xxyy,
                                     double pa_y,
                                     double pa_yy,
                                     double pb_yyyz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * fx * pa_y * pb_z

                     + 1.125 * fx * fx * fx * pb_yz

                     + 1.5 * pa_xxy * fx * fx * pb_z

                     + 2.25 * pa_xx * fx * fx * pb_yz

                     + 0.75 * fx * fx * pa_yy * pb_yz

                     + 1.5 * fx * fx * pa_y * pb_yyz

                     + 1.5 * pa_xxyy * pb_yz * fx

                     + 3.0 * pa_xxy * fx * pb_yyz

                     + 0.25 * fx * fx * pb_yyyz

                     + 0.5 * pa_xx * fx * pb_yyyz

                     + 0.5 * fx * pa_yy * pb_yyyz

                     + pa_xxyy * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_yyyz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xx,
                                     double pa_xxy,
                                     double pa_xxyy,
                                     double pa_y,
                                     double pa_yy,
                                     double pb_yyyz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_y * fz * fgb * pb_z

                     - 1.5 * fz * fga * pa_y * fx * fx * pb_z

                     - 3.0 * fz * fga * fx * fx * pb_yz

                     - 3.0 * pa_xxy * fx * fz * fgb * pb_z

                     + 7.5 * fx * fx * fx * pa_y * fz * pb_z

                     + 11.25 * fx * fx * fx * fz * pb_yz

                     + 18.0 * pa_xxy * fz * fx * fx * pb_z

                     + 27.0 * pa_xx * fx * fx * fz * pb_yz

                     - 0.75 * fx * fx * pb_yz * fz * fgb

                     - 1.5 * pa_xx * fx * pb_yz * fz * fgb

                     - 1.5 * pa_xx * fz * fga * pb_yz * fx

                     - 1.5 * fx * pa_yy * pb_yz * fz * fgb

                     - 1.5 * fz * fga * pa_yy * pb_yz * fx

                     - 3.0 * fz * fga * pa_y * fx * pb_yyz

                     - 3.0 * pa_xxyy * pb_yz * fz * fgb

                     + 9.0 * fx * fx * pa_yy * fz * pb_yz

                     + 18.0 * fx * fx * pa_y * fz * pb_yyz

                     + 21.0 * pa_xxyy * fz * pb_yz * fx

                     + 42.0 * pa_xxy * fz * fx * pb_yyz

                     - fx * fz * fga * pb_yyyz

                     - pa_xx * fz * fga * pb_yyyz

                     + 3.0 * fx * fx * fz * pb_yyyz

                     - fz * fga * pa_yy * pb_yyyz

                     + 7.0 * pa_xx * fz * fx * pb_yyyz

                     + 7.0 * fx * pa_yy * fz * pb_yyyz

                     + 16.0 * pa_xxyy * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_yyzz_s_0(double fx,
                                     double pa_xx,
                                     double pa_xxy,
                                     double pa_xxyy,
                                     double pa_y,
                                     double pa_yy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyzz,
                                     double pb_yzz,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * pa_xx * fx * fx * fx

                     + 0.125 * fx * fx * fx * pa_yy

                     + 0.5 * fx * fx * fx * pa_y * pb_y

                     + 0.375 * fx * fx * fx * pb_zz

                     + 0.25 * pa_xxyy * fx * fx

                     + pa_xxy * fx * fx * pb_y

                     + 0.75 * pa_xx * fx * fx * pb_zz

                     + 0.125 * fx * fx * fx * pb_yy

                     + 0.25 * pa_xx * fx * fx * pb_yy

                     + 0.25 * fx * fx * pa_yy * pb_yy

                     + 0.25 * fx * fx * pa_yy * pb_zz

                     + fx * fx * pa_y * pb_yzz

                     + 0.5 * pa_xxyy * pb_yy * fx

                     + 0.5 * pa_xxyy * fx * pb_zz

                     + 2.0 * pa_xxy * fx * pb_yzz

                     + 0.25 * fx * fx * pb_yyzz

                     + 0.5 * pa_xx * fx * pb_yyzz

                     + 0.5 * fx * pa_yy * pb_yyzz

                     + pa_xxyy * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_yyzz_r_0(double fga,
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
                                     double pb_yyzz,
                                     double pb_yzz,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * fx * fz * fgb

                     - 0.5 * fz * fga * fx * fx * fx

                     - pa_xx * fx * fx * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     + 3.75 * pa_xx * fx * fx * fx * fz

                     - 0.25 * pa_xx * fz * fga * fx * fx

                     - 0.5 * fx * fx * pa_yy * fz * fgb

                     - fx * fx * pa_y * pb_y * fz * fgb

                     - 0.25 * fz * fga * pa_yy * fx * fx

                     - fz * fga * pa_y * fx * fx * pb_y

                     - fz * fga * fx * fx * pb_zz

                     - pa_xxyy * fx * fz * fgb

                     - 2.0 * pa_xxy * fx * pb_y * fz * fgb

                     + 1.25 * fx * fx * fx * pa_yy * fz

                     + 5.0 * fx * fx * fx * pa_y * fz * pb_y

                     + 3.75 * fx * fx * fx * fz * pb_zz

                     + 3.0 * pa_xxyy * fz * fx * fx

                     + 12.0 * pa_xxy * fz * fx * fx * pb_y

                     + 9.0 * pa_xx * fx * fx * fz * pb_zz

                     - 0.25 * fx * fx * pb_yy * fz * fgb

                     - 0.25 * fx * fx * fz * fgb * pb_zz

                     - 0.5 * fx * fx * fz * fga * pb_yy

                     - 0.5 * pa_xx * fx * pb_yy * fz * fgb

                     - 0.5 * pa_xx * fx * fz * fgb * pb_zz

                     - 0.5 * pa_xx * fz * fga * pb_yy * fx

                     - 0.5 * pa_xx * fz * fga * fx * pb_zz

                     - 0.5 * fx * pa_yy * pb_yy * fz * fgb

                     - 0.5 * fx * pa_yy * fz * fgb * pb_zz

                     + 1.25 * fx * fx * fx * fz * pb_yy

                     - 0.5 * fz * fga * pa_yy * pb_yy * fx

                     - 0.5 * fz * fga * pa_yy * fx * pb_zz

                     - 2.0 * fz * fga * pa_y * fx * pb_yzz

                     - pa_xxyy * pb_yy * fz * fgb

                     - pa_xxyy * fz * fgb * pb_zz

                     + 3.0 * pa_xx * fz * fx * fx * pb_yy

                     + 3.0 * fx * fx * pa_yy * fz * pb_yy

                     + 3.0 * fx * fx * pa_yy * fz * pb_zz

                     + 12.0 * fx * fx * pa_y * fz * pb_yzz

                     + 7.0 * pa_xxyy * fz * pb_yy * fx

                     + 7.0 * pa_xxyy * fz * fx * pb_zz

                     + 28.0 * pa_xxy * fz * fx * pb_yzz

                     - fx * fz * fga * pb_yyzz

                     - pa_xx * fz * fga * pb_yyzz

                     + 3.0 * fx * fx * fz * pb_yyzz

                     - fz * fga * pa_yy * pb_yyzz

                     + 7.0 * pa_xx * fz * fx * pb_yyzz

                     + 7.0 * fx * pa_yy * fz * pb_yyzz

                     + 16.0 * pa_xxyy * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_yzzz_s_0(double fx,
                                     double pa_xx,
                                     double pa_xxy,
                                     double pa_xxyy,
                                     double pa_y,
                                     double pa_yy,
                                     double pb_yz,
                                     double pb_yzzz,
                                     double pb_z,
                                     double pb_zzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * fx * pa_y * pb_z

                     + 1.5 * pa_xxy * fx * fx * pb_z

                     + 0.375 * fx * fx * fx * pb_yz

                     + 0.75 * pa_xx * fx * fx * pb_yz

                     + 0.75 * fx * fx * pa_yy * pb_yz

                     + 0.5 * fx * fx * pa_y * pb_zzz

                     + 1.5 * pa_xxyy * pb_yz * fx

                     + pa_xxy * fx * pb_zzz

                     + 0.25 * fx * fx * pb_yzzz

                     + 0.5 * pa_xx * fx * pb_yzzz

                     + 0.5 * fx * pa_yy * pb_yzzz

                     + pa_xxyy * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_yzzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xx,
                                     double pa_xxy,
                                     double pa_xxyy,
                                     double pa_y,
                                     double pa_yy,
                                     double pb_yz,
                                     double pb_yzzz,
                                     double pb_z,
                                     double pb_zzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_y * pb_z * fz * fgb

                     - 1.5 * fz * fga * pa_y * fx * fx * pb_z

                     - 3.0 * pa_xxy * fx * pb_z * fz * fgb

                     + 7.5 * fx * fx * fx * pa_y * fz * pb_z

                     + 18.0 * pa_xxy * fz * fx * fx * pb_z

                     - 0.75 * fx * fx * pb_yz * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_yz

                     - 1.5 * pa_xx * fx * pb_yz * fz * fgb

                     - 1.5 * pa_xx * fz * fga * pb_yz * fx

                     - 1.5 * fx * pa_yy * pb_yz * fz * fgb

                     + 3.75 * fx * fx * fx * fz * pb_yz

                     - 1.5 * fz * fga * pa_yy * pb_yz * fx

                     - fz * fga * pa_y * fx * pb_zzz

                     - 3.0 * pa_xxyy * pb_yz * fz * fgb

                     + 9.0 * pa_xx * fz * fx * fx * pb_yz

                     + 9.0 * fx * fx * pa_yy * fz * pb_yz

                     + 6.0 * fx * fx * pa_y * fz * pb_zzz

                     + 21.0 * pa_xxyy * fz * pb_yz * fx

                     + 14.0 * pa_xxy * fz * fx * pb_zzz

                     - fx * fz * fga * pb_yzzz

                     - pa_xx * fz * fga * pb_yzzz

                     + 3.0 * fx * fx * fz * pb_yzzz

                     - fz * fga * pa_yy * pb_yzzz

                     + 7.0 * pa_xx * fz * fx * pb_yzzz

                     + 7.0 * fx * pa_yy * fz * pb_yzzz

                     + 16.0 * pa_xxyy * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_zzzz_s_0(double fx,
                                     double pa_xx,
                                     double pa_xxyy,
                                     double pa_yy,
                                     double pb_zz,
                                     double pb_zzzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * pa_xx * fx * fx * fx

                     + 0.375 * fx * fx * fx * pa_yy

                     + 0.75 * pa_xxyy * fx * fx

                     + 0.75 * fx * fx * fx * pb_zz

                     + 1.5 * pa_xx * fx * fx * pb_zz

                     + 1.5 * fx * fx * pa_yy * pb_zz

                     + 3.0 * pa_xxyy * pb_zz * fx

                     + 0.25 * fx * fx * pb_zzzz

                     + 0.5 * pa_xx * fx * pb_zzzz

                     + 0.5 * fx * pa_yy * pb_zzzz

                     + pa_xxyy * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_zzzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xx,
                                     double pa_xxyy,
                                     double pa_yy,
                                     double pb_zz,
                                     double pb_zzzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fx * fz * fgb

                     - 0.75 * fx * fx * fx * fz * fga

                     - 1.5 * pa_xx * fx * fx * fz * fgb

                     - 0.75 * pa_xx * fz * fga * fx * fx

                     - 1.5 * fx * fx * pa_yy * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     - 0.75 * fz * fga * pa_yy * fx * fx

                     - 3.0 * pa_xxyy * fx * fz * fgb

                     + 3.75 * pa_xx * fz * fx * fx * fx

                     + 3.75 * fx * fx * fx * pa_yy * fz

                     + 9.0 * pa_xxyy * fz * fx * fx

                     - 1.5 * fx * fx * pb_zz * fz * fgb

                     - 3.0 * fx * fx * fz * fga * pb_zz

                     - 3.0 * pa_xx * fx * pb_zz * fz * fgb

                     - 3.0 * pa_xx * fz * fga * pb_zz * fx

                     - 3.0 * fx * pa_yy * pb_zz * fz * fgb

                     + 7.5 * fx * fx * fx * fz * pb_zz

                     - 3.0 * fz * fga * pa_yy * pb_zz * fx

                     - 6.0 * pa_xxyy * pb_zz * fz * fgb

                     + 18.0 * pa_xx * fz * fx * fx * pb_zz

                     + 18.0 * fx * fx * pa_yy * fz * pb_zz

                     + 42.0 * pa_xxyy * fz * pb_zz * fx

                     - fx * fz * fga * pb_zzzz

                     - pa_xx * fz * fga * pb_zzzz

                     + 3.0 * fx * fx * fz * pb_zzzz

                     - fz * fga * pa_yy * pb_zzzz

                     + 7.0 * pa_xx * fz * fx * pb_zzzz

                     + 7.0 * fx * pa_yy * fz * pb_zzzz

                     + 16.0 * pa_xxyy * fz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xxxx_s_0(double fx,
                                     double pa_xxyz,
                                     double pa_xyz,
                                     double pa_yz,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxx,
                                     double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pa_yz

                     + 0.75 * pa_xxyz * fx * fx

                     + 6.0 * pa_xyz * fx * fx * pb_x

                     + 4.5 * fx * fx * pa_yz * pb_xx

                     + 3.0 * pa_xxyz * pb_xx * fx

                     + 4.0 * pa_xyz * fx * pb_xxx

                     + 0.5 * fx * pa_yz * pb_xxxx

                     + pa_xxyz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xxxx_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xxyz,
                                     double pa_xyz,
                                     double pa_yz,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxx,
                                     double r_0_0)
    {
        return r_0_0 * (- 4.5 * fx * fx * pa_yz * fz * fgb

                     + 18.75 * fx * fx * fx * pa_yz * fz

                     - 0.75 * fz * fga * pa_yz * fx * fx

                     - 3.0 * pa_xxyz * fx * fz * fgb

                     - 12.0 * pa_xyz * fx * pb_x * fz * fgb

                     + 9.0 * pa_xxyz * fz * fx * fx

                     + 72.0 * pa_xyz * fx * fx * fz * pb_x

                     + 54.0 * fx * fx * pa_yz * fz * pb_xx

                     - 3.0 * fx * pa_yz * pb_xx * fz * fgb

                     - 3.0 * fz * fga * pa_yz * pb_xx * fx

                     - 6.0 * pa_xxyz * pb_xx * fz * fgb

                     + 42.0 * pa_xxyz * fz * pb_xx * fx

                     + 56.0 * pa_xyz * fx * fz * pb_xxx

                     - fz * fga * pa_yz * pb_xxxx

                     + 7.0 * fx * pa_yz * fz * pb_xxxx

                     + 16.0 * pa_xxyz * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xxxy_s_0(double fx,
                                     double pa_xxyz,
                                     double pa_xxz,
                                     double pa_xyz,
                                     double pa_xz,
                                     double pa_yz,
                                     double pa_z,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxy,
                                     double pb_xxy,
                                     double pb_xy,
                                     double pb_y,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xz * fx * fx * fx

                     + 1.125 * fx * fx * fx * pa_z * pb_x

                     + 0.75 * pa_xxz * fx * fx * pb_x

                     + 1.5 * pa_xyz * fx * fx * pb_y

                     + 1.5 * pa_xz * fx * fx * pb_xx

                     + 2.25 * fx * fx * pa_yz * pb_xy

                     + 0.25 * fx * fx * pa_z * pb_xxx

                     + 1.5 * pa_xxyz * pb_xy * fx

                     + 0.5 * pa_xxz * fx * pb_xxx

                     + 3.0 * pa_xyz * fx * pb_xxy

                     + 0.5 * fx * pa_yz * pb_xxxy

                     + pa_xxyz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xxxy_r_0(double fga,
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
                                     double pb_xxx,
                                     double pb_xxxy,
                                     double pb_xxy,
                                     double pb_xy,
                                     double pb_y,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xz * fx * fx * fz * fgb

                     + 7.5 * pa_xz * fx * fx * fx * fz

                     + 11.25 * fx * fx * fx * fz * pa_z * pb_x

                     - 0.75 * fx * fx * pa_z * pb_x * fz * fgb

                     - 0.75 * fz * fga * fx * fx * pa_z * pb_x

                     - 1.5 * pa_xxz * fx * pb_x * fz * fgb

                     - 3.0 * pa_xyz * fx * fz * fgb * pb_y

                     + 9.0 * pa_xxz * fx * fx * fz * pb_x

                     + 18.0 * pa_xyz * fx * fx * fz * pb_y

                     + 18.0 * pa_xz * fx * fx * fz * pb_xx

                     + 27.0 * fx * fx * pa_yz * fz * pb_xy

                     - 1.5 * fx * pa_yz * pb_xy * fz * fgb

                     - 1.5 * fz * fga * pa_yz * pb_xy * fx

                     - 0.5 * fz * fga * fx * pa_z * pb_xxx

                     - 3.0 * pa_xxyz * pb_xy * fz * fgb

                     + 3.0 * fx * fx * fz * pa_z * pb_xxx

                     + 21.0 * pa_xxyz * fz * pb_xy * fx

                     + 7.0 * pa_xxz * fx * fz * pb_xxx

                     + 42.0 * pa_xyz * fx * fz * pb_xxy

                     - fz * fga * pa_yz * pb_xxxy

                     + 7.0 * fx * pa_yz * fz * pb_xxxy

                     + 16.0 * pa_xxyz * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xxxz_s_0(double fx,
                                     double pa_xxy,
                                     double pa_xxyz,
                                     double pa_xy,
                                     double pa_xyz,
                                     double pa_y,
                                     double pa_yz,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxz,
                                     double pb_xxz,
                                     double pb_xz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xy * fx * fx * fx

                     + 1.125 * fx * fx * fx * pa_y * pb_x

                     + 0.75 * pa_xxy * fx * fx * pb_x

                     + 1.5 * pa_xyz * fx * fx * pb_z

                     + 1.5 * pa_xy * fx * fx * pb_xx

                     + 2.25 * fx * fx * pa_yz * pb_xz

                     + 0.25 * fx * fx * pa_y * pb_xxx

                     + 1.5 * pa_xxyz * pb_xz * fx

                     + 0.5 * pa_xxy * fx * pb_xxx

                     + 3.0 * pa_xyz * fx * pb_xxz

                     + 0.5 * fx * pa_yz * pb_xxxz

                     + pa_xxyz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xxxz_r_0(double fga,
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
                                     double pb_xxx,
                                     double pb_xxxz,
                                     double pb_xxz,
                                     double pb_xz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xy * fx * fx * fz * fgb

                     + 7.5 * pa_xy * fx * fx * fx * fz

                     + 11.25 * fx * fx * fx * pa_y * fz * pb_x

                     - 0.75 * fx * fx * pa_y * pb_x * fz * fgb

                     - 0.75 * fz * fga * pa_y * fx * fx * pb_x

                     - 1.5 * pa_xxy * fx * pb_x * fz * fgb

                     - 3.0 * pa_xyz * fx * fz * fgb * pb_z

                     + 9.0 * pa_xxy * fz * fx * fx * pb_x

                     + 18.0 * pa_xyz * fx * fx * fz * pb_z

                     + 18.0 * pa_xy * fx * fx * fz * pb_xx

                     + 27.0 * fx * fx * pa_yz * fz * pb_xz

                     - 1.5 * fx * pa_yz * pb_xz * fz * fgb

                     - 1.5 * fz * fga * pa_yz * pb_xz * fx

                     - 0.5 * fz * fga * pa_y * fx * pb_xxx

                     - 3.0 * pa_xxyz * pb_xz * fz * fgb

                     + 3.0 * fx * fx * pa_y * fz * pb_xxx

                     + 21.0 * pa_xxyz * fz * pb_xz * fx

                     + 7.0 * pa_xxy * fz * fx * pb_xxx

                     + 42.0 * pa_xyz * fx * fz * pb_xxz

                     - fz * fga * pa_yz * pb_xxxz

                     + 7.0 * fx * pa_yz * fz * pb_xxxz

                     + 16.0 * pa_xxyz * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xxyy_s_0(double fx,
                                     double pa_xxyz,
                                     double pa_xxz,
                                     double pa_xyz,
                                     double pa_xz,
                                     double pa_yz,
                                     double pa_z,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyy,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_y,
                                     double pb_yy,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_yz

                     + 0.75 * fx * fx * fx * pa_z * pb_y

                     + 0.25 * pa_xxyz * fx * fx

                     + 0.5 * pa_xxz * fx * fx * pb_y

                     + pa_xyz * fx * fx * pb_x

                     + 2.0 * pa_xz * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_yz * pb_yy

                     + 0.25 * fx * fx * pa_yz * pb_xx

                     + 0.5 * fx * fx * pa_z * pb_xxy

                     + 0.5 * pa_xxyz * pb_xx * fx

                     + 0.5 * pa_xxyz * fx * pb_yy

                     + pa_xxz * fx * pb_xxy

                     + 2.0 * pa_xyz * fx * pb_xyy

                     + 0.5 * fx * pa_yz * pb_xxyy

                     + pa_xxyz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xxyy_r_0(double fga,
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
                                     double pb_xxyy,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_y,
                                     double pb_yy,
                                     double r_0_0)
    {
        return r_0_0 * (- fx * fx * pa_yz * fz * fgb

                     + 3.75 * fx * fx * fx * pa_yz * fz

                     + 7.5 * fx * fx * fx * fz * pa_z * pb_y

                     - 0.5 * fx * fx * pa_z * fz * fgb * pb_y

                     - 0.25 * fz * fga * pa_yz * fx * fx

                     - 0.5 * fz * fga * fx * fx * pa_z * pb_y

                     - pa_xxyz * fx * fz * fgb

                     - pa_xxz * fx * fz * fgb * pb_y

                     - 2.0 * pa_xyz * fx * pb_x * fz * fgb

                     + 3.0 * pa_xxyz * fz * fx * fx

                     + 6.0 * pa_xxz * fx * fx * fz * pb_y

                     + 12.0 * pa_xyz * fx * fx * fz * pb_x

                     + 24.0 * pa_xz * fx * fx * fz * pb_xy

                     + 9.0 * fx * fx * pa_yz * fz * pb_yy

                     - 0.5 * fx * pa_yz * pb_xx * fz * fgb

                     - 0.5 * fx * pa_yz * fz * fgb * pb_yy

                     - 0.5 * fz * fga * pa_yz * pb_xx * fx

                     - 0.5 * fz * fga * pa_yz * fx * pb_yy

                     - fz * fga * fx * pa_z * pb_xxy

                     - pa_xxyz * pb_xx * fz * fgb

                     - pa_xxyz * fz * fgb * pb_yy

                     + 3.0 * fx * fx * pa_yz * fz * pb_xx

                     + 6.0 * fx * fx * fz * pa_z * pb_xxy

                     + 7.0 * pa_xxyz * fz * pb_xx * fx

                     + 7.0 * pa_xxyz * fz * fx * pb_yy

                     + 14.0 * pa_xxz * fx * fz * pb_xxy

                     + 28.0 * pa_xyz * fx * fz * pb_xyy

                     - fz * fga * pa_yz * pb_xxyy

                     + 7.0 * fx * pa_yz * fz * pb_xxyy

                     + 16.0 * pa_xxyz * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xxyz_s_0(double fx,
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
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyz,
                                     double pb_xxz,
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.125 * pa_xx * fx * fx * fx

                     + 0.5 * pa_x * fx * fx * fx * pb_x

                     + 0.375 * fx * fx * fx * pa_y * pb_y

                     + 0.375 * fx * fx * fx * pa_z * pb_z

                     + 0.125 * fx * fx * fx * pb_xx

                     + 0.25 * pa_xxy * fx * fx * pb_y

                     + 0.25 * pa_xxz * fx * fx * pb_z

                     + 0.25 * pa_xx * fx * fx * pb_xx

                     + pa_xy * fx * fx * pb_xy

                     + pa_xz * fx * fx * pb_xz

                     + 0.75 * fx * fx * pa_yz * pb_yz

                     + 0.25 * fx * fx * pa_y * pb_xxy

                     + 0.25 * fx * fx * pa_z * pb_xxz

                     + 0.5 * pa_xxyz * fx * pb_yz

                     + 0.5 * pa_xxy * fx * pb_xxy

                     + 0.5 * pa_xxz * fx * pb_xxz

                     + 2.0 * pa_xyz * fx * pb_xyz

                     + 0.5 * fx * pa_yz * pb_xxyz

                     + pa_xxyz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xxyz_r_0(double fga,
                                     double fgb,
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
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyz,
                                     double pb_xxz,
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (1.5 * fx * fx * fx * fx * fz

                     - 0.125 * fx * fx * fx * fz * fgb

                     - 0.125 * fz * fga * fx * fx * fx

                     - 0.25 * pa_xx * fx * fx * fz * fgb

                     + 1.25 * pa_xx * fx * fx * fx * fz

                     + 5.0 * pa_x * fx * fx * fx * fz * pb_x

                     + 3.75 * fx * fx * fx * pa_y * fz * pb_y

                     + 3.75 * fx * fx * fx * fz * pa_z * pb_z

                     - 0.25 * fx * fx * pa_y * fz * fgb * pb_y

                     - 0.25 * fx * fx * pa_z * fz * fgb * pb_z

                     - 0.25 * fz * fga * pa_y * fx * fx * pb_y

                     - 0.25 * fz * fga * fx * fx * pa_z * pb_z

                     - 0.25 * fz * fga * fx * fx * pb_xx

                     - 0.5 * pa_xxy * fx * fz * fgb * pb_y

                     - 0.5 * pa_xxz * fx * fz * fgb * pb_z

                     + 1.25 * fx * fx * fx * fz * pb_xx

                     + 3.0 * pa_xxy * fz * fx * fx * pb_y

                     + 3.0 * pa_xxz * fx * fx * fz * pb_z

                     + 3.0 * pa_xx * fx * fx * fz * pb_xx

                     + 12.0 * pa_xy * fx * fx * fz * pb_xy

                     + 12.0 * pa_xz * fx * fx * fz * pb_xz

                     + 9.0 * fx * fx * pa_yz * fz * pb_yz

                     - 0.5 * fx * pa_yz * fz * fgb * pb_yz

                     - 0.5 * fz * fga * pa_yz * fx * pb_yz

                     - 0.5 * fz * fga * pa_y * fx * pb_xxy

                     - 0.5 * fz * fga * fx * pa_z * pb_xxz

                     - pa_xxyz * fz * fgb * pb_yz

                     + 3.0 * fx * fx * pa_y * fz * pb_xxy

                     + 3.0 * fx * fx * fz * pa_z * pb_xxz

                     + 7.0 * pa_xxyz * fz * fx * pb_yz

                     + 7.0 * pa_xxy * fz * fx * pb_xxy

                     + 7.0 * pa_xxz * fx * fz * pb_xxz

                     + 28.0 * pa_xyz * fx * fz * pb_xyz

                     - fz * fga * pa_yz * pb_xxyz

                     + 7.0 * fx * pa_yz * fz * pb_xxyz

                     + 16.0 * pa_xxyz * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xxzz_s_0(double fx,
                                     double pa_xxy,
                                     double pa_xxyz,
                                     double pa_xy,
                                     double pa_xyz,
                                     double pa_y,
                                     double pa_yz,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxz,
                                     double pb_xxzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_yz

                     + 0.75 * fx * fx * fx * pa_y * pb_z

                     + 0.25 * pa_xxyz * fx * fx

                     + 0.5 * pa_xxy * fx * fx * pb_z

                     + pa_xyz * fx * fx * pb_x

                     + 2.0 * pa_xy * fx * fx * pb_xz

                     + 0.75 * fx * fx * pa_yz * pb_zz

                     + 0.25 * fx * fx * pa_yz * pb_xx

                     + 0.5 * fx * fx * pa_y * pb_xxz

                     + 0.5 * pa_xxyz * pb_xx * fx

                     + 0.5 * pa_xxyz * fx * pb_zz

                     + pa_xxy * fx * pb_xxz

                     + 2.0 * pa_xyz * fx * pb_xzz

                     + 0.5 * fx * pa_yz * pb_xxzz

                     + pa_xxyz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xxzz_r_0(double fga,
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
                                     double pb_xxzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- fx * fx * pa_yz * fz * fgb

                     + 3.75 * fx * fx * fx * pa_yz * fz

                     + 7.5 * fx * fx * fx * pa_y * fz * pb_z

                     - 0.5 * fx * fx * pa_y * fz * fgb * pb_z

                     - 0.25 * fz * fga * pa_yz * fx * fx

                     - 0.5 * fz * fga * pa_y * fx * fx * pb_z

                     - pa_xxyz * fx * fz * fgb

                     - pa_xxy * fx * fz * fgb * pb_z

                     - 2.0 * pa_xyz * fx * pb_x * fz * fgb

                     + 3.0 * pa_xxyz * fz * fx * fx

                     + 6.0 * pa_xxy * fz * fx * fx * pb_z

                     + 12.0 * pa_xyz * fx * fx * fz * pb_x

                     + 24.0 * pa_xy * fx * fx * fz * pb_xz

                     + 9.0 * fx * fx * pa_yz * fz * pb_zz

                     - 0.5 * fx * pa_yz * pb_xx * fz * fgb

                     - 0.5 * fx * pa_yz * fz * fgb * pb_zz

                     - 0.5 * fz * fga * pa_yz * pb_xx * fx

                     - 0.5 * fz * fga * pa_yz * fx * pb_zz

                     - fz * fga * pa_y * fx * pb_xxz

                     - pa_xxyz * pb_xx * fz * fgb

                     - pa_xxyz * fz * fgb * pb_zz

                     + 3.0 * fx * fx * pa_yz * fz * pb_xx

                     + 6.0 * fx * fx * pa_y * fz * pb_xxz

                     + 7.0 * pa_xxyz * fz * pb_xx * fx

                     + 7.0 * pa_xxyz * fz * fx * pb_zz

                     + 14.0 * pa_xxy * fz * fx * pb_xxz

                     + 28.0 * pa_xyz * fx * fz * pb_xzz

                     - fz * fga * pa_yz * pb_xxzz

                     + 7.0 * fx * pa_yz * fz * pb_xxzz

                     + 16.0 * pa_xxyz * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xyyy_s_0(double fx,
                                     double pa_xxyz,
                                     double pa_xxz,
                                     double pa_xyz,
                                     double pa_xz,
                                     double pa_yz,
                                     double pa_z,
                                     double pb_x,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_xyyy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xz * fx * fx * fx

                     + 0.375 * fx * fx * fx * pa_z * pb_x

                     + 0.75 * pa_xxz * fx * fx * pb_x

                     + 1.5 * pa_xyz * fx * fx * pb_y

                     + 1.5 * pa_xz * fx * fx * pb_yy

                     + 0.75 * fx * fx * pa_yz * pb_xy

                     + 0.75 * fx * fx * pa_z * pb_xyy

                     + 1.5 * pa_xxyz * pb_xy * fx

                     + 1.5 * pa_xxz * fx * pb_xyy

                     + pa_xyz * fx * pb_yyy

                     + 0.5 * fx * pa_yz * pb_xyyy

                     + pa_xxyz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xyyy_r_0(double fga,
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
                                     double pb_xyyy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xz * fx * fx * fz * fgb

                     + 7.5 * pa_xz * fx * fx * fx * fz

                     - 0.75 * fx * fx * pa_z * pb_x * fz * fgb

                     - 0.75 * fz * fga * fx * fx * pa_z * pb_x

                     - 1.5 * pa_xxz * fx * pb_x * fz * fgb

                     - 3.0 * pa_xyz * fx * pb_y * fz * fgb

                     + 3.75 * fx * fx * fx * fz * pa_z * pb_x

                     + 9.0 * pa_xxz * fx * fx * fz * pb_x

                     + 18.0 * pa_xyz * fx * fx * fz * pb_y

                     + 18.0 * pa_xz * fx * fx * fz * pb_yy

                     - 1.5 * fx * pa_yz * pb_xy * fz * fgb

                     - 1.5 * fz * fga * pa_yz * pb_xy * fx

                     - 1.5 * fz * fga * fx * pa_z * pb_xyy

                     - 3.0 * pa_xxyz * pb_xy * fz * fgb

                     + 9.0 * fx * fx * pa_yz * fz * pb_xy

                     + 9.0 * fx * fx * fz * pa_z * pb_xyy

                     + 21.0 * pa_xxyz * fz * pb_xy * fx

                     + 21.0 * pa_xxz * fx * fz * pb_xyy

                     + 14.0 * pa_xyz * fx * fz * pb_yyy

                     - fz * fga * pa_yz * pb_xyyy

                     + 7.0 * fx * pa_yz * fz * pb_xyyy

                     + 16.0 * pa_xxyz * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xyyz_s_0(double fx,
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
                                     double pb_xyy,
                                     double pb_xyyz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.25 * pa_xy * fx * fx * fx

                     + 0.5 * pa_x * fx * fx * fx * pb_y

                     + 0.125 * fx * fx * fx * pa_y * pb_x

                     + 0.25 * fx * fx * fx * pb_xy

                     + 0.25 * pa_xxy * fx * fx * pb_x

                     + 0.5 * pa_xx * fx * fx * pb_xy

                     + 0.5 * pa_xyz * fx * fx * pb_z

                     + 0.5 * pa_xy * fx * fx * pb_yy

                     + pa_xz * fx * fx * pb_yz

                     + 0.25 * fx * fx * pa_yz * pb_xz

                     + 0.25 * fx * fx * pa_y * pb_xyy

                     + 0.5 * fx * fx * pa_z * pb_xyz

                     + 0.5 * pa_xxyz * pb_xz * fx

                     + 0.5 * pa_xxy * fx * pb_xyy

                     + pa_xxz * fx * pb_xyz

                     + pa_xyz * fx * pb_yyz

                     + 0.5 * fx * pa_yz * pb_xyyz

                     + pa_xxyz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xyyz_r_0(double fga,
                                     double fgb,
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
                                     double pb_xyy,
                                     double pb_xyyz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_xy * fx * fx * fz * fgb

                     + 2.5 * pa_xy * fx * fx * fx * fz

                     + 5.0 * pa_x * fx * fx * fx * fz * pb_y

                     - 0.25 * fx * fx * pa_y * pb_x * fz * fgb

                     - 0.25 * fz * fga * pa_y * fx * fx * pb_x

                     - 0.5 * fz * fga * fx * fx * pb_xy

                     - 0.5 * pa_xxy * fx * pb_x * fz * fgb

                     - pa_xyz * fx * fz * fgb * pb_z

                     + 1.25 * fx * fx * fx * pa_y * fz * pb_x

                     + 2.5 * fx * fx * fx * fz * pb_xy

                     + 3.0 * pa_xxy * fz * fx * fx * pb_x

                     + 6.0 * pa_xx * fx * fx * fz * pb_xy

                     + 6.0 * pa_xyz * fx * fx * fz * pb_z

                     + 6.0 * pa_xy * fx * fx * fz * pb_yy

                     + 12.0 * pa_xz * fx * fx * fz * pb_yz

                     - 0.5 * fx * pa_yz * pb_xz * fz * fgb

                     - 0.5 * fz * fga * pa_yz * pb_xz * fx

                     - 0.5 * fz * fga * pa_y * fx * pb_xyy

                     - fz * fga * fx * pa_z * pb_xyz

                     - pa_xxyz * pb_xz * fz * fgb

                     + 3.0 * fx * fx * pa_yz * fz * pb_xz

                     + 3.0 * fx * fx * pa_y * fz * pb_xyy

                     + 6.0 * fx * fx * fz * pa_z * pb_xyz

                     + 7.0 * pa_xxyz * fz * pb_xz * fx

                     + 7.0 * pa_xxy * fz * fx * pb_xyy

                     + 14.0 * pa_xxz * fx * fz * pb_xyz

                     + 14.0 * pa_xyz * fx * fz * pb_yyz

                     - fz * fga * pa_yz * pb_xyyz

                     + 7.0 * fx * pa_yz * fz * pb_xyyz

                     + 16.0 * pa_xxyz * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xyzz_s_0(double fx,
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
                                     double pb_xyzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.25 * pa_xz * fx * fx * fx

                     + 0.5 * pa_x * fx * fx * fx * pb_z

                     + 0.125 * fx * fx * fx * pa_z * pb_x

                     + 0.25 * fx * fx * fx * pb_xz

                     + 0.25 * pa_xxz * fx * fx * pb_x

                     + 0.5 * pa_xx * fx * fx * pb_xz

                     + 0.5 * pa_xyz * fx * fx * pb_y

                     + pa_xy * fx * fx * pb_yz

                     + 0.5 * pa_xz * fx * fx * pb_zz

                     + 0.25 * fx * fx * pa_yz * pb_xy

                     + 0.5 * fx * fx * pa_y * pb_xyz

                     + 0.25 * fx * fx * pa_z * pb_xzz

                     + 0.5 * pa_xxyz * pb_xy * fx

                     + pa_xxy * fx * pb_xyz

                     + 0.5 * pa_xxz * fx * pb_xzz

                     + pa_xyz * fx * pb_yzz

                     + 0.5 * fx * pa_yz * pb_xyzz

                     + pa_xxyz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xyzz_r_0(double fga,
                                     double fgb,
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
                                     double pb_xyzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_xz * fx * fx * fz * fgb

                     + 2.5 * pa_xz * fx * fx * fx * fz

                     + 5.0 * pa_x * fx * fx * fx * fz * pb_z

                     - 0.25 * fx * fx * pa_z * pb_x * fz * fgb

                     - 0.25 * fz * fga * fx * fx * pa_z * pb_x

                     - 0.5 * fz * fga * fx * fx * pb_xz

                     - 0.5 * pa_xxz * fx * pb_x * fz * fgb

                     - pa_xyz * fx * pb_y * fz * fgb

                     + 1.25 * fx * fx * fx * fz * pa_z * pb_x

                     + 2.5 * fx * fx * fx * fz * pb_xz

                     + 3.0 * pa_xxz * fx * fx * fz * pb_x

                     + 6.0 * pa_xx * fx * fx * fz * pb_xz

                     + 6.0 * pa_xyz * fx * fx * fz * pb_y

                     + 12.0 * pa_xy * fx * fx * fz * pb_yz

                     + 6.0 * pa_xz * fx * fx * fz * pb_zz

                     - 0.5 * fx * pa_yz * pb_xy * fz * fgb

                     - 0.5 * fz * fga * pa_yz * pb_xy * fx

                     - fz * fga * pa_y * fx * pb_xyz

                     - 0.5 * fz * fga * fx * pa_z * pb_xzz

                     - pa_xxyz * pb_xy * fz * fgb

                     + 3.0 * fx * fx * pa_yz * fz * pb_xy

                     + 6.0 * fx * fx * pa_y * fz * pb_xyz

                     + 3.0 * fx * fx * fz * pa_z * pb_xzz

                     + 7.0 * pa_xxyz * fz * pb_xy * fx

                     + 14.0 * pa_xxy * fz * fx * pb_xyz

                     + 7.0 * pa_xxz * fx * fz * pb_xzz

                     + 14.0 * pa_xyz * fx * fz * pb_yzz

                     - fz * fga * pa_yz * pb_xyzz

                     + 7.0 * fx * pa_yz * fz * pb_xyzz

                     + 16.0 * pa_xxyz * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xzzz_s_0(double fx,
                                     double pa_xxy,
                                     double pa_xxyz,
                                     double pa_xy,
                                     double pa_xyz,
                                     double pa_y,
                                     double pa_yz,
                                     double pb_x,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_xzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xy * fx * fx * fx

                     + 0.375 * fx * fx * fx * pa_y * pb_x

                     + 0.75 * pa_xxy * fx * fx * pb_x

                     + 1.5 * pa_xyz * fx * fx * pb_z

                     + 1.5 * pa_xy * fx * fx * pb_zz

                     + 0.75 * fx * fx * pa_yz * pb_xz

                     + 0.75 * fx * fx * pa_y * pb_xzz

                     + 1.5 * pa_xxyz * pb_xz * fx

                     + 1.5 * pa_xxy * fx * pb_xzz

                     + pa_xyz * fx * pb_zzz

                     + 0.5 * fx * pa_yz * pb_xzzz

                     + pa_xxyz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xzzz_r_0(double fga,
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
                                     double pb_xzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xy * fx * fx * fz * fgb

                     + 7.5 * pa_xy * fx * fx * fx * fz

                     - 0.75 * fx * fx * pa_y * pb_x * fz * fgb

                     - 0.75 * fz * fga * pa_y * fx * fx * pb_x

                     - 1.5 * pa_xxy * fx * pb_x * fz * fgb

                     - 3.0 * pa_xyz * fx * pb_z * fz * fgb

                     + 3.75 * fx * fx * fx * pa_y * fz * pb_x

                     + 9.0 * pa_xxy * fz * fx * fx * pb_x

                     + 18.0 * pa_xyz * fx * fx * fz * pb_z

                     + 18.0 * pa_xy * fx * fx * fz * pb_zz

                     - 1.5 * fx * pa_yz * pb_xz * fz * fgb

                     - 1.5 * fz * fga * pa_yz * pb_xz * fx

                     - 1.5 * fz * fga * pa_y * fx * pb_xzz

                     - 3.0 * pa_xxyz * pb_xz * fz * fgb

                     + 9.0 * fx * fx * pa_yz * fz * pb_xz

                     + 9.0 * fx * fx * pa_y * fz * pb_xzz

                     + 21.0 * pa_xxyz * fz * pb_xz * fx

                     + 21.0 * pa_xxy * fz * fx * pb_xzz

                     + 14.0 * pa_xyz * fx * fz * pb_zzz

                     - fz * fga * pa_yz * pb_xzzz

                     + 7.0 * fx * pa_yz * fz * pb_xzzz

                     + 16.0 * pa_xxyz * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_yyyy_s_0(double fx,
                                     double pa_xxyz,
                                     double pa_xxz,
                                     double pa_yz,
                                     double pa_z,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double pb_yyyy,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_yz

                     + 1.5 * fx * fx * fx * pa_z * pb_y

                     + 0.75 * pa_xxyz * fx * fx

                     + 3.0 * pa_xxz * fx * fx * pb_y

                     + 1.5 * fx * fx * pa_yz * pb_yy

                     + fx * fx * pa_z * pb_yyy

                     + 3.0 * pa_xxyz * pb_yy * fx

                     + 2.0 * pa_xxz * fx * pb_yyy

                     + 0.5 * fx * pa_yz * pb_yyyy

                     + pa_xxyz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_yyyy_r_0(double fga,
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
                                     double pb_yyyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_yz * fz * fgb

                     - 3.0 * fx * fx * pa_z * pb_y * fz * fgb

                     - 0.75 * fz * fga * pa_yz * fx * fx

                     - 3.0 * fz * fga * fx * fx * pa_z * pb_y

                     - 3.0 * pa_xxyz * fx * fz * fgb

                     - 6.0 * pa_xxz * fx * pb_y * fz * fgb

                     + 3.75 * fx * fx * fx * pa_yz * fz

                     + 15.0 * fx * fx * fx * fz * pa_z * pb_y

                     + 9.0 * pa_xxyz * fz * fx * fx

                     + 36.0 * pa_xxz * fx * fx * fz * pb_y

                     - 3.0 * fx * pa_yz * pb_yy * fz * fgb

                     - 3.0 * fz * fga * pa_yz * pb_yy * fx

                     - 2.0 * fz * fga * fx * pa_z * pb_yyy

                     - 6.0 * pa_xxyz * pb_yy * fz * fgb

                     + 18.0 * fx * fx * pa_yz * fz * pb_yy

                     + 12.0 * fx * fx * fz * pa_z * pb_yyy

                     + 42.0 * pa_xxyz * fz * pb_yy * fx

                     + 28.0 * pa_xxz * fx * fz * pb_yyy

                     - fz * fga * pa_yz * pb_yyyy

                     + 7.0 * fx * pa_yz * fz * pb_yyyy

                     + 16.0 * pa_xxyz * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_yyyz_s_0(double fx,
                                     double pa_xx,
                                     double pa_xxy,
                                     double pa_xxyz,
                                     double pa_xxz,
                                     double pa_y,
                                     double pa_yz,
                                     double pa_z,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double pb_yyyz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * pa_xx * fx * fx * fx

                     + 0.375 * fx * fx * fx * pa_y * pb_y

                     + 0.375 * fx * fx * fx * pa_z * pb_z

                     + 0.375 * fx * fx * fx * pb_yy

                     + 0.75 * pa_xxy * fx * fx * pb_y

                     + 0.75 * pa_xxz * fx * fx * pb_z

                     + 0.75 * pa_xx * fx * fx * pb_yy

                     + 0.75 * fx * fx * pa_yz * pb_yz

                     + 0.25 * fx * fx * pa_y * pb_yyy

                     + 0.75 * fx * fx * pa_z * pb_yyz

                     + 1.5 * pa_xxyz * pb_yz * fx

                     + 0.5 * pa_xxy * fx * pb_yyy

                     + 1.5 * pa_xxz * fx * pb_yyz

                     + 0.5 * fx * pa_yz * pb_yyyz

                     + pa_xxyz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_yyyz_r_0(double fga,
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
                                     double pb_yyy,
                                     double pb_yyyz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.375 * fx * fx * fx * fz * fgb

                     - 0.375 * fz * fga * fx * fx * fx

                     - 0.75 * pa_xx * fx * fx * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     + 3.75 * pa_xx * fx * fx * fx * fz

                     - 0.75 * fx * fx * pa_y * pb_y * fz * fgb

                     - 0.75 * fx * fx * pa_z * fz * fgb * pb_z

                     - 0.75 * fz * fga * pa_y * fx * fx * pb_y

                     - 0.75 * fz * fga * fx * fx * pa_z * pb_z

                     - 0.75 * fz * fga * fx * fx * pb_yy

                     - 1.5 * pa_xxy * fx * pb_y * fz * fgb

                     - 1.5 * pa_xxz * fx * fz * fgb * pb_z

                     + 3.75 * fx * fx * fx * pa_y * fz * pb_y

                     + 3.75 * fx * fx * fx * fz * pa_z * pb_z

                     + 3.75 * fx * fx * fx * fz * pb_yy

                     + 9.0 * pa_xxy * fz * fx * fx * pb_y

                     + 9.0 * pa_xxz * fx * fx * fz * pb_z

                     + 9.0 * pa_xx * fx * fx * fz * pb_yy

                     - 1.5 * fx * pa_yz * pb_yz * fz * fgb

                     - 1.5 * fz * fga * pa_yz * pb_yz * fx

                     - 0.5 * fz * fga * pa_y * fx * pb_yyy

                     - 1.5 * fz * fga * fx * pa_z * pb_yyz

                     - 3.0 * pa_xxyz * pb_yz * fz * fgb

                     + 9.0 * fx * fx * pa_yz * fz * pb_yz

                     + 3.0 * fx * fx * pa_y * fz * pb_yyy

                     + 9.0 * fx * fx * fz * pa_z * pb_yyz

                     + 21.0 * pa_xxyz * fz * pb_yz * fx

                     + 7.0 * pa_xxy * fz * fx * pb_yyy

                     + 21.0 * pa_xxz * fx * fz * pb_yyz

                     - fz * fga * pa_yz * pb_yyyz

                     + 7.0 * fx * pa_yz * fz * pb_yyyz

                     + 16.0 * pa_xxyz * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_yyzz_s_0(double fx,
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
                                     double pb_yyzz,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx * pa_yz

                     + 0.25 * fx * fx * fx * pa_y * pb_z

                     + 0.25 * fx * fx * fx * pa_z * pb_y

                     + 0.5 * fx * fx * fx * pb_yz

                     + 0.25 * pa_xxyz * fx * fx

                     + 0.5 * pa_xxy * fx * fx * pb_z

                     + 0.5 * pa_xxz * fx * fx * pb_y

                     + pa_xx * fx * fx * pb_yz

                     + 0.25 * fx * fx * pa_yz * pb_yy

                     + 0.25 * fx * fx * pa_yz * pb_zz

                     + 0.5 * fx * fx * pa_y * pb_yyz

                     + 0.5 * fx * fx * pa_z * pb_yzz

                     + 0.5 * pa_xxyz * pb_yy * fx

                     + 0.5 * pa_xxyz * fx * pb_zz

                     + pa_xxy * fx * pb_yyz

                     + pa_xxz * fx * pb_yzz

                     + 0.5 * fx * pa_yz * pb_yyzz

                     + pa_xxyz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_yyzz_r_0(double fga,
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
                                     double pb_yyzz,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * pa_yz * fz * fgb

                     - 0.5 * fx * fx * pa_y * fz * fgb * pb_z

                     - 0.5 * fx * fx * pa_z * pb_y * fz * fgb

                     - 0.25 * fz * fga * pa_yz * fx * fx

                     - 0.5 * fz * fga * pa_y * fx * fx * pb_z

                     - 0.5 * fz * fga * fx * fx * pa_z * pb_y

                     - fz * fga * fx * fx * pb_yz

                     - pa_xxyz * fx * fz * fgb

                     - pa_xxy * fx * fz * fgb * pb_z

                     - pa_xxz * fx * pb_y * fz * fgb

                     + 1.25 * fx * fx * fx * pa_yz * fz

                     + 2.5 * fx * fx * fx * pa_y * fz * pb_z

                     + 2.5 * fx * fx * fx * fz * pa_z * pb_y

                     + 5.0 * fx * fx * fx * fz * pb_yz

                     + 3.0 * pa_xxyz * fz * fx * fx

                     + 6.0 * pa_xxy * fz * fx * fx * pb_z

                     + 6.0 * pa_xxz * fx * fx * fz * pb_y

                     + 12.0 * pa_xx * fx * fx * fz * pb_yz

                     - 0.5 * fx * pa_yz * pb_yy * fz * fgb

                     - 0.5 * fx * pa_yz * fz * fgb * pb_zz

                     - 0.5 * fz * fga * pa_yz * pb_yy * fx

                     - 0.5 * fz * fga * pa_yz * fx * pb_zz

                     - fz * fga * pa_y * fx * pb_yyz

                     - fz * fga * fx * pa_z * pb_yzz

                     - pa_xxyz * pb_yy * fz * fgb

                     - pa_xxyz * fz * fgb * pb_zz

                     + 3.0 * fx * fx * pa_yz * fz * pb_yy

                     + 3.0 * fx * fx * pa_yz * fz * pb_zz

                     + 6.0 * fx * fx * pa_y * fz * pb_yyz

                     + 6.0 * fx * fx * fz * pa_z * pb_yzz

                     + 7.0 * pa_xxyz * fz * pb_yy * fx

                     + 7.0 * pa_xxyz * fz * fx * pb_zz

                     + 14.0 * pa_xxy * fz * fx * pb_yyz

                     + 14.0 * pa_xxz * fx * fz * pb_yzz

                     - fz * fga * pa_yz * pb_yyzz

                     + 7.0 * fx * pa_yz * fz * pb_yyzz

                     + 16.0 * pa_xxyz * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_yzzz_s_0(double fx,
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
                                     double pb_yzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * pa_xx * fx * fx * fx

                     + 0.375 * fx * fx * fx * pa_y * pb_y

                     + 0.375 * fx * fx * fx * pa_z * pb_z

                     + 0.375 * fx * fx * fx * pb_zz

                     + 0.75 * pa_xxy * fx * fx * pb_y

                     + 0.75 * pa_xxz * fx * fx * pb_z

                     + 0.75 * pa_xx * fx * fx * pb_zz

                     + 0.75 * fx * fx * pa_yz * pb_yz

                     + 0.75 * fx * fx * pa_y * pb_yzz

                     + 0.25 * fx * fx * pa_z * pb_zzz

                     + 1.5 * pa_xxyz * pb_yz * fx

                     + 1.5 * pa_xxy * fx * pb_yzz

                     + 0.5 * pa_xxz * fx * pb_zzz

                     + 0.5 * fx * pa_yz * pb_yzzz

                     + pa_xxyz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_yzzz_r_0(double fga,
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
                                     double pb_yzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.375 * fx * fx * fx * fz * fgb

                     - 0.375 * fz * fga * fx * fx * fx

                     - 0.75 * pa_xx * fx * fx * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     + 3.75 * pa_xx * fx * fx * fx * fz

                     - 0.75 * fx * fx * pa_y * pb_y * fz * fgb

                     - 0.75 * fx * fx * pa_z * pb_z * fz * fgb

                     - 0.75 * fz * fga * pa_y * fx * fx * pb_y

                     - 0.75 * fz * fga * fx * fx * pa_z * pb_z

                     - 0.75 * fz * fga * fx * fx * pb_zz

                     - 1.5 * pa_xxy * fx * pb_y * fz * fgb

                     - 1.5 * pa_xxz * fx * pb_z * fz * fgb

                     + 3.75 * fx * fx * fx * pa_y * fz * pb_y

                     + 3.75 * fx * fx * fx * fz * pa_z * pb_z

                     + 3.75 * fx * fx * fx * fz * pb_zz

                     + 9.0 * pa_xxy * fz * fx * fx * pb_y

                     + 9.0 * pa_xxz * fx * fx * fz * pb_z

                     + 9.0 * pa_xx * fx * fx * fz * pb_zz

                     - 1.5 * fx * pa_yz * pb_yz * fz * fgb

                     - 1.5 * fz * fga * pa_yz * pb_yz * fx

                     - 1.5 * fz * fga * pa_y * fx * pb_yzz

                     - 0.5 * fz * fga * fx * pa_z * pb_zzz

                     - 3.0 * pa_xxyz * pb_yz * fz * fgb

                     + 9.0 * fx * fx * pa_yz * fz * pb_yz

                     + 9.0 * fx * fx * pa_y * fz * pb_yzz

                     + 3.0 * fx * fx * fz * pa_z * pb_zzz

                     + 21.0 * pa_xxyz * fz * pb_yz * fx

                     + 21.0 * pa_xxy * fz * fx * pb_yzz

                     + 7.0 * pa_xxz * fx * fz * pb_zzz

                     - fz * fga * pa_yz * pb_yzzz

                     + 7.0 * fx * pa_yz * fz * pb_yzzz

                     + 16.0 * pa_xxyz * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_zzzz_s_0(double fx,
                                     double pa_xxy,
                                     double pa_xxyz,
                                     double pa_y,
                                     double pa_yz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double pb_zzzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_yz

                     + 1.5 * fx * fx * fx * pa_y * pb_z

                     + 0.75 * pa_xxyz * fx * fx

                     + 3.0 * pa_xxy * fx * fx * pb_z

                     + 1.5 * fx * fx * pa_yz * pb_zz

                     + fx * fx * pa_y * pb_zzz

                     + 3.0 * pa_xxyz * pb_zz * fx

                     + 2.0 * pa_xxy * fx * pb_zzz

                     + 0.5 * fx * pa_yz * pb_zzzz

                     + pa_xxyz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_zzzz_r_0(double fga,
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
                                     double pb_zzzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_yz * fz * fgb

                     - 3.0 * fx * fx * pa_y * pb_z * fz * fgb

                     - 0.75 * fz * fga * pa_yz * fx * fx

                     - 3.0 * fz * fga * pa_y * fx * fx * pb_z

                     - 3.0 * pa_xxyz * fx * fz * fgb

                     - 6.0 * pa_xxy * fx * pb_z * fz * fgb

                     + 3.75 * fx * fx * fx * pa_yz * fz

                     + 15.0 * fx * fx * fx * pa_y * fz * pb_z

                     + 9.0 * pa_xxyz * fz * fx * fx

                     + 36.0 * pa_xxy * fz * fx * fx * pb_z

                     - 3.0 * fx * pa_yz * pb_zz * fz * fgb

                     - 3.0 * fz * fga * pa_yz * pb_zz * fx

                     - 2.0 * fz * fga * pa_y * fx * pb_zzz

                     - 6.0 * pa_xxyz * pb_zz * fz * fgb

                     + 18.0 * fx * fx * pa_yz * fz * pb_zz

                     + 12.0 * fx * fx * pa_y * fz * pb_zzz

                     + 42.0 * pa_xxyz * fz * pb_zz * fx

                     + 28.0 * pa_xxy * fz * fx * pb_zzz

                     - fz * fga * pa_yz * pb_zzzz

                     + 7.0 * fx * pa_yz * fz * pb_zzzz

                     + 16.0 * pa_xxyz * fz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xxxx_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxzz,
                                     double pa_xzz,
                                     double pa_zz,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxx,
                                     double s_0_0)
    {
        return s_0_0 * (0.9375 * fx * fx * fx * fx

                     + 1.875 * fx * fx * fx * pa_zz

                     + 0.375 * pa_xx * fx * fx * fx

                     + 3.0 * pa_x * fx * fx * fx * pb_x

                     + 2.25 * fx * fx * fx * pb_xx

                     + 0.75 * pa_xxzz * fx * fx

                     + 6.0 * pa_xzz * fx * fx * pb_x

                     + 4.5 * fx * fx * pa_zz * pb_xx

                     + 1.5 * pa_xx * fx * fx * pb_xx

                     + 2.0 * pa_x * fx * fx * pb_xxx

                     + 3.0 * pa_xxzz * pb_xx * fx

                     + 4.0 * pa_xzz * fx * pb_xxx

                     + 0.25 * fx * fx * pb_xxxx

                     + 0.5 * pa_xx * fx * pb_xxxx

                     + 0.5 * fx * pa_zz * pb_xxxx

                     + pa_xxzz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xxxx_r_0(double fga,
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
                                     double pb_xxxx,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * fx * fz * fgb

                     - 2.25 * fx * fx * fx * fz * fga

                     - 4.5 * fx * fx * pa_zz * fz * fgb

                     + 7.5 * fx * fx * fx * fx * fz

                     + 18.75 * fx * fx * fx * pa_zz * fz

                     - 1.5 * pa_xx * fx * fx * fz * fgb

                     - 0.75 * pa_xx * fz * fga * fx * fx

                     - 6.0 * pa_x * fx * fx * pb_x * fz * fgb

                     - 6.0 * pa_x * fx * fx * fz * fga * pb_x

                     - 6.0 * fx * fx * fz * fga * pb_xx

                     - 0.75 * fz * fga * pa_zz * fx * fx

                     - 3.0 * pa_xxzz * fx * fz * fgb

                     + 3.75 * pa_xx * fz * fx * fx * fx

                     - 12.0 * pa_xzz * fx * pb_x * fz * fgb

                     + 30.0 * pa_x * fx * fx * fx * fz * pb_x

                     + 22.5 * fx * fx * fx * fz * pb_xx

                     + 9.0 * pa_xxzz * fz * fx * fx

                     + 72.0 * pa_xzz * fx * fx * fz * pb_x

                     + 54.0 * fx * fx * pa_zz * fz * pb_xx

                     - 1.5 * fx * fx * pb_xx * fz * fgb

                     - 3.0 * pa_xx * fx * pb_xx * fz * fgb

                     - 3.0 * pa_xx * fz * fga * pb_xx * fx

                     - 4.0 * pa_x * fx * fz * fga * pb_xxx

                     - 3.0 * fx * pa_zz * pb_xx * fz * fgb

                     - 3.0 * fz * fga * pa_zz * pb_xx * fx

                     - 6.0 * pa_xxzz * pb_xx * fz * fgb

                     + 18.0 * pa_xx * fz * fx * fx * pb_xx

                     + 24.0 * pa_x * fx * fx * fz * pb_xxx

                     + 42.0 * pa_xxzz * fz * pb_xx * fx

                     + 56.0 * pa_xzz * fx * fz * pb_xxx

                     - fx * fz * fga * pb_xxxx

                     - pa_xx * fz * fga * pb_xxxx

                     + 3.0 * fx * fx * fz * pb_xxxx

                     - fz * fga * pa_zz * pb_xxxx

                     + 7.0 * pa_xx * fz * fx * pb_xxxx

                     + 7.0 * fx * pa_zz * fz * pb_xxxx

                     + 16.0 * pa_xxzz * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xxxy_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxzz,
                                     double pa_xzz,
                                     double pa_zz,
                                     double pb_xxxy,
                                     double pb_xxy,
                                     double pb_xy,
                                     double pb_y,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * fx * pb_y

                     + 1.125 * fx * fx * fx * pb_xy

                     + 1.5 * pa_xzz * fx * fx * pb_y

                     + 2.25 * fx * fx * pa_zz * pb_xy

                     + 0.75 * pa_xx * fx * fx * pb_xy

                     + 1.5 * pa_x * fx * fx * pb_xxy

                     + 1.5 * pa_xxzz * pb_xy * fx

                     + 3.0 * pa_xzz * fx * pb_xxy

                     + 0.25 * fx * fx * pb_xxxy

                     + 0.5 * pa_xx * fx * pb_xxxy

                     + 0.5 * fx * pa_zz * pb_xxxy

                     + pa_xxzz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xxxy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxzz,
                                     double pa_xzz,
                                     double pa_zz,
                                     double pb_xxxy,
                                     double pb_xxy,
                                     double pb_xy,
                                     double pb_y,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fx * fz * fgb * pb_y

                     - 1.5 * pa_x * fx * fx * fz * fga * pb_y

                     - 3.0 * fx * fx * fz * fga * pb_xy

                     - 3.0 * pa_xzz * fx * fz * fgb * pb_y

                     + 7.5 * pa_x * fx * fx * fx * fz * pb_y

                     + 11.25 * fx * fx * fx * fz * pb_xy

                     + 18.0 * pa_xzz * fx * fx * fz * pb_y

                     + 27.0 * fx * fx * pa_zz * fz * pb_xy

                     - 0.75 * fx * fx * pb_xy * fz * fgb

                     - 1.5 * pa_xx * fx * pb_xy * fz * fgb

                     - 1.5 * pa_xx * fz * fga * pb_xy * fx

                     - 3.0 * pa_x * fx * fz * fga * pb_xxy

                     - 1.5 * fx * pa_zz * pb_xy * fz * fgb

                     - 1.5 * fz * fga * pa_zz * pb_xy * fx

                     - 3.0 * pa_xxzz * pb_xy * fz * fgb

                     + 9.0 * pa_xx * fz * fx * fx * pb_xy

                     + 18.0 * pa_x * fx * fx * fz * pb_xxy

                     + 21.0 * pa_xxzz * fz * pb_xy * fx

                     + 42.0 * pa_xzz * fx * fz * pb_xxy

                     - fx * fz * fga * pb_xxxy

                     - pa_xx * fz * fga * pb_xxxy

                     + 3.0 * fx * fx * fz * pb_xxxy

                     - fz * fga * pa_zz * pb_xxxy

                     + 7.0 * pa_xx * fz * fx * pb_xxxy

                     + 7.0 * fx * pa_zz * fz * pb_xxxy

                     + 16.0 * pa_xxzz * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xxxz_s_0(double fx,
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
                                     double pb_xxx,
                                     double pb_xxxz,
                                     double pb_xxz,
                                     double pb_xz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xz * fx * fx * fx

                     + 2.25 * fx * fx * fx * pa_z * pb_x

                     + 0.75 * pa_x * fx * fx * fx * pb_z

                     + 1.125 * fx * fx * fx * pb_xz

                     + 1.5 * pa_xxz * fx * fx * pb_x

                     + 1.5 * pa_xzz * fx * fx * pb_z

                     + 3.0 * pa_xz * fx * fx * pb_xx

                     + 2.25 * fx * fx * pa_zz * pb_xz

                     + 0.75 * pa_xx * fx * fx * pb_xz

                     + 1.5 * pa_x * fx * fx * pb_xxz

                     + 0.5 * fx * fx * pa_z * pb_xxx

                     + 1.5 * pa_xxzz * pb_xz * fx

                     + pa_xxz * fx * pb_xxx

                     + 3.0 * pa_xzz * fx * pb_xxz

                     + 0.25 * fx * fx * pb_xxxz

                     + 0.5 * pa_xx * fx * pb_xxxz

                     + 0.5 * fx * pa_zz * pb_xxxz

                     + pa_xxzz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xxxz_r_0(double fga,
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
                                     double pb_xxx,
                                     double pb_xxxz,
                                     double pb_xxz,
                                     double pb_xz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xz * fx * fx * fz * fgb

                     + 15.0 * pa_xz * fx * fx * fx * fz

                     + 22.5 * fx * fx * fx * pa_z * fz * pb_x

                     - 1.5 * pa_x * fx * fx * fz * fgb * pb_z

                     - 1.5 * pa_x * fx * fx * fz * fga * pb_z

                     - 1.5 * fx * fx * pa_z * pb_x * fz * fgb

                     - 3.0 * fx * fx * fz * fga * pb_xz

                     - 1.5 * fz * fga * pa_z * fx * fx * pb_x

                     - 3.0 * pa_xxz * fx * pb_x * fz * fgb

                     - 3.0 * pa_xzz * fx * fz * fgb * pb_z

                     + 7.5 * pa_x * fx * fx * fx * fz * pb_z

                     + 11.25 * fx * fx * fx * fz * pb_xz

                     + 18.0 * pa_xxz * fz * fx * fx * pb_x

                     + 18.0 * pa_xzz * fx * fx * fz * pb_z

                     + 36.0 * pa_xz * fx * fx * fz * pb_xx

                     + 27.0 * fx * fx * pa_zz * fz * pb_xz

                     - 0.75 * fx * fx * pb_xz * fz * fgb

                     - 1.5 * pa_xx * fx * pb_xz * fz * fgb

                     - 1.5 * pa_xx * fz * fga * pb_xz * fx

                     - 3.0 * pa_x * fx * fz * fga * pb_xxz

                     - 1.5 * fx * pa_zz * pb_xz * fz * fgb

                     - 1.5 * fz * fga * pa_zz * pb_xz * fx

                     - fz * fga * pa_z * fx * pb_xxx

                     - 3.0 * pa_xxzz * pb_xz * fz * fgb

                     + 9.0 * pa_xx * fz * fx * fx * pb_xz

                     + 18.0 * pa_x * fx * fx * fz * pb_xxz

                     + 6.0 * fx * fx * pa_z * fz * pb_xxx

                     + 21.0 * pa_xxzz * fz * pb_xz * fx

                     + 14.0 * pa_xxz * fz * fx * pb_xxx

                     + 42.0 * pa_xzz * fx * fz * pb_xxz

                     - fx * fz * fga * pb_xxxz

                     - pa_xx * fz * fga * pb_xxxz

                     + 3.0 * fx * fx * fz * pb_xxxz

                     - fz * fga * pa_zz * pb_xxxz

                     + 7.0 * pa_xx * fz * fx * pb_xxxz

                     + 7.0 * fx * pa_zz * fz * pb_xxxz

                     + 16.0 * pa_xxzz * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xxyy_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxzz,
                                     double pa_xzz,
                                     double pa_zz,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxyy,
                                     double pb_xyy,
                                     double pb_yy,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * fx * fx * fx * pa_zz

                     + 0.125 * pa_xx * fx * fx * fx

                     + 0.5 * pa_x * fx * fx * fx * pb_x

                     + 0.375 * fx * fx * fx * pb_yy

                     + 0.25 * pa_xxzz * fx * fx

                     + pa_xzz * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_zz * pb_yy

                     + 0.125 * fx * fx * fx * pb_xx

                     + 0.25 * pa_xx * fx * fx * pb_xx

                     + 0.25 * pa_xx * fx * fx * pb_yy

                     + pa_x * fx * fx * pb_xyy

                     + 0.25 * fx * fx * pa_zz * pb_xx

                     + 0.5 * pa_xxzz * pb_xx * fx

                     + 0.5 * pa_xxzz * fx * pb_yy

                     + 2.0 * pa_xzz * fx * pb_xyy

                     + 0.25 * fx * fx * pb_xxyy

                     + 0.5 * pa_xx * fx * pb_xxyy

                     + 0.5 * fx * pa_zz * pb_xxyy

                     + pa_xxzz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xxyy_r_0(double fga,
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
                                     double pb_xxyy,
                                     double pb_xyy,
                                     double pb_yy,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * fx * fz * fgb

                     - 0.5 * fx * fx * fx * fz * fga

                     - fx * fx * pa_zz * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     + 3.75 * fx * fx * fx * pa_zz * fz

                     - 0.5 * pa_xx * fx * fx * fz * fgb

                     - 0.25 * pa_xx * fz * fga * fx * fx

                     - pa_x * fx * fx * pb_x * fz * fgb

                     - pa_x * fx * fx * fz * fga * pb_x

                     - fx * fx * fz * fga * pb_yy

                     - 0.25 * fz * fga * pa_zz * fx * fx

                     - pa_xxzz * fx * fz * fgb

                     + 1.25 * pa_xx * fz * fx * fx * fx

                     - 2.0 * pa_xzz * fx * pb_x * fz * fgb

                     + 5.0 * pa_x * fx * fx * fx * fz * pb_x

                     + 3.75 * fx * fx * fx * fz * pb_yy

                     + 3.0 * pa_xxzz * fz * fx * fx

                     + 12.0 * pa_xzz * fx * fx * fz * pb_x

                     + 9.0 * fx * fx * pa_zz * fz * pb_yy

                     - 0.25 * fx * fx * pb_xx * fz * fgb

                     - 0.25 * fx * fx * fz * fgb * pb_yy

                     - 0.5 * fx * fx * fz * fga * pb_xx

                     - 0.5 * pa_xx * fx * pb_xx * fz * fgb

                     - 0.5 * pa_xx * fx * fz * fgb * pb_yy

                     - 0.5 * pa_xx * fz * fga * pb_xx * fx

                     - 0.5 * pa_xx * fz * fga * fx * pb_yy

                     - 2.0 * pa_x * fx * fz * fga * pb_xyy

                     - 0.5 * fx * pa_zz * pb_xx * fz * fgb

                     - 0.5 * fx * pa_zz * fz * fgb * pb_yy

                     + 1.25 * fx * fx * fx * fz * pb_xx

                     - 0.5 * fz * fga * pa_zz * pb_xx * fx

                     - 0.5 * fz * fga * pa_zz * fx * pb_yy

                     - pa_xxzz * pb_xx * fz * fgb

                     - pa_xxzz * fz * fgb * pb_yy

                     + 3.0 * pa_xx * fz * fx * fx * pb_xx

                     + 3.0 * pa_xx * fz * fx * fx * pb_yy

                     + 12.0 * pa_x * fx * fx * fz * pb_xyy

                     + 3.0 * fx * fx * pa_zz * fz * pb_xx

                     + 7.0 * pa_xxzz * fz * pb_xx * fx

                     + 7.0 * pa_xxzz * fz * fx * pb_yy

                     + 28.0 * pa_xzz * fx * fz * pb_xyy

                     - fx * fz * fga * pb_xxyy

                     - pa_xx * fz * fga * pb_xxyy

                     + 3.0 * fx * fx * fz * pb_xxyy

                     - fz * fga * pa_zz * pb_xxyy

                     + 7.0 * pa_xx * fz * fx * pb_xxyy

                     + 7.0 * fx * pa_zz * fz * pb_xxyy

                     + 16.0 * pa_xxzz * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xxyz_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxz,
                                     double pa_xxzz,
                                     double pa_xz,
                                     double pa_xzz,
                                     double pa_z,
                                     double pa_zz,
                                     double pb_xxy,
                                     double pb_xxyz,
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_y,
                                     double pb_yz,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * fx * pa_z * pb_y

                     + 0.375 * fx * fx * fx * pb_yz

                     + 0.5 * pa_xxz * fx * fx * pb_y

                     + 2.0 * pa_xz * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_zz * pb_yz

                     + 0.25 * pa_xx * fx * fx * pb_yz

                     + pa_x * fx * fx * pb_xyz

                     + 0.5 * fx * fx * pa_z * pb_xxy

                     + 0.5 * pa_xxzz * fx * pb_yz

                     + pa_xxz * fx * pb_xxy

                     + 2.0 * pa_xzz * fx * pb_xyz

                     + 0.25 * fx * fx * pb_xxyz

                     + 0.5 * pa_xx * fx * pb_xxyz

                     + 0.5 * fx * pa_zz * pb_xxyz

                     + pa_xxzz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xxyz_r_0(double fga,
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
                                     double pb_xxy,
                                     double pb_xxyz,
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_y,
                                     double pb_yz,
                                     double r_0_0)
    {
        return r_0_0 * (7.5 * fx * fx * fx * pa_z * fz * pb_y

                     - 0.5 * fx * fx * pa_z * fz * fgb * pb_y

                     - fx * fx * fz * fga * pb_yz

                     - 0.5 * fz * fga * pa_z * fx * fx * pb_y

                     - pa_xxz * fx * fz * fgb * pb_y

                     + 3.75 * fx * fx * fx * fz * pb_yz

                     + 6.0 * pa_xxz * fz * fx * fx * pb_y

                     + 24.0 * pa_xz * fx * fx * fz * pb_xy

                     + 9.0 * fx * fx * pa_zz * fz * pb_yz

                     - 0.25 * fx * fx * fz * fgb * pb_yz

                     - 0.5 * pa_xx * fx * fz * fgb * pb_yz

                     - 0.5 * pa_xx * fz * fga * fx * pb_yz

                     - 2.0 * pa_x * fx * fz * fga * pb_xyz

                     - 0.5 * fx * pa_zz * fz * fgb * pb_yz

                     - 0.5 * fz * fga * pa_zz * fx * pb_yz

                     - fz * fga * pa_z * fx * pb_xxy

                     - pa_xxzz * fz * fgb * pb_yz

                     + 3.0 * pa_xx * fz * fx * fx * pb_yz

                     + 12.0 * pa_x * fx * fx * fz * pb_xyz

                     + 6.0 * fx * fx * pa_z * fz * pb_xxy

                     + 7.0 * pa_xxzz * fz * fx * pb_yz

                     + 14.0 * pa_xxz * fz * fx * pb_xxy

                     + 28.0 * pa_xzz * fx * fz * pb_xyz

                     - fx * fz * fga * pb_xxyz

                     - pa_xx * fz * fga * pb_xxyz

                     + 3.0 * fx * fx * fz * pb_xxyz

                     - fz * fga * pa_zz * pb_xxyz

                     + 7.0 * pa_xx * fz * fx * pb_xxyz

                     + 7.0 * fx * pa_zz * fz * pb_xxyz

                     + 16.0 * pa_xxzz * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xxzz_s_0(double fx,
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
                                     double pb_xxzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.5625 * fx * fx * fx * fx

                     + 0.375 * pa_xx * fx * fx * fx

                     + 1.5 * pa_x * fx * fx * fx * pb_x

                     + 0.375 * fx * fx * fx * pa_zz

                     + 1.5 * fx * fx * fx * pa_z * pb_z

                     + 0.375 * fx * fx * fx * pb_xx

                     + 0.375 * fx * fx * fx * pb_zz

                     + 0.25 * pa_xxzz * fx * fx

                     + pa_xxz * fx * fx * pb_z

                     + 0.75 * pa_xx * fx * fx * pb_xx

                     + pa_xzz * fx * fx * pb_x

                     + 4.0 * pa_xz * fx * fx * pb_xz

                     + 0.75 * fx * fx * pa_zz * pb_zz

                     + 0.25 * pa_xx * fx * fx * pb_zz

                     + pa_x * fx * fx * pb_xzz

                     + 0.25 * fx * fx * pa_zz * pb_xx

                     + fx * fx * pa_z * pb_xxz

                     + 0.5 * pa_xxzz * pb_xx * fx

                     + 0.5 * pa_xxzz * fx * pb_zz

                     + 2.0 * pa_xxz * fx * pb_xxz

                     + 2.0 * pa_xzz * fx * pb_xzz

                     + 0.25 * fx * fx * pb_xxzz

                     + 0.5 * pa_xx * fx * pb_xxzz

                     + 0.5 * fx * pa_zz * pb_xxzz

                     + pa_xxzz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xxzz_r_0(double fga,
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
                                     double pb_xxzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (4.5 * fx * fx * fx * fx * fz

                     - 0.75 * fx * fx * fx * fz * fgb

                     - 0.75 * fx * fx * fx * fz * fga

                     - pa_xx * fx * fx * fz * fgb

                     - fx * fx * pa_zz * fz * fgb

                     + 3.75 * pa_xx * fx * fx * fx * fz

                     + 15.0 * pa_x * fx * fx * fx * fz * pb_x

                     + 3.75 * fx * fx * fx * pa_zz * fz

                     + 15.0 * fx * fx * fx * pa_z * fz * pb_z

                     - 0.25 * pa_xx * fz * fga * fx * fx

                     - pa_x * fx * fx * pb_x * fz * fgb

                     - pa_x * fx * fx * fz * fga * pb_x

                     - fx * fx * pa_z * fz * fgb * pb_z

                     - fx * fx * fz * fga * pb_zz

                     - 0.25 * fz * fga * pa_zz * fx * fx

                     - fz * fga * pa_z * fx * fx * pb_z

                     - fz * fga * fx * fx * pb_xx

                     - pa_xxzz * fx * fz * fgb

                     - 2.0 * pa_xxz * fx * fz * fgb * pb_z

                     - 2.0 * pa_xzz * fx * pb_x * fz * fgb

                     + 3.75 * fx * fx * fx * fz * pb_xx

                     + 3.75 * fx * fx * fx * fz * pb_zz

                     + 3.0 * pa_xxzz * fz * fx * fx

                     + 12.0 * pa_xxz * fz * fx * fx * pb_z

                     + 9.0 * pa_xx * fx * fx * fz * pb_xx

                     + 12.0 * pa_xzz * fx * fx * fz * pb_x

                     + 48.0 * pa_xz * fx * fx * fz * pb_xz

                     + 9.0 * fx * fx * pa_zz * fz * pb_zz

                     - 0.25 * fx * fx * pb_xx * fz * fgb

                     - 0.25 * fx * fx * fz * fgb * pb_zz

                     - 0.5 * pa_xx * fx * pb_xx * fz * fgb

                     - 0.5 * pa_xx * fx * fz * fgb * pb_zz

                     - 0.5 * pa_xx * fz * fga * pb_xx * fx

                     - 0.5 * pa_xx * fz * fga * fx * pb_zz

                     - 2.0 * pa_x * fx * fz * fga * pb_xzz

                     - 0.5 * fx * pa_zz * pb_xx * fz * fgb

                     - 0.5 * fx * pa_zz * fz * fgb * pb_zz

                     - 0.5 * fz * fga * pa_zz * pb_xx * fx

                     - 0.5 * fz * fga * pa_zz * fx * pb_zz

                     - 2.0 * fz * fga * pa_z * fx * pb_xxz

                     - pa_xxzz * pb_xx * fz * fgb

                     - pa_xxzz * fz * fgb * pb_zz

                     + 3.0 * pa_xx * fz * fx * fx * pb_zz

                     + 12.0 * pa_x * fx * fx * fz * pb_xzz

                     + 3.0 * fx * fx * pa_zz * fz * pb_xx

                     + 12.0 * fx * fx * pa_z * fz * pb_xxz

                     + 7.0 * pa_xxzz * fz * pb_xx * fx

                     + 7.0 * pa_xxzz * fz * fx * pb_zz

                     + 28.0 * pa_xxz * fz * fx * pb_xxz

                     + 28.0 * pa_xzz * fx * fz * pb_xzz

                     - fx * fz * fga * pb_xxzz

                     - pa_xx * fz * fga * pb_xxzz

                     + 3.0 * fx * fx * fz * pb_xxzz

                     - fz * fga * pa_zz * pb_xxzz

                     + 7.0 * pa_xx * fz * fx * pb_xxzz

                     + 7.0 * fx * pa_zz * fz * pb_xxzz

                     + 16.0 * pa_xxzz * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xyyy_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxzz,
                                     double pa_xzz,
                                     double pa_zz,
                                     double pb_xy,
                                     double pb_xyyy,
                                     double pb_y,
                                     double pb_yyy,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * fx * pb_y

                     + 1.5 * pa_xzz * fx * fx * pb_y

                     + 0.375 * fx * fx * fx * pb_xy

                     + 0.75 * pa_xx * fx * fx * pb_xy

                     + 0.5 * pa_x * fx * fx * pb_yyy

                     + 0.75 * fx * fx * pa_zz * pb_xy

                     + 1.5 * pa_xxzz * pb_xy * fx

                     + pa_xzz * fx * pb_yyy

                     + 0.25 * fx * fx * pb_xyyy

                     + 0.5 * pa_xx * fx * pb_xyyy

                     + 0.5 * fx * pa_zz * pb_xyyy

                     + pa_xxzz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xyyy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxzz,
                                     double pa_xzz,
                                     double pa_zz,
                                     double pb_xy,
                                     double pb_xyyy,
                                     double pb_y,
                                     double pb_yyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fx * pb_y * fz * fgb

                     - 1.5 * pa_x * fx * fx * fz * fga * pb_y

                     - 3.0 * pa_xzz * fx * pb_y * fz * fgb

                     + 7.5 * pa_x * fx * fx * fx * fz * pb_y

                     + 18.0 * pa_xzz * fx * fx * fz * pb_y

                     - 0.75 * fx * fx * pb_xy * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_xy

                     - 1.5 * pa_xx * fx * pb_xy * fz * fgb

                     - 1.5 * pa_xx * fz * fga * pb_xy * fx

                     - pa_x * fx * fz * fga * pb_yyy

                     - 1.5 * fx * pa_zz * pb_xy * fz * fgb

                     + 3.75 * fx * fx * fx * fz * pb_xy

                     - 1.5 * fz * fga * pa_zz * pb_xy * fx

                     - 3.0 * pa_xxzz * pb_xy * fz * fgb

                     + 9.0 * pa_xx * fz * fx * fx * pb_xy

                     + 6.0 * pa_x * fx * fx * fz * pb_yyy

                     + 9.0 * fx * fx * pa_zz * fz * pb_xy

                     + 21.0 * pa_xxzz * fz * pb_xy * fx

                     + 14.0 * pa_xzz * fx * fz * pb_yyy

                     - fx * fz * fga * pb_xyyy

                     - pa_xx * fz * fga * pb_xyyy

                     + 3.0 * fx * fx * fz * pb_xyyy

                     - fz * fga * pa_zz * pb_xyyy

                     + 7.0 * pa_xx * fz * fx * pb_xyyy

                     + 7.0 * fx * pa_zz * fz * pb_xyyy

                     + 16.0 * pa_xxzz * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xyyz_s_0(double fx,
                                     double pa_x,
                                     double pa_xx,
                                     double pa_xxz,
                                     double pa_xxzz,
                                     double pa_xz,
                                     double pa_xzz,
                                     double pa_z,
                                     double pa_zz,
                                     double pb_x,
                                     double pb_xyy,
                                     double pb_xyyz,
                                     double pb_xz,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xz * fx * fx * fx

                     + 0.25 * pa_x * fx * fx * fx * pb_z

                     + 0.25 * fx * fx * fx * pa_z * pb_x

                     + 0.5 * pa_xxz * fx * fx * pb_x

                     + 0.5 * pa_xzz * fx * fx * pb_z

                     + pa_xz * fx * fx * pb_yy

                     + 0.125 * fx * fx * fx * pb_xz

                     + 0.25 * pa_xx * fx * fx * pb_xz

                     + 0.5 * pa_x * fx * fx * pb_yyz

                     + 0.25 * fx * fx * pa_zz * pb_xz

                     + 0.5 * fx * fx * pa_z * pb_xyy

                     + 0.5 * pa_xxzz * pb_xz * fx

                     + pa_xxz * fx * pb_xyy

                     + pa_xzz * fx * pb_yyz

                     + 0.25 * fx * fx * pb_xyyz

                     + 0.5 * pa_xx * fx * pb_xyyz

                     + 0.5 * fx * pa_zz * pb_xyyz

                     + pa_xxzz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xyyz_r_0(double fga,
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
                                     double pb_xyy,
                                     double pb_xyyz,
                                     double pb_xz,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- pa_xz * fx * fx * fz * fgb

                     + 5.0 * pa_xz * fx * fx * fx * fz

                     - 0.5 * pa_x * fx * fx * fz * fgb * pb_z

                     - 0.5 * pa_x * fx * fx * fz * fga * pb_z

                     - 0.5 * fx * fx * pa_z * pb_x * fz * fgb

                     - 0.5 * fz * fga * pa_z * fx * fx * pb_x

                     - pa_xxz * fx * pb_x * fz * fgb

                     - pa_xzz * fx * fz * fgb * pb_z

                     + 2.5 * pa_x * fx * fx * fx * fz * pb_z

                     + 2.5 * fx * fx * fx * pa_z * fz * pb_x

                     + 6.0 * pa_xxz * fz * fx * fx * pb_x

                     + 6.0 * pa_xzz * fx * fx * fz * pb_z

                     + 12.0 * pa_xz * fx * fx * fz * pb_yy

                     - 0.25 * fx * fx * pb_xz * fz * fgb

                     - 0.5 * fx * fx * fz * fga * pb_xz

                     - 0.5 * pa_xx * fx * pb_xz * fz * fgb

                     - 0.5 * pa_xx * fz * fga * pb_xz * fx

                     - pa_x * fx * fz * fga * pb_yyz

                     - 0.5 * fx * pa_zz * pb_xz * fz * fgb

                     + 1.25 * fx * fx * fx * fz * pb_xz

                     - 0.5 * fz * fga * pa_zz * pb_xz * fx

                     - fz * fga * pa_z * fx * pb_xyy

                     - pa_xxzz * pb_xz * fz * fgb

                     + 3.0 * pa_xx * fz * fx * fx * pb_xz

                     + 6.0 * pa_x * fx * fx * fz * pb_yyz

                     + 3.0 * fx * fx * pa_zz * fz * pb_xz

                     + 6.0 * fx * fx * pa_z * fz * pb_xyy

                     + 7.0 * pa_xxzz * fz * pb_xz * fx

                     + 14.0 * pa_xxz * fz * fx * pb_xyy

                     + 14.0 * pa_xzz * fx * fz * pb_yyz

                     - fx * fz * fga * pb_xyyz

                     - pa_xx * fz * fga * pb_xyyz

                     + 3.0 * fx * fx * fz * pb_xyyz

                     - fz * fga * pa_zz * pb_xyyz

                     + 7.0 * pa_xx * fz * fx * pb_xyyz

                     + 7.0 * fx * pa_zz * fz * pb_xyyz

                     + 16.0 * pa_xxzz * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xyzz_s_0(double fx,
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
                                     double pb_xyzz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_yzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * fx * pb_y

                     + 0.375 * fx * fx * fx * pb_xy

                     + 0.75 * pa_xx * fx * fx * pb_xy

                     + 0.5 * pa_xzz * fx * fx * pb_y

                     + 2.0 * pa_xz * fx * fx * pb_yz

                     + 0.5 * pa_x * fx * fx * pb_yzz

                     + 0.25 * fx * fx * pa_zz * pb_xy

                     + fx * fx * pa_z * pb_xyz

                     + 0.5 * pa_xxzz * pb_xy * fx

                     + 2.0 * pa_xxz * fx * pb_xyz

                     + pa_xzz * fx * pb_yzz

                     + 0.25 * fx * fx * pb_xyzz

                     + 0.5 * pa_xx * fx * pb_xyzz

                     + 0.5 * fx * pa_zz * pb_xyzz

                     + pa_xxzz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xyzz_r_0(double fga,
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
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_xyzz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_yzz,
                                     double r_0_0)
    {
        return r_0_0 * (7.5 * pa_x * fx * fx * fx * fz * pb_y

                     - 0.5 * pa_x * fx * fx * pb_y * fz * fgb

                     - 0.5 * pa_x * fx * fx * fz * fga * pb_y

                     - fz * fga * fx * fx * pb_xy

                     - pa_xzz * fx * pb_y * fz * fgb

                     + 3.75 * fx * fx * fx * fz * pb_xy

                     + 9.0 * pa_xx * fx * fx * fz * pb_xy

                     + 6.0 * pa_xzz * fx * fx * fz * pb_y

                     + 24.0 * pa_xz * fx * fx * fz * pb_yz

                     - 0.25 * fx * fx * pb_xy * fz * fgb

                     - 0.5 * pa_xx * fx * pb_xy * fz * fgb

                     - 0.5 * pa_xx * fz * fga * pb_xy * fx

                     - pa_x * fx * fz * fga * pb_yzz

                     - 0.5 * fx * pa_zz * pb_xy * fz * fgb

                     - 0.5 * fz * fga * pa_zz * pb_xy * fx

                     - 2.0 * fz * fga * pa_z * fx * pb_xyz

                     - pa_xxzz * pb_xy * fz * fgb

                     + 6.0 * pa_x * fx * fx * fz * pb_yzz

                     + 3.0 * fx * fx * pa_zz * fz * pb_xy

                     + 12.0 * fx * fx * pa_z * fz * pb_xyz

                     + 7.0 * pa_xxzz * fz * pb_xy * fx

                     + 28.0 * pa_xxz * fz * fx * pb_xyz

                     + 14.0 * pa_xzz * fx * fz * pb_yzz

                     - fx * fz * fga * pb_xyzz

                     - pa_xx * fz * fga * pb_xyzz

                     + 3.0 * fx * fx * fz * pb_xyzz

                     - fz * fga * pa_zz * pb_xyzz

                     + 7.0 * pa_xx * fz * fx * pb_xyzz

                     + 7.0 * fx * pa_zz * fz * pb_xyzz

                     + 16.0 * pa_xxzz * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xzzz_s_0(double fx,
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
                                     double pb_xzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xz * fx * fx * fx

                     + 2.25 * pa_x * fx * fx * fx * pb_z

                     + 0.75 * fx * fx * fx * pa_z * pb_x

                     + 1.125 * fx * fx * fx * pb_xz

                     + 1.5 * pa_xxz * fx * fx * pb_x

                     + 2.25 * pa_xx * fx * fx * pb_xz

                     + 1.5 * pa_xzz * fx * fx * pb_z

                     + 3.0 * pa_xz * fx * fx * pb_zz

                     + 0.5 * pa_x * fx * fx * pb_zzz

                     + 0.75 * fx * fx * pa_zz * pb_xz

                     + 1.5 * fx * fx * pa_z * pb_xzz

                     + 1.5 * pa_xxzz * pb_xz * fx

                     + 3.0 * pa_xxz * fx * pb_xzz

                     + pa_xzz * fx * pb_zzz

                     + 0.25 * fx * fx * pb_xzzz

                     + 0.5 * pa_xx * fx * pb_xzzz

                     + 0.5 * fx * pa_zz * pb_xzzz

                     + pa_xxzz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xzzz_r_0(double fga,
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
                                     double pb_xzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xz * fx * fx * fz * fgb

                     + 15.0 * pa_xz * fx * fx * fx * fz

                     + 22.5 * pa_x * fx * fx * fx * fz * pb_z

                     - 1.5 * pa_x * fx * fx * pb_z * fz * fgb

                     - 1.5 * pa_x * fx * fx * fz * fga * pb_z

                     - 1.5 * fx * fx * pa_z * pb_x * fz * fgb

                     - 1.5 * fz * fga * pa_z * fx * fx * pb_x

                     - 3.0 * fz * fga * fx * fx * pb_xz

                     - 3.0 * pa_xxz * fx * pb_x * fz * fgb

                     - 3.0 * pa_xzz * fx * pb_z * fz * fgb

                     + 7.5 * fx * fx * fx * pa_z * fz * pb_x

                     + 11.25 * fx * fx * fx * fz * pb_xz

                     + 18.0 * pa_xxz * fz * fx * fx * pb_x

                     + 27.0 * pa_xx * fx * fx * fz * pb_xz

                     + 18.0 * pa_xzz * fx * fx * fz * pb_z

                     + 36.0 * pa_xz * fx * fx * fz * pb_zz

                     - 0.75 * fx * fx * pb_xz * fz * fgb

                     - 1.5 * pa_xx * fx * pb_xz * fz * fgb

                     - 1.5 * pa_xx * fz * fga * pb_xz * fx

                     - pa_x * fx * fz * fga * pb_zzz

                     - 1.5 * fx * pa_zz * pb_xz * fz * fgb

                     - 1.5 * fz * fga * pa_zz * pb_xz * fx

                     - 3.0 * fz * fga * pa_z * fx * pb_xzz

                     - 3.0 * pa_xxzz * pb_xz * fz * fgb

                     + 6.0 * pa_x * fx * fx * fz * pb_zzz

                     + 9.0 * fx * fx * pa_zz * fz * pb_xz

                     + 18.0 * fx * fx * pa_z * fz * pb_xzz

                     + 21.0 * pa_xxzz * fz * pb_xz * fx

                     + 42.0 * pa_xxz * fz * fx * pb_xzz

                     + 14.0 * pa_xzz * fx * fz * pb_zzz

                     - fx * fz * fga * pb_xzzz

                     - pa_xx * fz * fga * pb_xzzz

                     + 3.0 * fx * fx * fz * pb_xzzz

                     - fz * fga * pa_zz * pb_xzzz

                     + 7.0 * pa_xx * fz * fx * pb_xzzz

                     + 7.0 * fx * pa_zz * fz * pb_xzzz

                     + 16.0 * pa_xxzz * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_yyyy_s_0(double fx,
                                     double pa_xx,
                                     double pa_xxzz,
                                     double pa_zz,
                                     double pb_yy,
                                     double pb_yyyy,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * pa_xx * fx * fx * fx

                     + 0.375 * fx * fx * fx * pa_zz

                     + 0.75 * pa_xxzz * fx * fx

                     + 0.75 * fx * fx * fx * pb_yy

                     + 1.5 * pa_xx * fx * fx * pb_yy

                     + 1.5 * fx * fx * pa_zz * pb_yy

                     + 3.0 * pa_xxzz * pb_yy * fx

                     + 0.25 * fx * fx * pb_yyyy

                     + 0.5 * pa_xx * fx * pb_yyyy

                     + 0.5 * fx * pa_zz * pb_yyyy

                     + pa_xxzz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_yyyy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xx,
                                     double pa_xxzz,
                                     double pa_zz,
                                     double pb_yy,
                                     double pb_yyyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fx * fz * fgb

                     - 0.75 * fx * fx * fx * fz * fga

                     - 1.5 * pa_xx * fx * fx * fz * fgb

                     - 0.75 * pa_xx * fz * fga * fx * fx

                     - 1.5 * fx * fx * pa_zz * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     - 0.75 * fz * fga * pa_zz * fx * fx

                     - 3.0 * pa_xxzz * fx * fz * fgb

                     + 3.75 * pa_xx * fz * fx * fx * fx

                     + 3.75 * fx * fx * fx * pa_zz * fz

                     + 9.0 * pa_xxzz * fz * fx * fx

                     - 1.5 * fx * fx * pb_yy * fz * fgb

                     - 3.0 * fx * fx * fz * fga * pb_yy

                     - 3.0 * pa_xx * fx * pb_yy * fz * fgb

                     - 3.0 * pa_xx * fz * fga * pb_yy * fx

                     - 3.0 * fx * pa_zz * pb_yy * fz * fgb

                     + 7.5 * fx * fx * fx * fz * pb_yy

                     - 3.0 * fz * fga * pa_zz * pb_yy * fx

                     - 6.0 * pa_xxzz * pb_yy * fz * fgb

                     + 18.0 * pa_xx * fz * fx * fx * pb_yy

                     + 18.0 * fx * fx * pa_zz * fz * pb_yy

                     + 42.0 * pa_xxzz * fz * pb_yy * fx

                     - fx * fz * fga * pb_yyyy

                     - pa_xx * fz * fga * pb_yyyy

                     + 3.0 * fx * fx * fz * pb_yyyy

                     - fz * fga * pa_zz * pb_yyyy

                     + 7.0 * pa_xx * fz * fx * pb_yyyy

                     + 7.0 * fx * pa_zz * fz * pb_yyyy

                     + 16.0 * pa_xxzz * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_yyyz_s_0(double fx,
                                     double pa_xx,
                                     double pa_xxz,
                                     double pa_xxzz,
                                     double pa_z,
                                     double pa_zz,
                                     double pb_y,
                                     double pb_yyy,
                                     double pb_yyyz,
                                     double pb_yz,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * fx * pa_z * pb_y

                     + 1.5 * pa_xxz * fx * fx * pb_y

                     + 0.375 * fx * fx * fx * pb_yz

                     + 0.75 * pa_xx * fx * fx * pb_yz

                     + 0.75 * fx * fx * pa_zz * pb_yz

                     + 0.5 * fx * fx * pa_z * pb_yyy

                     + 1.5 * pa_xxzz * pb_yz * fx

                     + pa_xxz * fx * pb_yyy

                     + 0.25 * fx * fx * pb_yyyz

                     + 0.5 * pa_xx * fx * pb_yyyz

                     + 0.5 * fx * pa_zz * pb_yyyz

                     + pa_xxzz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_yyyz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xx,
                                     double pa_xxz,
                                     double pa_xxzz,
                                     double pa_z,
                                     double pa_zz,
                                     double pb_y,
                                     double pb_yyy,
                                     double pb_yyyz,
                                     double pb_yz,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_z * pb_y * fz * fgb

                     - 1.5 * fz * fga * pa_z * fx * fx * pb_y

                     - 3.0 * pa_xxz * fx * pb_y * fz * fgb

                     + 7.5 * fx * fx * fx * pa_z * fz * pb_y

                     + 18.0 * pa_xxz * fz * fx * fx * pb_y

                     - 0.75 * fx * fx * pb_yz * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_yz

                     - 1.5 * pa_xx * fx * pb_yz * fz * fgb

                     - 1.5 * pa_xx * fz * fga * pb_yz * fx

                     - 1.5 * fx * pa_zz * pb_yz * fz * fgb

                     + 3.75 * fx * fx * fx * fz * pb_yz

                     - 1.5 * fz * fga * pa_zz * pb_yz * fx

                     - fz * fga * pa_z * fx * pb_yyy

                     - 3.0 * pa_xxzz * pb_yz * fz * fgb

                     + 9.0 * pa_xx * fz * fx * fx * pb_yz

                     + 9.0 * fx * fx * pa_zz * fz * pb_yz

                     + 6.0 * fx * fx * pa_z * fz * pb_yyy

                     + 21.0 * pa_xxzz * fz * pb_yz * fx

                     + 14.0 * pa_xxz * fz * fx * pb_yyy

                     - fx * fz * fga * pb_yyyz

                     - pa_xx * fz * fga * pb_yyyz

                     + 3.0 * fx * fx * fz * pb_yyyz

                     - fz * fga * pa_zz * pb_yyyz

                     + 7.0 * pa_xx * fz * fx * pb_yyyz

                     + 7.0 * fx * pa_zz * fz * pb_yyyz

                     + 16.0 * pa_xxzz * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_yyzz_s_0(double fx,
                                     double pa_xx,
                                     double pa_xxz,
                                     double pa_xxzz,
                                     double pa_z,
                                     double pa_zz,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_yyzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * pa_xx * fx * fx * fx

                     + 0.125 * fx * fx * fx * pa_zz

                     + 0.5 * fx * fx * fx * pa_z * pb_z

                     + 0.375 * fx * fx * fx * pb_yy

                     + 0.25 * pa_xxzz * fx * fx

                     + pa_xxz * fx * fx * pb_z

                     + 0.75 * pa_xx * fx * fx * pb_yy

                     + 0.125 * fx * fx * fx * pb_zz

                     + 0.25 * pa_xx * fx * fx * pb_zz

                     + 0.25 * fx * fx * pa_zz * pb_yy

                     + 0.25 * fx * fx * pa_zz * pb_zz

                     + fx * fx * pa_z * pb_yyz

                     + 0.5 * pa_xxzz * pb_yy * fx

                     + 0.5 * pa_xxzz * fx * pb_zz

                     + 2.0 * pa_xxz * fx * pb_yyz

                     + 0.25 * fx * fx * pb_yyzz

                     + 0.5 * pa_xx * fx * pb_yyzz

                     + 0.5 * fx * pa_zz * pb_yyzz

                     + pa_xxzz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_yyzz_r_0(double fga,
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
                                     double pb_yyzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * fx * fz * fgb

                     - 0.5 * fz * fga * fx * fx * fx

                     - pa_xx * fx * fx * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     + 3.75 * pa_xx * fx * fx * fx * fz

                     - 0.25 * pa_xx * fz * fga * fx * fx

                     - 0.5 * fx * fx * pa_zz * fz * fgb

                     - fx * fx * pa_z * fz * fgb * pb_z

                     - 0.25 * fz * fga * pa_zz * fx * fx

                     - fz * fga * pa_z * fx * fx * pb_z

                     - fz * fga * fx * fx * pb_yy

                     - pa_xxzz * fx * fz * fgb

                     - 2.0 * pa_xxz * fx * fz * fgb * pb_z

                     + 1.25 * fx * fx * fx * pa_zz * fz

                     + 5.0 * fx * fx * fx * pa_z * fz * pb_z

                     + 3.75 * fx * fx * fx * fz * pb_yy

                     + 3.0 * pa_xxzz * fz * fx * fx

                     + 12.0 * pa_xxz * fz * fx * fx * pb_z

                     + 9.0 * pa_xx * fx * fx * fz * pb_yy

                     - 0.25 * fx * fx * pb_yy * fz * fgb

                     - 0.25 * fx * fx * fz * fgb * pb_zz

                     - 0.5 * fx * fx * fz * fga * pb_zz

                     - 0.5 * pa_xx * fx * pb_yy * fz * fgb

                     - 0.5 * pa_xx * fx * fz * fgb * pb_zz

                     - 0.5 * pa_xx * fz * fga * pb_yy * fx

                     - 0.5 * pa_xx * fz * fga * fx * pb_zz

                     - 0.5 * fx * pa_zz * pb_yy * fz * fgb

                     - 0.5 * fx * pa_zz * fz * fgb * pb_zz

                     + 1.25 * fx * fx * fx * fz * pb_zz

                     - 0.5 * fz * fga * pa_zz * pb_yy * fx

                     - 0.5 * fz * fga * pa_zz * fx * pb_zz

                     - 2.0 * fz * fga * pa_z * fx * pb_yyz

                     - pa_xxzz * pb_yy * fz * fgb

                     - pa_xxzz * fz * fgb * pb_zz

                     + 3.0 * pa_xx * fz * fx * fx * pb_zz

                     + 3.0 * fx * fx * pa_zz * fz * pb_yy

                     + 3.0 * fx * fx * pa_zz * fz * pb_zz

                     + 12.0 * fx * fx * pa_z * fz * pb_yyz

                     + 7.0 * pa_xxzz * fz * pb_yy * fx

                     + 7.0 * pa_xxzz * fz * fx * pb_zz

                     + 28.0 * pa_xxz * fz * fx * pb_yyz

                     - fx * fz * fga * pb_yyzz

                     - pa_xx * fz * fga * pb_yyzz

                     + 3.0 * fx * fx * fz * pb_yyzz

                     - fz * fga * pa_zz * pb_yyzz

                     + 7.0 * pa_xx * fz * fx * pb_yyzz

                     + 7.0 * fx * pa_zz * fz * pb_yyzz

                     + 16.0 * pa_xxzz * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_yzzz_s_0(double fx,
                                     double pa_xx,
                                     double pa_xxz,
                                     double pa_xxzz,
                                     double pa_z,
                                     double pa_zz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_yzzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * fx * pa_z * pb_y

                     + 1.125 * fx * fx * fx * pb_yz

                     + 1.5 * pa_xxz * fx * fx * pb_y

                     + 2.25 * pa_xx * fx * fx * pb_yz

                     + 0.75 * fx * fx * pa_zz * pb_yz

                     + 1.5 * fx * fx * pa_z * pb_yzz

                     + 1.5 * pa_xxzz * pb_yz * fx

                     + 3.0 * pa_xxz * fx * pb_yzz

                     + 0.25 * fx * fx * pb_yzzz

                     + 0.5 * pa_xx * fx * pb_yzzz

                     + 0.5 * fx * pa_zz * pb_yzzz

                     + pa_xxzz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_yzzz_r_0(double fga,
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
                                     double pb_yzzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_z * pb_y * fz * fgb

                     - 1.5 * fz * fga * pa_z * fx * fx * pb_y

                     - 3.0 * fz * fga * fx * fx * pb_yz

                     - 3.0 * pa_xxz * fx * pb_y * fz * fgb

                     + 7.5 * fx * fx * fx * pa_z * fz * pb_y

                     + 11.25 * fx * fx * fx * fz * pb_yz

                     + 18.0 * pa_xxz * fz * fx * fx * pb_y

                     + 27.0 * pa_xx * fx * fx * fz * pb_yz

                     - 0.75 * fx * fx * pb_yz * fz * fgb

                     - 1.5 * pa_xx * fx * pb_yz * fz * fgb

                     - 1.5 * pa_xx * fz * fga * pb_yz * fx

                     - 1.5 * fx * pa_zz * pb_yz * fz * fgb

                     - 1.5 * fz * fga * pa_zz * pb_yz * fx

                     - 3.0 * fz * fga * pa_z * fx * pb_yzz

                     - 3.0 * pa_xxzz * pb_yz * fz * fgb

                     + 9.0 * fx * fx * pa_zz * fz * pb_yz

                     + 18.0 * fx * fx * pa_z * fz * pb_yzz

                     + 21.0 * pa_xxzz * fz * pb_yz * fx

                     + 42.0 * pa_xxz * fz * fx * pb_yzz

                     - fx * fz * fga * pb_yzzz

                     - pa_xx * fz * fga * pb_yzzz

                     + 3.0 * fx * fx * fz * pb_yzzz

                     - fz * fga * pa_zz * pb_yzzz

                     + 7.0 * pa_xx * fz * fx * pb_yzzz

                     + 7.0 * fx * pa_zz * fz * pb_yzzz

                     + 16.0 * pa_xxzz * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_zzzz_s_0(double fx,
                                     double pa_xx,
                                     double pa_xxz,
                                     double pa_xxzz,
                                     double pa_z,
                                     double pa_zz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double pb_zzzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.9375 * fx * fx * fx * fx

                     + 1.875 * pa_xx * fx * fx * fx

                     + 0.375 * fx * fx * fx * pa_zz

                     + 3.0 * fx * fx * fx * pa_z * pb_z

                     + 2.25 * fx * fx * fx * pb_zz

                     + 0.75 * pa_xxzz * fx * fx

                     + 6.0 * pa_xxz * fx * fx * pb_z

                     + 4.5 * pa_xx * fx * fx * pb_zz

                     + 1.5 * fx * fx * pa_zz * pb_zz

                     + 2.0 * fx * fx * pa_z * pb_zzz

                     + 3.0 * pa_xxzz * pb_zz * fx

                     + 4.0 * pa_xxz * fx * pb_zzz

                     + 0.25 * fx * fx * pb_zzzz

                     + 0.5 * pa_xx * fx * pb_zzzz

                     + 0.5 * fx * pa_zz * pb_zzzz

                     + pa_xxzz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_zzzz_r_0(double fga,
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
                                     double pb_zzzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * fx * fz * fgb

                     - 2.25 * fz * fga * fx * fx * fx

                     - 4.5 * pa_xx * fx * fx * fz * fgb

                     + 7.5 * fx * fx * fx * fx * fz

                     + 18.75 * pa_xx * fx * fx * fx * fz

                     - 0.75 * pa_xx * fz * fga * fx * fx

                     - 1.5 * fx * fx * pa_zz * fz * fgb

                     - 6.0 * fx * fx * pa_z * pb_z * fz * fgb

                     - 0.75 * fz * fga * pa_zz * fx * fx

                     - 6.0 * fz * fga * pa_z * fx * fx * pb_z

                     - 6.0 * fz * fga * fx * fx * pb_zz

                     - 3.0 * pa_xxzz * fx * fz * fgb

                     - 12.0 * pa_xxz * fx * pb_z * fz * fgb

                     + 3.75 * fx * fx * fx * pa_zz * fz

                     + 30.0 * fx * fx * fx * pa_z * fz * pb_z

                     + 22.5 * fx * fx * fx * fz * pb_zz

                     + 9.0 * pa_xxzz * fz * fx * fx

                     + 72.0 * pa_xxz * fz * fx * fx * pb_z

                     + 54.0 * pa_xx * fx * fx * fz * pb_zz

                     - 1.5 * fx * fx * pb_zz * fz * fgb

                     - 3.0 * pa_xx * fx * pb_zz * fz * fgb

                     - 3.0 * pa_xx * fz * fga * pb_zz * fx

                     - 3.0 * fx * pa_zz * pb_zz * fz * fgb

                     - 3.0 * fz * fga * pa_zz * pb_zz * fx

                     - 4.0 * fz * fga * pa_z * fx * pb_zzz

                     - 6.0 * pa_xxzz * pb_zz * fz * fgb

                     + 18.0 * fx * fx * pa_zz * fz * pb_zz

                     + 24.0 * fx * fx * pa_z * fz * pb_zzz

                     + 42.0 * pa_xxzz * fz * pb_zz * fx

                     + 56.0 * pa_xxz * fz * fx * pb_zzz

                     - fx * fz * fga * pb_zzzz

                     - pa_xx * fz * fga * pb_zzzz

                     + 3.0 * fx * fx * fz * pb_zzzz

                     - fz * fga * pa_zz * pb_zzzz

                     + 7.0 * pa_xx * fz * fx * pb_zzzz

                     + 7.0 * fx * pa_zz * fz * pb_zzzz

                     + 16.0 * pa_xxzz * fz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xxxx_s_0(double fx,
                                     double pa_xy,
                                     double pa_xyyy,
                                     double pa_y,
                                     double pa_yyy,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxx,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_xy * fx * fx * fx

                     + 4.5 * fx * fx * fx * pa_y * pb_x

                     + 0.75 * pa_xyyy * fx * fx

                     + 3.0 * fx * fx * pa_yyy * pb_x

                     + 4.5 * pa_xy * fx * fx * pb_xx

                     + 3.0 * fx * fx * pa_y * pb_xxx

                     + 3.0 * pa_xyyy * pb_xx * fx

                     + 2.0 * fx * pa_yyy * pb_xxx

                     + 1.5 * pa_xy * fx * pb_xxxx

                     + pa_xyyy * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xxxx_r_0(double fga,
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
                                     double pb_xxxx,
                                     double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_xy * fx * fx * fz * fgb

                     - 2.25 * pa_xy * fz * fga * fx * fx

                     - 9.0 * fx * fx * pa_y * pb_x * fz * fgb

                     - 9.0 * fx * fx * pa_y * fz * fga * pb_x

                     - 3.0 * pa_xyyy * fx * fz * fgb

                     + 11.25 * pa_xy * fz * fx * fx * fx

                     - 6.0 * fx * pa_yyy * pb_x * fz * fgb

                     + 45.0 * fx * fx * fx * pa_y * fz * pb_x

                     + 9.0 * pa_xyyy * fz * fx * fx

                     + 36.0 * fx * fx * pa_yyy * fz * pb_x

                     - 9.0 * pa_xy * fx * pb_xx * fz * fgb

                     - 9.0 * pa_xy * fz * fga * pb_xx * fx

                     - 6.0 * fx * pa_y * fz * fga * pb_xxx

                     - 6.0 * pa_xyyy * pb_xx * fz * fgb

                     + 54.0 * pa_xy * fz * fx * fx * pb_xx

                     + 36.0 * fx * fx * pa_y * fz * pb_xxx

                     + 42.0 * pa_xyyy * fz * pb_xx * fx

                     + 28.0 * fx * pa_yyy * fz * pb_xxx

                     - 3.0 * pa_xy * fz * fga * pb_xxxx

                     + 21.0 * pa_xy * fz * fx * pb_xxxx

                     + 16.0 * pa_xyyy * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xxxy_s_0(double fx,
                                     double pa_x,
                                     double pa_xy,
                                     double pa_xyy,
                                     double pa_xyyy,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxy,
                                     double pb_xxy,
                                     double pb_xy,
                                     double pb_y,
                                     double s_0_0)
    {
        return s_0_0 * (0.5625 * fx * fx * fx * fx

                     + 1.125 * fx * fx * fx * pa_yy

                     + 1.125 * pa_x * fx * fx * fx * pb_x

                     + 1.125 * fx * fx * fx * pa_y * pb_y

                     + 1.125 * fx * fx * fx * pb_xx

                     + 2.25 * pa_xyy * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_yyy * pb_y

                     + 2.25 * fx * fx * pa_yy * pb_xx

                     + 2.25 * pa_xy * fx * fx * pb_xy

                     + 0.75 * pa_x * fx * fx * pb_xxx

                     + 2.25 * fx * fx * pa_y * pb_xxy

                     + 1.5 * pa_xyyy * pb_xy * fx

                     + 1.5 * pa_xyy * fx * pb_xxx

                     + 1.5 * fx * pa_yyy * pb_xxy

                     + 1.5 * pa_xy * fx * pb_xxxy

                     + pa_xyyy * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xxxy_r_0(double fga,
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
                                     double pb_xxx,
                                     double pb_xxxy,
                                     double pb_xxy,
                                     double pb_xy,
                                     double pb_y,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.125 * fx * fx * fx * fz * fgb

                     - 1.125 * fx * fx * fx * fz * fga

                     - 2.25 * fx * fx * pa_yy * fz * fgb

                     + 4.5 * fx * fx * fx * fx * fz

                     + 11.25 * fx * fx * fx * pa_yy * fz

                     - 2.25 * pa_x * fx * fx * pb_x * fz * fgb

                     - 2.25 * pa_x * fx * fx * fz * fga * pb_x

                     - 2.25 * fx * fx * pa_y * fz * fgb * pb_y

                     - 2.25 * fx * fx * pa_y * fz * fga * pb_y

                     - 2.25 * fx * fx * fz * fga * pb_xx

                     - 4.5 * pa_xyy * fx * pb_x * fz * fgb

                     + 11.25 * pa_x * fx * fx * fx * fz * pb_x

                     - 1.5 * fx * pa_yyy * fz * fgb * pb_y

                     + 11.25 * fx * fx * fx * pa_y * fz * pb_y

                     + 11.25 * fx * fx * fx * fz * pb_xx

                     + 27.0 * pa_xyy * fz * fx * fx * pb_x

                     + 9.0 * fx * fx * pa_yyy * fz * pb_y

                     + 27.0 * fx * fx * pa_yy * fz * pb_xx

                     - 4.5 * pa_xy * fx * pb_xy * fz * fgb

                     - 4.5 * pa_xy * fz * fga * pb_xy * fx

                     - 1.5 * pa_x * fx * fz * fga * pb_xxx

                     - 4.5 * fx * pa_y * fz * fga * pb_xxy

                     - 3.0 * pa_xyyy * pb_xy * fz * fgb

                     + 27.0 * pa_xy * fz * fx * fx * pb_xy

                     + 9.0 * pa_x * fx * fx * fz * pb_xxx

                     + 27.0 * fx * fx * pa_y * fz * pb_xxy

                     + 21.0 * pa_xyyy * fz * pb_xy * fx

                     + 21.0 * pa_xyy * fz * fx * pb_xxx

                     + 21.0 * fx * pa_yyy * fz * pb_xxy

                     - 3.0 * pa_xy * fz * fga * pb_xxxy

                     + 21.0 * pa_xy * fz * fx * pb_xxxy

                     + 16.0 * pa_xyyy * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xxxz_s_0(double fx,
                                     double pa_xy,
                                     double pa_xyyy,
                                     double pa_y,
                                     double pa_yyy,
                                     double pb_xxxz,
                                     double pb_xxz,
                                     double pb_xz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pa_y * pb_z

                     + 0.75 * fx * fx * pa_yyy * pb_z

                     + 2.25 * pa_xy * fx * fx * pb_xz

                     + 2.25 * fx * fx * pa_y * pb_xxz

                     + 1.5 * pa_xyyy * pb_xz * fx

                     + 1.5 * fx * pa_yyy * pb_xxz

                     + 1.5 * pa_xy * fx * pb_xxxz

                     + pa_xyyy * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xxxz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xy,
                                     double pa_xyyy,
                                     double pa_y,
                                     double pa_yyy,
                                     double pb_xxxz,
                                     double pb_xxz,
                                     double pb_xz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pa_y * fz * fgb * pb_z

                     - 2.25 * fx * fx * pa_y * fz * fga * pb_z

                     - 1.5 * fx * pa_yyy * fz * fgb * pb_z

                     + 11.25 * fx * fx * fx * pa_y * fz * pb_z

                     + 9.0 * fx * fx * pa_yyy * fz * pb_z

                     - 4.5 * pa_xy * fx * pb_xz * fz * fgb

                     - 4.5 * pa_xy * fz * fga * pb_xz * fx

                     - 4.5 * fx * pa_y * fz * fga * pb_xxz

                     - 3.0 * pa_xyyy * pb_xz * fz * fgb

                     + 27.0 * pa_xy * fz * fx * fx * pb_xz

                     + 27.0 * fx * fx * pa_y * fz * pb_xxz

                     + 21.0 * pa_xyyy * fz * pb_xz * fx

                     + 21.0 * fx * pa_yyy * fz * pb_xxz

                     - 3.0 * pa_xy * fz * fga * pb_xxxz

                     + 21.0 * pa_xy * fz * fx * pb_xxxz

                     + 16.0 * pa_xyyy * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xxyy_s_0(double fx,
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
                                     double pb_xxyy,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_y,
                                     double pb_yy,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_xy * fx * fx * fx

                     + 2.25 * fx * fx * fx * pa_y * pb_x

                     + 0.75 * pa_x * fx * fx * fx * pb_y

                     + 1.5 * fx * fx * fx * pb_xy

                     + 0.25 * pa_xyyy * fx * fx

                     + 1.5 * pa_xyy * fx * fx * pb_y

                     + 2.25 * pa_xy * fx * fx * pb_xx

                     + 0.5 * fx * fx * pa_yyy * pb_x

                     + 3.0 * fx * fx * pa_yy * pb_xy

                     + 0.75 * pa_xy * fx * fx * pb_yy

                     + 1.5 * pa_x * fx * fx * pb_xxy

                     + 1.5 * fx * fx * pa_y * pb_xyy

                     + 0.5 * pa_xyyy * pb_xx * fx

                     + 0.5 * pa_xyyy * fx * pb_yy

                     + 3.0 * pa_xyy * fx * pb_xxy

                     + fx * pa_yyy * pb_xyy

                     + 1.5 * pa_xy * fx * pb_xxyy

                     + pa_xyyy * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xxyy_r_0(double fga,
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
                                     double pb_xxyy,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_y,
                                     double pb_yy,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xy * fx * fx * fz * fgb

                     + 11.25 * pa_xy * fx * fx * fx * fz

                     + 22.5 * fx * fx * fx * pa_y * fz * pb_x

                     - 0.75 * pa_xy * fz * fga * fx * fx

                     - 1.5 * pa_x * fx * fx * fz * fgb * pb_y

                     - 1.5 * pa_x * fx * fx * fz * fga * pb_y

                     - 1.5 * fx * fx * pa_y * pb_x * fz * fgb

                     - 1.5 * fx * fx * pa_y * fz * fga * pb_x

                     - 3.0 * fx * fx * fz * fga * pb_xy

                     - pa_xyyy * fx * fz * fgb

                     - 3.0 * pa_xyy * fx * fz * fgb * pb_y

                     + 7.5 * pa_x * fx * fx * fx * fz * pb_y

                     - fx * pa_yyy * pb_x * fz * fgb

                     + 15.0 * fx * fx * fx * fz * pb_xy

                     + 3.0 * pa_xyyy * fz * fx * fx

                     + 18.0 * pa_xyy * fz * fx * fx * pb_y

                     + 27.0 * pa_xy * fx * fx * fz * pb_xx

                     + 6.0 * fx * fx * pa_yyy * fz * pb_x

                     + 36.0 * fx * fx * pa_yy * fz * pb_xy

                     - 1.5 * pa_xy * fx * pb_xx * fz * fgb

                     - 1.5 * pa_xy * fx * fz * fgb * pb_yy

                     - 1.5 * pa_xy * fz * fga * pb_xx * fx

                     - 1.5 * pa_xy * fz * fga * fx * pb_yy

                     - 3.0 * pa_x * fx * fz * fga * pb_xxy

                     - 3.0 * fx * pa_y * fz * fga * pb_xyy

                     - pa_xyyy * pb_xx * fz * fgb

                     - pa_xyyy * fz * fgb * pb_yy

                     + 9.0 * pa_xy * fz * fx * fx * pb_yy

                     + 18.0 * pa_x * fx * fx * fz * pb_xxy

                     + 18.0 * fx * fx * pa_y * fz * pb_xyy

                     + 7.0 * pa_xyyy * fz * pb_xx * fx

                     + 7.0 * pa_xyyy * fz * fx * pb_yy

                     + 42.0 * pa_xyy * fz * fx * pb_xxy

                     + 14.0 * fx * pa_yyy * fz * pb_xyy

                     - 3.0 * pa_xy * fz * fga * pb_xxyy

                     + 21.0 * pa_xy * fz * fx * pb_xxyy

                     + 16.0 * pa_xyyy * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xxyz_s_0(double fx,
                                     double pa_x,
                                     double pa_xy,
                                     double pa_xyy,
                                     double pa_xyyy,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pb_xxyz,
                                     double pb_xxz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * pa_x * fx * fx * fx * pb_z

                     + 0.75 * fx * fx * fx * pb_xz

                     + 0.75 * pa_xyy * fx * fx * pb_z

                     + 1.5 * fx * fx * pa_yy * pb_xz

                     + 0.75 * pa_xy * fx * fx * pb_yz

                     + 0.75 * pa_x * fx * fx * pb_xxz

                     + 1.5 * fx * fx * pa_y * pb_xyz

                     + 0.5 * pa_xyyy * fx * pb_yz

                     + 1.5 * pa_xyy * fx * pb_xxz

                     + fx * pa_yyy * pb_xyz

                     + 1.5 * pa_xy * fx * pb_xxyz

                     + pa_xyyy * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xxyz_r_0(double fga,
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
                                     double pb_xxyz,
                                     double pb_xxz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.75 * pa_x * fx * fx * fz * fgb * pb_z

                     - 0.75 * pa_x * fx * fx * fz * fga * pb_z

                     - 1.5 * fx * fx * fz * fga * pb_xz

                     - 1.5 * pa_xyy * fx * fz * fgb * pb_z

                     + 3.75 * pa_x * fx * fx * fx * fz * pb_z

                     + 7.5 * fx * fx * fx * fz * pb_xz

                     + 9.0 * pa_xyy * fz * fx * fx * pb_z

                     + 18.0 * fx * fx * pa_yy * fz * pb_xz

                     - 1.5 * pa_xy * fx * fz * fgb * pb_yz

                     - 1.5 * pa_xy * fz * fga * fx * pb_yz

                     - 1.5 * pa_x * fx * fz * fga * pb_xxz

                     - 3.0 * fx * pa_y * fz * fga * pb_xyz

                     - pa_xyyy * fz * fgb * pb_yz

                     + 9.0 * pa_xy * fz * fx * fx * pb_yz

                     + 9.0 * pa_x * fx * fx * fz * pb_xxz

                     + 18.0 * fx * fx * pa_y * fz * pb_xyz

                     + 7.0 * pa_xyyy * fz * fx * pb_yz

                     + 21.0 * pa_xyy * fz * fx * pb_xxz

                     + 14.0 * fx * pa_yyy * fz * pb_xyz

                     - 3.0 * pa_xy * fz * fga * pb_xxyz

                     + 21.0 * pa_xy * fz * fx * pb_xxyz

                     + 16.0 * pa_xyyy * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xxzz_s_0(double fx,
                                     double pa_xy,
                                     double pa_xyyy,
                                     double pa_y,
                                     double pa_yyy,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxzz,
                                     double pb_xzz,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * pa_xy * fx * fx * fx

                     + 0.75 * fx * fx * fx * pa_y * pb_x

                     + 0.25 * pa_xyyy * fx * fx

                     + 0.5 * fx * fx * pa_yyy * pb_x

                     + 0.75 * pa_xy * fx * fx * pb_xx

                     + 0.75 * pa_xy * fx * fx * pb_zz

                     + 1.5 * fx * fx * pa_y * pb_xzz

                     + 0.5 * pa_xyyy * pb_xx * fx

                     + 0.5 * pa_xyyy * fx * pb_zz

                     + fx * pa_yyy * pb_xzz

                     + 1.5 * pa_xy * fx * pb_xxzz

                     + pa_xyyy * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xxzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xy,
                                     double pa_xyyy,
                                     double pa_y,
                                     double pa_yyy,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxzz,
                                     double pb_xzz,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xy * fx * fx * fz * fgb

                     - 0.75 * pa_xy * fz * fga * fx * fx

                     - 1.5 * fx * fx * pa_y * pb_x * fz * fgb

                     - 1.5 * fx * fx * pa_y * fz * fga * pb_x

                     - pa_xyyy * fx * fz * fgb

                     + 3.75 * pa_xy * fz * fx * fx * fx

                     - fx * pa_yyy * pb_x * fz * fgb

                     + 7.5 * fx * fx * fx * pa_y * fz * pb_x

                     + 3.0 * pa_xyyy * fz * fx * fx

                     + 6.0 * fx * fx * pa_yyy * fz * pb_x

                     - 1.5 * pa_xy * fx * pb_xx * fz * fgb

                     - 1.5 * pa_xy * fx * fz * fgb * pb_zz

                     - 1.5 * pa_xy * fz * fga * pb_xx * fx

                     - 1.5 * pa_xy * fz * fga * fx * pb_zz

                     - 3.0 * fx * pa_y * fz * fga * pb_xzz

                     - pa_xyyy * pb_xx * fz * fgb

                     - pa_xyyy * fz * fgb * pb_zz

                     + 9.0 * pa_xy * fz * fx * fx * pb_xx

                     + 9.0 * pa_xy * fz * fx * fx * pb_zz

                     + 18.0 * fx * fx * pa_y * fz * pb_xzz

                     + 7.0 * pa_xyyy * fz * pb_xx * fx

                     + 7.0 * pa_xyyy * fz * fx * pb_zz

                     + 14.0 * fx * pa_yyy * fz * pb_xzz

                     - 3.0 * pa_xy * fz * fga * pb_xxzz

                     + 21.0 * pa_xy * fz * fx * pb_xxzz

                     + 16.0 * pa_xyyy * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xyyy_s_0(double fx,
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
                                     double pb_xyyy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double s_0_0)
    {
        return s_0_0 * (0.9375 * fx * fx * fx * fx

                     + 1.875 * pa_x * fx * fx * fx * pb_x

                     + 1.125 * fx * fx * fx * pa_yy

                     + 3.375 * fx * fx * fx * pa_y * pb_y

                     + 1.125 * fx * fx * fx * pb_yy

                     + 2.25 * pa_xyy * fx * fx * pb_x

                     + 6.75 * pa_xy * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_yyy * pb_y

                     + 2.25 * fx * fx * pa_yy * pb_yy

                     + 2.25 * pa_x * fx * fx * pb_xyy

                     + 0.75 * fx * fx * pa_y * pb_yyy

                     + 1.5 * pa_xyyy * pb_xy * fx

                     + 4.5 * pa_xyy * fx * pb_xyy

                     + 0.5 * fx * pa_yyy * pb_yyy

                     + 1.5 * pa_xy * fx * pb_xyyy

                     + pa_xyyy * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xyyy_r_0(double fga,
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
                                     double pb_xyyy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double r_0_0)
    {
        return r_0_0 * (7.5 * fx * fx * fx * fx * fz

                     - 1.125 * fx * fx * fx * fz * fgb

                     - 1.125 * fx * fx * fx * fz * fga

                     - 2.25 * fx * fx * pa_yy * fz * fgb

                     + 18.75 * pa_x * fx * fx * fx * fz * pb_x

                     + 11.25 * fx * fx * fx * pa_yy * fz

                     + 33.75 * fx * fx * fx * pa_y * fz * pb_y

                     - 2.25 * pa_x * fx * fx * pb_x * fz * fgb

                     - 2.25 * pa_x * fx * fx * fz * fga * pb_x

                     - 2.25 * fx * fx * pa_y * pb_y * fz * fgb

                     - 2.25 * fx * fx * pa_y * fz * fga * pb_y

                     - 2.25 * fx * fx * fz * fga * pb_yy

                     - 4.5 * pa_xyy * fx * pb_x * fz * fgb

                     - 1.5 * fx * pa_yyy * pb_y * fz * fgb

                     + 11.25 * fx * fx * fx * fz * pb_yy

                     + 27.0 * pa_xyy * fz * fx * fx * pb_x

                     + 81.0 * pa_xy * fx * fx * fz * pb_xy

                     + 9.0 * fx * fx * pa_yyy * fz * pb_y

                     + 27.0 * fx * fx * pa_yy * fz * pb_yy

                     - 4.5 * pa_xy * fx * pb_xy * fz * fgb

                     - 4.5 * pa_xy * fz * fga * pb_xy * fx

                     - 4.5 * pa_x * fx * fz * fga * pb_xyy

                     - 1.5 * fx * pa_y * fz * fga * pb_yyy

                     - 3.0 * pa_xyyy * pb_xy * fz * fgb

                     + 27.0 * pa_x * fx * fx * fz * pb_xyy

                     + 9.0 * fx * fx * pa_y * fz * pb_yyy

                     + 21.0 * pa_xyyy * fz * pb_xy * fx

                     + 63.0 * pa_xyy * fz * fx * pb_xyy

                     + 7.0 * fx * pa_yyy * fz * pb_yyy

                     - 3.0 * pa_xy * fz * fga * pb_xyyy

                     + 21.0 * pa_xy * fz * fx * pb_xyyy

                     + 16.0 * pa_xyyy * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xyyz_s_0(double fx,
                                     double pa_x,
                                     double pa_xy,
                                     double pa_xyy,
                                     double pa_xyyy,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pb_xyyz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pa_y * pb_z

                     + 0.75 * fx * fx * fx * pb_yz

                     + 2.25 * pa_xy * fx * fx * pb_xz

                     + 0.25 * fx * fx * pa_yyy * pb_z

                     + 1.5 * fx * fx * pa_yy * pb_yz

                     + 1.5 * pa_x * fx * fx * pb_xyz

                     + 0.75 * fx * fx * pa_y * pb_yyz

                     + 0.5 * pa_xyyy * pb_xz * fx

                     + 3.0 * pa_xyy * fx * pb_xyz

                     + 0.5 * fx * pa_yyy * pb_yyz

                     + 1.5 * pa_xy * fx * pb_xyyz

                     + pa_xyyy * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xyyz_r_0(double fga,
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
                                     double pb_xyyz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (11.25 * fx * fx * fx * pa_y * fz * pb_z

                     - 0.75 * fx * fx * pa_y * fz * fgb * pb_z

                     - 0.75 * fx * fx * pa_y * fz * fga * pb_z

                     - 1.5 * fx * fx * fz * fga * pb_yz

                     - 0.5 * fx * pa_yyy * fz * fgb * pb_z

                     + 7.5 * fx * fx * fx * fz * pb_yz

                     + 27.0 * pa_xy * fx * fx * fz * pb_xz

                     + 3.0 * fx * fx * pa_yyy * fz * pb_z

                     + 18.0 * fx * fx * pa_yy * fz * pb_yz

                     - 1.5 * pa_xy * fx * pb_xz * fz * fgb

                     - 1.5 * pa_xy * fz * fga * pb_xz * fx

                     - 3.0 * pa_x * fx * fz * fga * pb_xyz

                     - 1.5 * fx * pa_y * fz * fga * pb_yyz

                     - pa_xyyy * pb_xz * fz * fgb

                     + 18.0 * pa_x * fx * fx * fz * pb_xyz

                     + 9.0 * fx * fx * pa_y * fz * pb_yyz

                     + 7.0 * pa_xyyy * fz * pb_xz * fx

                     + 42.0 * pa_xyy * fz * fx * pb_xyz

                     + 7.0 * fx * pa_yyy * fz * pb_yyz

                     - 3.0 * pa_xy * fz * fga * pb_xyyz

                     + 21.0 * pa_xy * fz * fx * pb_xyyz

                     + 16.0 * pa_xyyy * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xyzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xy,
                                     double pa_xyy,
                                     double pa_xyyy,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pb_x,
                                     double pb_xy,
                                     double pb_xyzz,
                                     double pb_xzz,
                                     double pb_y,
                                     double pb_yzz,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * fx * fx * fx * pa_yy

                     + 0.375 * pa_x * fx * fx * fx * pb_x

                     + 0.375 * fx * fx * fx * pa_y * pb_y

                     + 0.375 * fx * fx * fx * pb_zz

                     + 0.75 * pa_xyy * fx * fx * pb_x

                     + 0.25 * fx * fx * pa_yyy * pb_y

                     + 0.75 * fx * fx * pa_yy * pb_zz

                     + 0.75 * pa_xy * fx * fx * pb_xy

                     + 0.75 * pa_x * fx * fx * pb_xzz

                     + 0.75 * fx * fx * pa_y * pb_yzz

                     + 0.5 * pa_xyyy * pb_xy * fx

                     + 1.5 * pa_xyy * fx * pb_xzz

                     + 0.5 * fx * pa_yyy * pb_yzz

                     + 1.5 * pa_xy * fx * pb_xyzz

                     + pa_xyyy * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xyzz_r_0(double fga,
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
                                     double pb_xyzz,
                                     double pb_xzz,
                                     double pb_y,
                                     double pb_yzz,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.375 * fx * fx * fx * fz * fgb

                     - 0.375 * fx * fx * fx * fz * fga

                     - 0.75 * fx * fx * pa_yy * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     + 3.75 * fx * fx * fx * pa_yy * fz

                     - 0.75 * pa_x * fx * fx * pb_x * fz * fgb

                     - 0.75 * pa_x * fx * fx * fz * fga * pb_x

                     - 0.75 * fx * fx * pa_y * pb_y * fz * fgb

                     - 0.75 * fx * fx * pa_y * fz * fga * pb_y

                     - 0.75 * fx * fx * fz * fga * pb_zz

                     - 1.5 * pa_xyy * fx * pb_x * fz * fgb

                     + 3.75 * pa_x * fx * fx * fx * fz * pb_x

                     - 0.5 * fx * pa_yyy * pb_y * fz * fgb

                     + 3.75 * fx * fx * fx * pa_y * fz * pb_y

                     + 3.75 * fx * fx * fx * fz * pb_zz

                     + 9.0 * pa_xyy * fz * fx * fx * pb_x

                     + 3.0 * fx * fx * pa_yyy * fz * pb_y

                     + 9.0 * fx * fx * pa_yy * fz * pb_zz

                     - 1.5 * pa_xy * fx * pb_xy * fz * fgb

                     - 1.5 * pa_xy * fz * fga * pb_xy * fx

                     - 1.5 * pa_x * fx * fz * fga * pb_xzz

                     - 1.5 * fx * pa_y * fz * fga * pb_yzz

                     - pa_xyyy * pb_xy * fz * fgb

                     + 9.0 * pa_xy * fz * fx * fx * pb_xy

                     + 9.0 * pa_x * fx * fx * fz * pb_xzz

                     + 9.0 * fx * fx * pa_y * fz * pb_yzz

                     + 7.0 * pa_xyyy * fz * pb_xy * fx

                     + 21.0 * pa_xyy * fz * fx * pb_xzz

                     + 7.0 * fx * pa_yyy * fz * pb_yzz

                     - 3.0 * pa_xy * fz * fga * pb_xyzz

                     + 21.0 * pa_xy * fz * fx * pb_xyzz

                     + 16.0 * pa_xyyy * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xzzz_s_0(double fx,
                                     double pa_xy,
                                     double pa_xyyy,
                                     double pa_y,
                                     double pa_yyy,
                                     double pb_xz,
                                     double pb_xzzz,
                                     double pb_z,
                                     double pb_zzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pa_y * pb_z

                     + 0.75 * fx * fx * pa_yyy * pb_z

                     + 2.25 * pa_xy * fx * fx * pb_xz

                     + 0.75 * fx * fx * pa_y * pb_zzz

                     + 1.5 * pa_xyyy * pb_xz * fx

                     + 0.5 * fx * pa_yyy * pb_zzz

                     + 1.5 * pa_xy * fx * pb_xzzz

                     + pa_xyyy * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xzzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xy,
                                     double pa_xyyy,
                                     double pa_y,
                                     double pa_yyy,
                                     double pb_xz,
                                     double pb_xzzz,
                                     double pb_z,
                                     double pb_zzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pa_y * pb_z * fz * fgb

                     - 2.25 * fx * fx * pa_y * fz * fga * pb_z

                     - 1.5 * fx * pa_yyy * pb_z * fz * fgb

                     + 11.25 * fx * fx * fx * pa_y * fz * pb_z

                     + 9.0 * fx * fx * pa_yyy * fz * pb_z

                     - 4.5 * pa_xy * fx * pb_xz * fz * fgb

                     - 4.5 * pa_xy * fz * fga * pb_xz * fx

                     - 1.5 * fx * pa_y * fz * fga * pb_zzz

                     - 3.0 * pa_xyyy * pb_xz * fz * fgb

                     + 27.0 * pa_xy * fz * fx * fx * pb_xz

                     + 9.0 * fx * fx * pa_y * fz * pb_zzz

                     + 21.0 * pa_xyyy * fz * pb_xz * fx

                     + 7.0 * fx * pa_yyy * fz * pb_zzz

                     - 3.0 * pa_xy * fz * fga * pb_xzzz

                     + 21.0 * pa_xy * fz * fx * pb_xzzz

                     + 16.0 * pa_xyyy * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_yyyy_s_0(double fx,
                                     double pa_x,
                                     double pa_xy,
                                     double pa_xyy,
                                     double pa_xyyy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double pb_yyyy,
                                     double s_0_0)
    {
        return s_0_0 * (5.625 * pa_xy * fx * fx * fx

                     + 7.5 * pa_x * fx * fx * fx * pb_y

                     + 0.75 * pa_xyyy * fx * fx

                     + 9.0 * pa_xyy * fx * fx * pb_y

                     + 13.5 * pa_xy * fx * fx * pb_yy

                     + 3.0 * pa_x * fx * fx * pb_yyy

                     + 3.0 * pa_xyyy * pb_yy * fx

                     + 6.0 * pa_xyy * fx * pb_yyy

                     + 1.5 * pa_xy * fx * pb_yyyy

                     + pa_xyyy * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_yyyy_r_0(double fga,
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
                                     double pb_yyyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 13.5 * pa_xy * fx * fx * fz * fgb

                     + 56.25 * pa_xy * fx * fx * fx * fz

                     + 75.0 * pa_x * fx * fx * fx * fz * pb_y

                     - 2.25 * pa_xy * fz * fga * fx * fx

                     - 9.0 * pa_x * fx * fx * pb_y * fz * fgb

                     - 9.0 * pa_x * fx * fx * fz * fga * pb_y

                     - 3.0 * pa_xyyy * fx * fz * fgb

                     - 18.0 * pa_xyy * fx * pb_y * fz * fgb

                     + 9.0 * pa_xyyy * fz * fx * fx

                     + 108.0 * pa_xyy * fz * fx * fx * pb_y

                     + 162.0 * pa_xy * fx * fx * fz * pb_yy

                     - 9.0 * pa_xy * fx * pb_yy * fz * fgb

                     - 9.0 * pa_xy * fz * fga * pb_yy * fx

                     - 6.0 * pa_x * fx * fz * fga * pb_yyy

                     - 6.0 * pa_xyyy * pb_yy * fz * fgb

                     + 36.0 * pa_x * fx * fx * fz * pb_yyy

                     + 42.0 * pa_xyyy * fz * pb_yy * fx

                     + 84.0 * pa_xyy * fz * fx * pb_yyy

                     - 3.0 * pa_xy * fz * fga * pb_yyyy

                     + 21.0 * pa_xy * fz * fx * pb_yyyy

                     + 16.0 * pa_xyyy * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_yyyz_s_0(double fx,
                                     double pa_x,
                                     double pa_xy,
                                     double pa_xyy,
                                     double pa_xyyy,
                                     double pb_yyyz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (1.875 * pa_x * fx * fx * fx * pb_z

                     + 2.25 * pa_xyy * fx * fx * pb_z

                     + 6.75 * pa_xy * fx * fx * pb_yz

                     + 2.25 * pa_x * fx * fx * pb_yyz

                     + 1.5 * pa_xyyy * pb_yz * fx

                     + 4.5 * pa_xyy * fx * pb_yyz

                     + 1.5 * pa_xy * fx * pb_yyyz

                     + pa_xyyy * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_yyyz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xy,
                                     double pa_xyy,
                                     double pa_xyyy,
                                     double pb_yyyz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (18.75 * pa_x * fx * fx * fx * fz * pb_z

                     - 2.25 * pa_x * fx * fx * fz * fgb * pb_z

                     - 2.25 * pa_x * fx * fx * fz * fga * pb_z

                     - 4.5 * pa_xyy * fx * fz * fgb * pb_z

                     + 27.0 * pa_xyy * fz * fx * fx * pb_z

                     + 81.0 * pa_xy * fx * fx * fz * pb_yz

                     - 4.5 * pa_xy * fx * pb_yz * fz * fgb

                     - 4.5 * pa_xy * fz * fga * pb_yz * fx

                     - 4.5 * pa_x * fx * fz * fga * pb_yyz

                     - 3.0 * pa_xyyy * pb_yz * fz * fgb

                     + 27.0 * pa_x * fx * fx * fz * pb_yyz

                     + 21.0 * pa_xyyy * fz * pb_yz * fx

                     + 63.0 * pa_xyy * fz * fx * pb_yyz

                     - 3.0 * pa_xy * fz * fga * pb_yyyz

                     + 21.0 * pa_xy * fz * fx * pb_yyyz

                     + 16.0 * pa_xyyy * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_yyzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xy,
                                     double pa_xyy,
                                     double pa_xyyy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyzz,
                                     double pb_yzz,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_xy * fx * fx * fx

                     + 0.75 * pa_x * fx * fx * fx * pb_y

                     + 0.25 * pa_xyyy * fx * fx

                     + 1.5 * pa_xyy * fx * fx * pb_y

                     + 2.25 * pa_xy * fx * fx * pb_zz

                     + 0.75 * pa_xy * fx * fx * pb_yy

                     + 1.5 * pa_x * fx * fx * pb_yzz

                     + 0.5 * pa_xyyy * pb_yy * fx

                     + 0.5 * pa_xyyy * fx * pb_zz

                     + 3.0 * pa_xyy * fx * pb_yzz

                     + 1.5 * pa_xy * fx * pb_yyzz

                     + pa_xyyy * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_yyzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xy,
                                     double pa_xyy,
                                     double pa_xyyy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyzz,
                                     double pb_yzz,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xy * fx * fx * fz * fgb

                     + 11.25 * pa_xy * fx * fx * fx * fz

                     - 0.75 * pa_xy * fz * fga * fx * fx

                     - 1.5 * pa_x * fx * fx * pb_y * fz * fgb

                     - 1.5 * pa_x * fx * fx * fz * fga * pb_y

                     - pa_xyyy * fx * fz * fgb

                     - 3.0 * pa_xyy * fx * pb_y * fz * fgb

                     + 7.5 * pa_x * fx * fx * fx * fz * pb_y

                     + 3.0 * pa_xyyy * fz * fx * fx

                     + 18.0 * pa_xyy * fz * fx * fx * pb_y

                     + 27.0 * pa_xy * fx * fx * fz * pb_zz

                     - 1.5 * pa_xy * fx * pb_yy * fz * fgb

                     - 1.5 * pa_xy * fx * fz * fgb * pb_zz

                     - 1.5 * pa_xy * fz * fga * pb_yy * fx

                     - 1.5 * pa_xy * fz * fga * fx * pb_zz

                     - 3.0 * pa_x * fx * fz * fga * pb_yzz

                     - pa_xyyy * pb_yy * fz * fgb

                     - pa_xyyy * fz * fgb * pb_zz

                     + 9.0 * pa_xy * fz * fx * fx * pb_yy

                     + 18.0 * pa_x * fx * fx * fz * pb_yzz

                     + 7.0 * pa_xyyy * fz * pb_yy * fx

                     + 7.0 * pa_xyyy * fz * fx * pb_zz

                     + 42.0 * pa_xyy * fz * fx * pb_yzz

                     - 3.0 * pa_xy * fz * fga * pb_yyzz

                     + 21.0 * pa_xy * fz * fx * pb_yyzz

                     + 16.0 * pa_xyyy * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_yzzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xy,
                                     double pa_xyy,
                                     double pa_xyyy,
                                     double pb_yz,
                                     double pb_yzzz,
                                     double pb_z,
                                     double pb_zzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_x * fx * fx * fx * pb_z

                     + 2.25 * pa_xyy * fx * fx * pb_z

                     + 2.25 * pa_xy * fx * fx * pb_yz

                     + 0.75 * pa_x * fx * fx * pb_zzz

                     + 1.5 * pa_xyyy * pb_yz * fx

                     + 1.5 * pa_xyy * fx * pb_zzz

                     + 1.5 * pa_xy * fx * pb_yzzz

                     + pa_xyyy * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_yzzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xy,
                                     double pa_xyy,
                                     double pa_xyyy,
                                     double pb_yz,
                                     double pb_yzzz,
                                     double pb_z,
                                     double pb_zzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * pa_x * fx * fx * pb_z * fz * fgb

                     - 2.25 * pa_x * fx * fx * fz * fga * pb_z

                     - 4.5 * pa_xyy * fx * pb_z * fz * fgb

                     + 11.25 * pa_x * fx * fx * fx * fz * pb_z

                     + 27.0 * pa_xyy * fz * fx * fx * pb_z

                     - 4.5 * pa_xy * fx * pb_yz * fz * fgb

                     - 4.5 * pa_xy * fz * fga * pb_yz * fx

                     - 1.5 * pa_x * fx * fz * fga * pb_zzz

                     - 3.0 * pa_xyyy * pb_yz * fz * fgb

                     + 27.0 * pa_xy * fz * fx * fx * pb_yz

                     + 9.0 * pa_x * fx * fx * fz * pb_zzz

                     + 21.0 * pa_xyyy * fz * pb_yz * fx

                     + 21.0 * pa_xyy * fz * fx * pb_zzz

                     - 3.0 * pa_xy * fz * fga * pb_yzzz

                     + 21.0 * pa_xy * fz * fx * pb_yzzz

                     + 16.0 * pa_xyyy * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_zzzz_s_0(double fx,
                                     double pa_xy,
                                     double pa_xyyy,
                                     double pb_zz,
                                     double pb_zzzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_xy * fx * fx * fx

                     + 0.75 * pa_xyyy * fx * fx

                     + 4.5 * pa_xy * fx * fx * pb_zz

                     + 3.0 * pa_xyyy * pb_zz * fx

                     + 1.5 * pa_xy * fx * pb_zzzz

                     + pa_xyyy * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_zzzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xy,
                                     double pa_xyyy,
                                     double pb_zz,
                                     double pb_zzzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_xy * fx * fx * fz * fgb

                     - 2.25 * pa_xy * fz * fga * fx * fx

                     - 3.0 * pa_xyyy * fx * fz * fgb

                     + 11.25 * pa_xy * fz * fx * fx * fx

                     + 9.0 * pa_xyyy * fz * fx * fx

                     - 9.0 * pa_xy * fx * pb_zz * fz * fgb

                     - 9.0 * pa_xy * fz * fga * pb_zz * fx

                     - 6.0 * pa_xyyy * pb_zz * fz * fgb

                     + 54.0 * pa_xy * fz * fx * fx * pb_zz

                     + 42.0 * pa_xyyy * fz * pb_zz * fx

                     - 3.0 * pa_xy * fz * fga * pb_zzzz

                     + 21.0 * pa_xy * fz * fx * pb_zzzz

                     + 16.0 * pa_xyyy * fz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xxxx_s_0(double fx,
                                     double pa_xyyz,
                                     double pa_xz,
                                     double pa_yyz,
                                     double pa_z,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxx,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * pa_xz * fx * fx * fx

                     + 1.5 * fx * fx * fx * pa_z * pb_x

                     + 0.75 * pa_xyyz * fx * fx

                     + 3.0 * fx * fx * pa_yyz * pb_x

                     + 1.5 * pa_xz * fx * fx * pb_xx

                     + fx * fx * pa_z * pb_xxx

                     + 3.0 * pa_xyyz * pb_xx * fx

                     + 2.0 * fx * pa_yyz * pb_xxx

                     + 0.5 * pa_xz * fx * pb_xxxx

                     + pa_xyyz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xxxx_r_0(double fga,
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
                                     double pb_xxxx,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xz * fx * fx * fz * fgb

                     - 0.75 * pa_xz * fz * fga * fx * fx

                     - 3.0 * fx * fx * pa_z * pb_x * fz * fgb

                     - 3.0 * fx * fx * fz * fga * pa_z * pb_x

                     - 3.0 * pa_xyyz * fx * fz * fgb

                     + 3.75 * pa_xz * fx * fx * fx * fz

                     - 6.0 * fx * pa_yyz * pb_x * fz * fgb

                     + 15.0 * fx * fx * fx * fz * pa_z * pb_x

                     + 9.0 * pa_xyyz * fz * fx * fx

                     + 36.0 * fx * fx * pa_yyz * fz * pb_x

                     - 3.0 * pa_xz * fx * pb_xx * fz * fgb

                     - 3.0 * pa_xz * fz * fga * pb_xx * fx

                     - 2.0 * fx * fz * fga * pa_z * pb_xxx

                     - 6.0 * pa_xyyz * pb_xx * fz * fgb

                     + 18.0 * pa_xz * fx * fx * fz * pb_xx

                     + 12.0 * fx * fx * fz * pa_z * pb_xxx

                     + 42.0 * pa_xyyz * fz * pb_xx * fx

                     + 28.0 * fx * pa_yyz * fz * pb_xxx

                     - pa_xz * fz * fga * pb_xxxx

                     + 7.0 * pa_xz * fx * fz * pb_xxxx

                     + 16.0 * pa_xyyz * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xxxy_s_0(double fx,
                                     double pa_xyyz,
                                     double pa_xyz,
                                     double pa_xz,
                                     double pa_yyz,
                                     double pa_yz,
                                     double pa_z,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxy,
                                     double pb_xxy,
                                     double pb_xy,
                                     double pb_y,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * fx * pa_yz

                     + 0.375 * fx * fx * fx * pa_z * pb_y

                     + 1.5 * pa_xyz * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_yyz * pb_y

                     + 1.5 * fx * fx * pa_yz * pb_xx

                     + 0.75 * pa_xz * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_z * pb_xxy

                     + 1.5 * pa_xyyz * pb_xy * fx

                     + pa_xyz * fx * pb_xxx

                     + 1.5 * fx * pa_yyz * pb_xxy

                     + 0.5 * pa_xz * fx * pb_xxxy

                     + pa_xyyz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xxxy_r_0(double fga,
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
                                     double pb_xxx,
                                     double pb_xxxy,
                                     double pb_xxy,
                                     double pb_xy,
                                     double pb_y,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_yz * fz * fgb

                     + 7.5 * fx * fx * fx * pa_yz * fz

                     - 0.75 * fx * fx * pa_z * fz * fgb * pb_y

                     - 0.75 * fx * fx * fz * fga * pa_z * pb_y

                     - 3.0 * pa_xyz * fx * pb_x * fz * fgb

                     - 1.5 * fx * pa_yyz * fz * fgb * pb_y

                     + 3.75 * fx * fx * fx * fz * pa_z * pb_y

                     + 18.0 * pa_xyz * fx * fx * fz * pb_x

                     + 9.0 * fx * fx * pa_yyz * fz * pb_y

                     + 18.0 * fx * fx * pa_yz * fz * pb_xx

                     - 1.5 * pa_xz * fx * pb_xy * fz * fgb

                     - 1.5 * pa_xz * fz * fga * pb_xy * fx

                     - 1.5 * fx * fz * fga * pa_z * pb_xxy

                     - 3.0 * pa_xyyz * pb_xy * fz * fgb

                     + 9.0 * pa_xz * fx * fx * fz * pb_xy

                     + 9.0 * fx * fx * fz * pa_z * pb_xxy

                     + 21.0 * pa_xyyz * fz * pb_xy * fx

                     + 14.0 * pa_xyz * fx * fz * pb_xxx

                     + 21.0 * fx * pa_yyz * fz * pb_xxy

                     - pa_xz * fz * fga * pb_xxxy

                     + 7.0 * pa_xz * fx * fz * pb_xxxy

                     + 16.0 * pa_xyyz * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xxxz_s_0(double fx,
                                     double pa_x,
                                     double pa_xyy,
                                     double pa_xyyz,
                                     double pa_xz,
                                     double pa_yy,
                                     double pa_yyz,
                                     double pa_z,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxz,
                                     double pb_xxz,
                                     double pb_xz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * fx * fx * fx * pa_yy

                     + 0.375 * pa_x * fx * fx * fx * pb_x

                     + 0.375 * fx * fx * fx * pa_z * pb_z

                     + 0.375 * fx * fx * fx * pb_xx

                     + 0.75 * pa_xyy * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_yyz * pb_z

                     + 0.75 * fx * fx * pa_yy * pb_xx

                     + 0.75 * pa_xz * fx * fx * pb_xz

                     + 0.25 * pa_x * fx * fx * pb_xxx

                     + 0.75 * fx * fx * pa_z * pb_xxz

                     + 1.5 * pa_xyyz * pb_xz * fx

                     + 0.5 * pa_xyy * fx * pb_xxx

                     + 1.5 * fx * pa_yyz * pb_xxz

                     + 0.5 * pa_xz * fx * pb_xxxz

                     + pa_xyyz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xxxz_r_0(double fga,
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
                                     double pb_xxx,
                                     double pb_xxxz,
                                     double pb_xxz,
                                     double pb_xz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.375 * fx * fx * fx * fz * fgb

                     - 0.375 * fx * fx * fx * fz * fga

                     - 0.75 * fx * fx * pa_yy * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     + 3.75 * fx * fx * fx * pa_yy * fz

                     - 0.75 * pa_x * fx * fx * pb_x * fz * fgb

                     - 0.75 * pa_x * fz * fga * fx * fx * pb_x

                     - 0.75 * fx * fx * pa_z * fz * fgb * pb_z

                     - 0.75 * fx * fx * fz * fga * pa_z * pb_z

                     - 0.75 * fx * fx * fz * fga * pb_xx

                     - 1.5 * pa_xyy * fx * pb_x * fz * fgb

                     + 3.75 * pa_x * fx * fx * fx * fz * pb_x

                     - 1.5 * fx * pa_yyz * fz * fgb * pb_z

                     + 3.75 * fx * fx * fx * fz * pa_z * pb_z

                     + 3.75 * fx * fx * fx * fz * pb_xx

                     + 9.0 * pa_xyy * fz * fx * fx * pb_x

                     + 9.0 * fx * fx * pa_yyz * fz * pb_z

                     + 9.0 * fx * fx * pa_yy * fz * pb_xx

                     - 1.5 * pa_xz * fx * pb_xz * fz * fgb

                     - 1.5 * pa_xz * fz * fga * pb_xz * fx

                     - 0.5 * pa_x * fz * fga * fx * pb_xxx

                     - 1.5 * fx * fz * fga * pa_z * pb_xxz

                     - 3.0 * pa_xyyz * pb_xz * fz * fgb

                     + 9.0 * pa_xz * fx * fx * fz * pb_xz

                     + 3.0 * pa_x * fx * fx * fz * pb_xxx

                     + 9.0 * fx * fx * fz * pa_z * pb_xxz

                     + 21.0 * pa_xyyz * fz * pb_xz * fx

                     + 7.0 * pa_xyy * fz * fx * pb_xxx

                     + 21.0 * fx * pa_yyz * fz * pb_xxz

                     - pa_xz * fz * fga * pb_xxxz

                     + 7.0 * pa_xz * fx * fz * pb_xxxz

                     + 16.0 * pa_xyyz * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xxyy_s_0(double fx,
                                     double pa_xyyz,
                                     double pa_xyz,
                                     double pa_xz,
                                     double pa_yyz,
                                     double pa_yz,
                                     double pa_z,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyy,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_y,
                                     double pb_yy,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * pa_xz * fx * fx * fx

                     + 0.75 * fx * fx * fx * pa_z * pb_x

                     + 0.25 * pa_xyyz * fx * fx

                     + pa_xyz * fx * fx * pb_y

                     + 0.75 * pa_xz * fx * fx * pb_xx

                     + 0.5 * fx * fx * pa_yyz * pb_x

                     + 2.0 * fx * fx * pa_yz * pb_xy

                     + 0.25 * pa_xz * fx * fx * pb_yy

                     + 0.5 * fx * fx * pa_z * pb_xyy

                     + 0.5 * pa_xyyz * pb_xx * fx

                     + 0.5 * pa_xyyz * fx * pb_yy

                     + 2.0 * pa_xyz * fx * pb_xxy

                     + fx * pa_yyz * pb_xyy

                     + 0.5 * pa_xz * fx * pb_xxyy

                     + pa_xyyz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xxyy_r_0(double fga,
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
                                     double pb_xxyy,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_y,
                                     double pb_yy,
                                     double r_0_0)
    {
        return r_0_0 * (- pa_xz * fx * fx * fz * fgb

                     + 3.75 * pa_xz * fx * fx * fx * fz

                     + 7.5 * fx * fx * fx * fz * pa_z * pb_x

                     - 0.25 * pa_xz * fz * fga * fx * fx

                     - 0.5 * fx * fx * pa_z * pb_x * fz * fgb

                     - 0.5 * fx * fx * fz * fga * pa_z * pb_x

                     - pa_xyyz * fx * fz * fgb

                     - 2.0 * pa_xyz * fx * fz * fgb * pb_y

                     - fx * pa_yyz * pb_x * fz * fgb

                     + 3.0 * pa_xyyz * fz * fx * fx

                     + 12.0 * pa_xyz * fx * fx * fz * pb_y

                     + 9.0 * pa_xz * fx * fx * fz * pb_xx

                     + 6.0 * fx * fx * pa_yyz * fz * pb_x

                     + 24.0 * fx * fx * pa_yz * fz * pb_xy

                     - 0.5 * pa_xz * fx * pb_xx * fz * fgb

                     - 0.5 * pa_xz * fx * fz * fgb * pb_yy

                     - 0.5 * pa_xz * fz * fga * pb_xx * fx

                     - 0.5 * pa_xz * fz * fga * fx * pb_yy

                     - fx * fz * fga * pa_z * pb_xyy

                     - pa_xyyz * pb_xx * fz * fgb

                     - pa_xyyz * fz * fgb * pb_yy

                     + 3.0 * pa_xz * fx * fx * fz * pb_yy

                     + 6.0 * fx * fx * fz * pa_z * pb_xyy

                     + 7.0 * pa_xyyz * fz * pb_xx * fx

                     + 7.0 * pa_xyyz * fz * fx * pb_yy

                     + 28.0 * pa_xyz * fx * fz * pb_xxy

                     + 14.0 * fx * pa_yyz * fz * pb_xyy

                     - pa_xz * fz * fga * pb_xxyy

                     + 7.0 * pa_xz * fx * fz * pb_xxyy

                     + 16.0 * pa_xyyz * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xxyz_s_0(double fx,
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
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyz,
                                     double pb_xxz,
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.25 * pa_xy * fx * fx * fx

                     + 0.5 * fx * fx * fx * pa_y * pb_x

                     + 0.125 * pa_x * fx * fx * fx * pb_y

                     + 0.25 * fx * fx * fx * pb_xy

                     + 0.25 * pa_xyy * fx * fx * pb_y

                     + 0.5 * pa_xyz * fx * fx * pb_z

                     + 0.5 * pa_xy * fx * fx * pb_xx

                     + 0.5 * fx * fx * pa_yy * pb_xy

                     + fx * fx * pa_yz * pb_xz

                     + 0.25 * pa_xz * fx * fx * pb_yz

                     + 0.25 * pa_x * fx * fx * pb_xxy

                     + 0.5 * fx * fx * pa_z * pb_xyz

                     + 0.5 * pa_xyyz * fx * pb_yz

                     + 0.5 * pa_xyy * fx * pb_xxy

                     + pa_xyz * fx * pb_xxz

                     + fx * pa_yyz * pb_xyz

                     + 0.5 * pa_xz * fx * pb_xxyz

                     + pa_xyyz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xxyz_r_0(double fga,
                                     double fgb,
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
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyz,
                                     double pb_xxz,
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_xy * fx * fx * fz * fgb

                     + 2.5 * pa_xy * fx * fx * fx * fz

                     + 5.0 * fx * fx * fx * pa_y * fz * pb_x

                     - 0.25 * pa_x * fx * fx * fz * fgb * pb_y

                     - 0.25 * pa_x * fz * fga * fx * fx * pb_y

                     - 0.5 * fx * fx * fz * fga * pb_xy

                     - 0.5 * pa_xyy * fx * fz * fgb * pb_y

                     - pa_xyz * fx * fz * fgb * pb_z

                     + 1.25 * pa_x * fx * fx * fx * fz * pb_y

                     + 2.5 * fx * fx * fx * fz * pb_xy

                     + 3.0 * pa_xyy * fz * fx * fx * pb_y

                     + 6.0 * pa_xyz * fx * fx * fz * pb_z

                     + 6.0 * pa_xy * fx * fx * fz * pb_xx

                     + 6.0 * fx * fx * pa_yy * fz * pb_xy

                     + 12.0 * fx * fx * pa_yz * fz * pb_xz

                     - 0.5 * pa_xz * fx * fz * fgb * pb_yz

                     - 0.5 * pa_xz * fz * fga * fx * pb_yz

                     - 0.5 * pa_x * fz * fga * fx * pb_xxy

                     - fx * fz * fga * pa_z * pb_xyz

                     - pa_xyyz * fz * fgb * pb_yz

                     + 3.0 * pa_xz * fx * fx * fz * pb_yz

                     + 3.0 * pa_x * fx * fx * fz * pb_xxy

                     + 6.0 * fx * fx * fz * pa_z * pb_xyz

                     + 7.0 * pa_xyyz * fz * fx * pb_yz

                     + 7.0 * pa_xyy * fz * fx * pb_xxy

                     + 14.0 * pa_xyz * fx * fz * pb_xxz

                     + 14.0 * fx * pa_yyz * fz * pb_xyz

                     - pa_xz * fz * fga * pb_xxyz

                     + 7.0 * pa_xz * fx * fz * pb_xxyz

                     + 16.0 * pa_xyyz * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xxzz_s_0(double fx,
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
                                     double pb_xxzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.125 * pa_xz * fx * fx * fx

                     + 0.25 * pa_x * fx * fx * fx * pb_z

                     + 0.25 * fx * fx * fx * pa_z * pb_x

                     + 0.5 * fx * fx * fx * pb_xz

                     + 0.25 * pa_xyyz * fx * fx

                     + 0.5 * pa_xyy * fx * fx * pb_z

                     + 0.5 * fx * fx * pa_yyz * pb_x

                     + fx * fx * pa_yy * pb_xz

                     + 0.25 * pa_xz * fx * fx * pb_xx

                     + 0.25 * pa_xz * fx * fx * pb_zz

                     + 0.5 * pa_x * fx * fx * pb_xxz

                     + 0.5 * fx * fx * pa_z * pb_xzz

                     + 0.5 * pa_xyyz * pb_xx * fx

                     + 0.5 * pa_xyyz * fx * pb_zz

                     + pa_xyy * fx * pb_xxz

                     + fx * pa_yyz * pb_xzz

                     + 0.5 * pa_xz * fx * pb_xxzz

                     + pa_xyyz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xxzz_r_0(double fga,
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
                                     double pb_xxzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_xz * fx * fx * fz * fgb

                     - 0.5 * pa_x * fx * fx * fz * fgb * pb_z

                     - 0.25 * pa_xz * fz * fga * fx * fx

                     - 0.5 * pa_x * fz * fga * fx * fx * pb_z

                     - 0.5 * fx * fx * pa_z * pb_x * fz * fgb

                     - 0.5 * fx * fx * fz * fga * pa_z * pb_x

                     - fx * fx * fz * fga * pb_xz

                     - pa_xyyz * fx * fz * fgb

                     - pa_xyy * fx * fz * fgb * pb_z

                     + 1.25 * pa_xz * fx * fx * fx * fz

                     + 2.5 * pa_x * fx * fx * fx * fz * pb_z

                     - fx * pa_yyz * pb_x * fz * fgb

                     + 2.5 * fx * fx * fx * fz * pa_z * pb_x

                     + 5.0 * fx * fx * fx * fz * pb_xz

                     + 3.0 * pa_xyyz * fz * fx * fx

                     + 6.0 * pa_xyy * fz * fx * fx * pb_z

                     + 6.0 * fx * fx * pa_yyz * fz * pb_x

                     + 12.0 * fx * fx * pa_yy * fz * pb_xz

                     - 0.5 * pa_xz * fx * pb_xx * fz * fgb

                     - 0.5 * pa_xz * fx * fz * fgb * pb_zz

                     - 0.5 * pa_xz * fz * fga * pb_xx * fx

                     - 0.5 * pa_xz * fz * fga * fx * pb_zz

                     - pa_x * fz * fga * fx * pb_xxz

                     - fx * fz * fga * pa_z * pb_xzz

                     - pa_xyyz * pb_xx * fz * fgb

                     - pa_xyyz * fz * fgb * pb_zz

                     + 3.0 * pa_xz * fx * fx * fz * pb_xx

                     + 3.0 * pa_xz * fx * fx * fz * pb_zz

                     + 6.0 * pa_x * fx * fx * fz * pb_xxz

                     + 6.0 * fx * fx * fz * pa_z * pb_xzz

                     + 7.0 * pa_xyyz * fz * pb_xx * fx

                     + 7.0 * pa_xyyz * fz * fx * pb_zz

                     + 14.0 * pa_xyy * fz * fx * pb_xxz

                     + 14.0 * fx * pa_yyz * fz * pb_xzz

                     - pa_xz * fz * fga * pb_xxzz

                     + 7.0 * pa_xz * fx * fz * pb_xxzz

                     + 16.0 * pa_xyyz * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xyyy_s_0(double fx,
                                     double pa_xyyz,
                                     double pa_xyz,
                                     double pa_xz,
                                     double pa_yyz,
                                     double pa_yz,
                                     double pa_z,
                                     double pb_x,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_xyyy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * fx * pa_yz

                     + 1.125 * fx * fx * fx * pa_z * pb_y

                     + 1.5 * pa_xyz * fx * fx * pb_x

                     + 2.25 * pa_xz * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_yyz * pb_y

                     + 1.5 * fx * fx * pa_yz * pb_yy

                     + 0.25 * fx * fx * pa_z * pb_yyy

                     + 1.5 * pa_xyyz * pb_xy * fx

                     + 3.0 * pa_xyz * fx * pb_xyy

                     + 0.5 * fx * pa_yyz * pb_yyy

                     + 0.5 * pa_xz * fx * pb_xyyy

                     + pa_xyyz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xyyy_r_0(double fga,
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
                                     double pb_xyyy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_yz * fz * fgb

                     + 7.5 * fx * fx * fx * pa_yz * fz

                     + 11.25 * fx * fx * fx * fz * pa_z * pb_y

                     - 0.75 * fx * fx * pa_z * pb_y * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pa_z * pb_y

                     - 3.0 * pa_xyz * fx * pb_x * fz * fgb

                     - 1.5 * fx * pa_yyz * pb_y * fz * fgb

                     + 18.0 * pa_xyz * fx * fx * fz * pb_x

                     + 27.0 * pa_xz * fx * fx * fz * pb_xy

                     + 9.0 * fx * fx * pa_yyz * fz * pb_y

                     + 18.0 * fx * fx * pa_yz * fz * pb_yy

                     - 1.5 * pa_xz * fx * pb_xy * fz * fgb

                     - 1.5 * pa_xz * fz * fga * pb_xy * fx

                     - 0.5 * fx * fz * fga * pa_z * pb_yyy

                     - 3.0 * pa_xyyz * pb_xy * fz * fgb

                     + 3.0 * fx * fx * fz * pa_z * pb_yyy

                     + 21.0 * pa_xyyz * fz * pb_xy * fx

                     + 42.0 * pa_xyz * fx * fz * pb_xyy

                     + 7.0 * fx * pa_yyz * fz * pb_yyy

                     - pa_xz * fz * fga * pb_xyyy

                     + 7.0 * pa_xz * fx * fz * pb_xyyy

                     + 16.0 * pa_xyyz * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xyyz_s_0(double fx,
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
                                     double pb_xyy,
                                     double pb_xyyz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * pa_x * fx * fx * fx * pb_x

                     + 0.125 * fx * fx * fx * pa_yy

                     + 0.5 * fx * fx * fx * pa_y * pb_y

                     + 0.375 * fx * fx * fx * pa_z * pb_z

                     + 0.125 * fx * fx * fx * pb_yy

                     + 0.25 * pa_xyy * fx * fx * pb_x

                     + pa_xy * fx * fx * pb_xy

                     + 0.75 * pa_xz * fx * fx * pb_xz

                     + 0.25 * fx * fx * pa_yyz * pb_z

                     + 0.25 * fx * fx * pa_yy * pb_yy

                     + fx * fx * pa_yz * pb_yz

                     + 0.25 * pa_x * fx * fx * pb_xyy

                     + 0.25 * fx * fx * pa_z * pb_yyz

                     + 0.5 * pa_xyyz * pb_xz * fx

                     + 0.5 * pa_xyy * fx * pb_xyy

                     + 2.0 * pa_xyz * fx * pb_xyz

                     + 0.5 * fx * pa_yyz * pb_yyz

                     + 0.5 * pa_xz * fx * pb_xyyz

                     + pa_xyyz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xyyz_r_0(double fga,
                                     double fgb,
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
                                     double pb_xyy,
                                     double pb_xyyz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (1.5 * fx * fx * fx * fx * fz

                     - 0.125 * fx * fx * fx * fz * fgb

                     - 0.125 * fx * fx * fx * fz * fga

                     - 0.25 * fx * fx * pa_yy * fz * fgb

                     + 3.75 * pa_x * fx * fx * fx * fz * pb_x

                     + 1.25 * fx * fx * fx * pa_yy * fz

                     + 5.0 * fx * fx * fx * pa_y * fz * pb_y

                     + 3.75 * fx * fx * fx * fz * pa_z * pb_z

                     - 0.25 * pa_x * fx * fx * pb_x * fz * fgb

                     - 0.25 * pa_x * fz * fga * fx * fx * pb_x

                     - 0.25 * fx * fx * pa_z * fz * fgb * pb_z

                     - 0.25 * fx * fx * fz * fga * pa_z * pb_z

                     - 0.25 * fx * fx * fz * fga * pb_yy

                     - 0.5 * pa_xyy * fx * pb_x * fz * fgb

                     - 0.5 * fx * pa_yyz * fz * fgb * pb_z

                     + 1.25 * fx * fx * fx * fz * pb_yy

                     + 3.0 * pa_xyy * fz * fx * fx * pb_x

                     + 12.0 * pa_xy * fx * fx * fz * pb_xy

                     + 9.0 * pa_xz * fx * fx * fz * pb_xz

                     + 3.0 * fx * fx * pa_yyz * fz * pb_z

                     + 3.0 * fx * fx * pa_yy * fz * pb_yy

                     + 12.0 * fx * fx * pa_yz * fz * pb_yz

                     - 0.5 * pa_xz * fx * pb_xz * fz * fgb

                     - 0.5 * pa_xz * fz * fga * pb_xz * fx

                     - 0.5 * pa_x * fz * fga * fx * pb_xyy

                     - 0.5 * fx * fz * fga * pa_z * pb_yyz

                     - pa_xyyz * pb_xz * fz * fgb

                     + 3.0 * pa_x * fx * fx * fz * pb_xyy

                     + 3.0 * fx * fx * fz * pa_z * pb_yyz

                     + 7.0 * pa_xyyz * fz * pb_xz * fx

                     + 7.0 * pa_xyy * fz * fx * pb_xyy

                     + 28.0 * pa_xyz * fx * fz * pb_xyz

                     + 7.0 * fx * pa_yyz * fz * pb_yyz

                     - pa_xz * fz * fga * pb_xyyz

                     + 7.0 * pa_xz * fx * fz * pb_xyyz

                     + 16.0 * pa_xyyz * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xyzz_s_0(double fx,
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
                                     double pb_xyzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * fx * pa_yz

                     + 0.5 * fx * fx * fx * pa_y * pb_z

                     + 0.125 * fx * fx * fx * pa_z * pb_y

                     + 0.25 * fx * fx * fx * pb_yz

                     + 0.5 * pa_xyz * fx * fx * pb_x

                     + pa_xy * fx * fx * pb_xz

                     + 0.25 * fx * fx * pa_yyz * pb_y

                     + 0.5 * fx * fx * pa_yy * pb_yz

                     + 0.5 * fx * fx * pa_yz * pb_zz

                     + 0.25 * pa_xz * fx * fx * pb_xy

                     + 0.5 * pa_x * fx * fx * pb_xyz

                     + 0.25 * fx * fx * pa_z * pb_yzz

                     + 0.5 * pa_xyyz * pb_xy * fx

                     + pa_xyy * fx * pb_xyz

                     + pa_xyz * fx * pb_xzz

                     + 0.5 * fx * pa_yyz * pb_yzz

                     + 0.5 * pa_xz * fx * pb_xyzz

                     + pa_xyyz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xyzz_r_0(double fga,
                                     double fgb,
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
                                     double pb_xyzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * pa_yz * fz * fgb

                     + 2.5 * fx * fx * fx * pa_yz * fz

                     + 5.0 * fx * fx * fx * pa_y * fz * pb_z

                     - 0.25 * fx * fx * pa_z * pb_y * fz * fgb

                     - 0.25 * fx * fx * fz * fga * pa_z * pb_y

                     - 0.5 * fx * fx * fz * fga * pb_yz

                     - pa_xyz * fx * pb_x * fz * fgb

                     - 0.5 * fx * pa_yyz * pb_y * fz * fgb

                     + 1.25 * fx * fx * fx * fz * pa_z * pb_y

                     + 2.5 * fx * fx * fx * fz * pb_yz

                     + 6.0 * pa_xyz * fx * fx * fz * pb_x

                     + 12.0 * pa_xy * fx * fx * fz * pb_xz

                     + 3.0 * fx * fx * pa_yyz * fz * pb_y

                     + 6.0 * fx * fx * pa_yy * fz * pb_yz

                     + 6.0 * fx * fx * pa_yz * fz * pb_zz

                     - 0.5 * pa_xz * fx * pb_xy * fz * fgb

                     - 0.5 * pa_xz * fz * fga * pb_xy * fx

                     - pa_x * fz * fga * fx * pb_xyz

                     - 0.5 * fx * fz * fga * pa_z * pb_yzz

                     - pa_xyyz * pb_xy * fz * fgb

                     + 3.0 * pa_xz * fx * fx * fz * pb_xy

                     + 6.0 * pa_x * fx * fx * fz * pb_xyz

                     + 3.0 * fx * fx * fz * pa_z * pb_yzz

                     + 7.0 * pa_xyyz * fz * pb_xy * fx

                     + 14.0 * pa_xyy * fz * fx * pb_xyz

                     + 14.0 * pa_xyz * fx * fz * pb_xzz

                     + 7.0 * fx * pa_yyz * fz * pb_yzz

                     - pa_xz * fz * fga * pb_xyzz

                     + 7.0 * pa_xz * fx * fz * pb_xyzz

                     + 16.0 * pa_xyyz * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xzzz_s_0(double fx,
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
                                     double pb_xzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * fx * fx * fx * pa_yy

                     + 0.375 * pa_x * fx * fx * fx * pb_x

                     + 0.375 * fx * fx * fx * pa_z * pb_z

                     + 0.375 * fx * fx * fx * pb_zz

                     + 0.75 * pa_xyy * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_yyz * pb_z

                     + 0.75 * fx * fx * pa_yy * pb_zz

                     + 0.75 * pa_xz * fx * fx * pb_xz

                     + 0.75 * pa_x * fx * fx * pb_xzz

                     + 0.25 * fx * fx * pa_z * pb_zzz

                     + 1.5 * pa_xyyz * pb_xz * fx

                     + 1.5 * pa_xyy * fx * pb_xzz

                     + 0.5 * fx * pa_yyz * pb_zzz

                     + 0.5 * pa_xz * fx * pb_xzzz

                     + pa_xyyz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xzzz_r_0(double fga,
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
                                     double pb_xzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.375 * fx * fx * fx * fz * fgb

                     - 0.375 * fx * fx * fx * fz * fga

                     - 0.75 * fx * fx * pa_yy * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     + 3.75 * fx * fx * fx * pa_yy * fz

                     - 0.75 * pa_x * fx * fx * pb_x * fz * fgb

                     - 0.75 * pa_x * fz * fga * fx * fx * pb_x

                     - 0.75 * fx * fx * pa_z * pb_z * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pa_z * pb_z

                     - 0.75 * fx * fx * fz * fga * pb_zz

                     - 1.5 * pa_xyy * fx * pb_x * fz * fgb

                     + 3.75 * pa_x * fx * fx * fx * fz * pb_x

                     - 1.5 * fx * pa_yyz * pb_z * fz * fgb

                     + 3.75 * fx * fx * fx * fz * pa_z * pb_z

                     + 3.75 * fx * fx * fx * fz * pb_zz

                     + 9.0 * pa_xyy * fz * fx * fx * pb_x

                     + 9.0 * fx * fx * pa_yyz * fz * pb_z

                     + 9.0 * fx * fx * pa_yy * fz * pb_zz

                     - 1.5 * pa_xz * fx * pb_xz * fz * fgb

                     - 1.5 * pa_xz * fz * fga * pb_xz * fx

                     - 1.5 * pa_x * fz * fga * fx * pb_xzz

                     - 0.5 * fx * fz * fga * pa_z * pb_zzz

                     - 3.0 * pa_xyyz * pb_xz * fz * fgb

                     + 9.0 * pa_xz * fx * fx * fz * pb_xz

                     + 9.0 * pa_x * fx * fx * fz * pb_xzz

                     + 3.0 * fx * fx * fz * pa_z * pb_zzz

                     + 21.0 * pa_xyyz * fz * pb_xz * fx

                     + 21.0 * pa_xyy * fz * fx * pb_xzz

                     + 7.0 * fx * pa_yyz * fz * pb_zzz

                     - pa_xz * fz * fga * pb_xzzz

                     + 7.0 * pa_xz * fx * fz * pb_xzzz

                     + 16.0 * pa_xyyz * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_yyyy_s_0(double fx,
                                     double pa_xyyz,
                                     double pa_xyz,
                                     double pa_xz,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double pb_yyyy,
                                     double s_0_0)
    {
        return s_0_0 * (1.875 * pa_xz * fx * fx * fx

                     + 0.75 * pa_xyyz * fx * fx

                     + 6.0 * pa_xyz * fx * fx * pb_y

                     + 4.5 * pa_xz * fx * fx * pb_yy

                     + 3.0 * pa_xyyz * pb_yy * fx

                     + 4.0 * pa_xyz * fx * pb_yyy

                     + 0.5 * pa_xz * fx * pb_yyyy

                     + pa_xyyz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_yyyy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xyyz,
                                     double pa_xyz,
                                     double pa_xz,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double pb_yyyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_xz * fx * fx * fz * fgb

                     + 18.75 * pa_xz * fx * fx * fx * fz

                     - 0.75 * pa_xz * fz * fga * fx * fx

                     - 3.0 * pa_xyyz * fx * fz * fgb

                     - 12.0 * pa_xyz * fx * pb_y * fz * fgb

                     + 9.0 * pa_xyyz * fz * fx * fx

                     + 72.0 * pa_xyz * fx * fx * fz * pb_y

                     + 54.0 * pa_xz * fx * fx * fz * pb_yy

                     - 3.0 * pa_xz * fx * pb_yy * fz * fgb

                     - 3.0 * pa_xz * fz * fga * pb_yy * fx

                     - 6.0 * pa_xyyz * pb_yy * fz * fgb

                     + 42.0 * pa_xyyz * fz * pb_yy * fx

                     + 56.0 * pa_xyz * fx * fz * pb_yyy

                     - pa_xz * fz * fga * pb_yyyy

                     + 7.0 * pa_xz * fx * fz * pb_yyyy

                     + 16.0 * pa_xyyz * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_yyyz_s_0(double fx,
                                     double pa_x,
                                     double pa_xy,
                                     double pa_xyy,
                                     double pa_xyyz,
                                     double pa_xyz,
                                     double pa_xz,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double pb_yyyz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xy * fx * fx * fx

                     + 1.125 * pa_x * fx * fx * fx * pb_y

                     + 0.75 * pa_xyy * fx * fx * pb_y

                     + 1.5 * pa_xyz * fx * fx * pb_z

                     + 1.5 * pa_xy * fx * fx * pb_yy

                     + 2.25 * pa_xz * fx * fx * pb_yz

                     + 0.25 * pa_x * fx * fx * pb_yyy

                     + 1.5 * pa_xyyz * pb_yz * fx

                     + 0.5 * pa_xyy * fx * pb_yyy

                     + 3.0 * pa_xyz * fx * pb_yyz

                     + 0.5 * pa_xz * fx * pb_yyyz

                     + pa_xyyz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_yyyz_r_0(double fga,
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
                                     double pb_yyy,
                                     double pb_yyyz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xy * fx * fx * fz * fgb

                     + 7.5 * pa_xy * fx * fx * fx * fz

                     + 11.25 * pa_x * fx * fx * fx * fz * pb_y

                     - 0.75 * pa_x * fx * fx * pb_y * fz * fgb

                     - 0.75 * pa_x * fz * fga * fx * fx * pb_y

                     - 1.5 * pa_xyy * fx * pb_y * fz * fgb

                     - 3.0 * pa_xyz * fx * fz * fgb * pb_z

                     + 9.0 * pa_xyy * fz * fx * fx * pb_y

                     + 18.0 * pa_xyz * fx * fx * fz * pb_z

                     + 18.0 * pa_xy * fx * fx * fz * pb_yy

                     + 27.0 * pa_xz * fx * fx * fz * pb_yz

                     - 1.5 * pa_xz * fx * pb_yz * fz * fgb

                     - 1.5 * pa_xz * fz * fga * pb_yz * fx

                     - 0.5 * pa_x * fz * fga * fx * pb_yyy

                     - 3.0 * pa_xyyz * pb_yz * fz * fgb

                     + 3.0 * pa_x * fx * fx * fz * pb_yyy

                     + 21.0 * pa_xyyz * fz * pb_yz * fx

                     + 7.0 * pa_xyy * fz * fx * pb_yyy

                     + 42.0 * pa_xyz * fx * fz * pb_yyz

                     - pa_xz * fz * fga * pb_yyyz

                     + 7.0 * pa_xz * fx * fz * pb_yyyz

                     + 16.0 * pa_xyyz * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_yyzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xy,
                                     double pa_xyy,
                                     double pa_xyyz,
                                     double pa_xyz,
                                     double pa_xz,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_yyzz,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * pa_xz * fx * fx * fx

                     + 0.75 * pa_x * fx * fx * fx * pb_z

                     + 0.25 * pa_xyyz * fx * fx

                     + 0.5 * pa_xyy * fx * fx * pb_z

                     + pa_xyz * fx * fx * pb_y

                     + 2.0 * pa_xy * fx * fx * pb_yz

                     + 0.75 * pa_xz * fx * fx * pb_zz

                     + 0.25 * pa_xz * fx * fx * pb_yy

                     + 0.5 * pa_x * fx * fx * pb_yyz

                     + 0.5 * pa_xyyz * pb_yy * fx

                     + 0.5 * pa_xyyz * fx * pb_zz

                     + pa_xyy * fx * pb_yyz

                     + 2.0 * pa_xyz * fx * pb_yzz

                     + 0.5 * pa_xz * fx * pb_yyzz

                     + pa_xyyz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_yyzz_r_0(double fga,
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
                                     double pb_yyzz,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- pa_xz * fx * fx * fz * fgb

                     + 3.75 * pa_xz * fx * fx * fx * fz

                     + 7.5 * pa_x * fx * fx * fx * fz * pb_z

                     - 0.5 * pa_x * fx * fx * fz * fgb * pb_z

                     - 0.25 * pa_xz * fz * fga * fx * fx

                     - 0.5 * pa_x * fz * fga * fx * fx * pb_z

                     - pa_xyyz * fx * fz * fgb

                     - pa_xyy * fx * fz * fgb * pb_z

                     - 2.0 * pa_xyz * fx * pb_y * fz * fgb

                     + 3.0 * pa_xyyz * fz * fx * fx

                     + 6.0 * pa_xyy * fz * fx * fx * pb_z

                     + 12.0 * pa_xyz * fx * fx * fz * pb_y

                     + 24.0 * pa_xy * fx * fx * fz * pb_yz

                     + 9.0 * pa_xz * fx * fx * fz * pb_zz

                     - 0.5 * pa_xz * fx * pb_yy * fz * fgb

                     - 0.5 * pa_xz * fx * fz * fgb * pb_zz

                     - 0.5 * pa_xz * fz * fga * pb_yy * fx

                     - 0.5 * pa_xz * fz * fga * fx * pb_zz

                     - pa_x * fz * fga * fx * pb_yyz

                     - pa_xyyz * pb_yy * fz * fgb

                     - pa_xyyz * fz * fgb * pb_zz

                     + 3.0 * pa_xz * fx * fx * fz * pb_yy

                     + 6.0 * pa_x * fx * fx * fz * pb_yyz

                     + 7.0 * pa_xyyz * fz * pb_yy * fx

                     + 7.0 * pa_xyyz * fz * fx * pb_zz

                     + 14.0 * pa_xyy * fz * fx * pb_yyz

                     + 28.0 * pa_xyz * fx * fz * pb_yzz

                     - pa_xz * fz * fga * pb_yyzz

                     + 7.0 * pa_xz * fx * fz * pb_yyzz

                     + 16.0 * pa_xyyz * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_yzzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xy,
                                     double pa_xyy,
                                     double pa_xyyz,
                                     double pa_xyz,
                                     double pa_xz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_yzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xy * fx * fx * fx

                     + 0.375 * pa_x * fx * fx * fx * pb_y

                     + 0.75 * pa_xyy * fx * fx * pb_y

                     + 1.5 * pa_xyz * fx * fx * pb_z

                     + 1.5 * pa_xy * fx * fx * pb_zz

                     + 0.75 * pa_xz * fx * fx * pb_yz

                     + 0.75 * pa_x * fx * fx * pb_yzz

                     + 1.5 * pa_xyyz * pb_yz * fx

                     + 1.5 * pa_xyy * fx * pb_yzz

                     + pa_xyz * fx * pb_zzz

                     + 0.5 * pa_xz * fx * pb_yzzz

                     + pa_xyyz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_yzzz_r_0(double fga,
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
                                     double pb_yzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xy * fx * fx * fz * fgb

                     + 7.5 * pa_xy * fx * fx * fx * fz

                     - 0.75 * pa_x * fx * fx * pb_y * fz * fgb

                     - 0.75 * pa_x * fz * fga * fx * fx * pb_y

                     - 1.5 * pa_xyy * fx * pb_y * fz * fgb

                     - 3.0 * pa_xyz * fx * pb_z * fz * fgb

                     + 3.75 * pa_x * fx * fx * fx * fz * pb_y

                     + 9.0 * pa_xyy * fz * fx * fx * pb_y

                     + 18.0 * pa_xyz * fx * fx * fz * pb_z

                     + 18.0 * pa_xy * fx * fx * fz * pb_zz

                     - 1.5 * pa_xz * fx * pb_yz * fz * fgb

                     - 1.5 * pa_xz * fz * fga * pb_yz * fx

                     - 1.5 * pa_x * fz * fga * fx * pb_yzz

                     - 3.0 * pa_xyyz * pb_yz * fz * fgb

                     + 9.0 * pa_xz * fx * fx * fz * pb_yz

                     + 9.0 * pa_x * fx * fx * fz * pb_yzz

                     + 21.0 * pa_xyyz * fz * pb_yz * fx

                     + 21.0 * pa_xyy * fz * fx * pb_yzz

                     + 14.0 * pa_xyz * fx * fz * pb_zzz

                     - pa_xz * fz * fga * pb_yzzz

                     + 7.0 * pa_xz * fx * fz * pb_yzzz

                     + 16.0 * pa_xyyz * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_zzzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xyy,
                                     double pa_xyyz,
                                     double pa_xz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double pb_zzzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * pa_xz * fx * fx * fx

                     + 1.5 * pa_x * fx * fx * fx * pb_z

                     + 0.75 * pa_xyyz * fx * fx

                     + 3.0 * pa_xyy * fx * fx * pb_z

                     + 1.5 * pa_xz * fx * fx * pb_zz

                     + pa_x * fx * fx * pb_zzz

                     + 3.0 * pa_xyyz * pb_zz * fx

                     + 2.0 * pa_xyy * fx * pb_zzz

                     + 0.5 * pa_xz * fx * pb_zzzz

                     + pa_xyyz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_zzzz_r_0(double fga,
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
                                     double pb_zzzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xz * fx * fx * fz * fgb

                     - 3.0 * pa_x * fx * fx * pb_z * fz * fgb

                     - 0.75 * pa_xz * fz * fga * fx * fx

                     - 3.0 * pa_x * fz * fga * fx * fx * pb_z

                     - 3.0 * pa_xyyz * fx * fz * fgb

                     - 6.0 * pa_xyy * fx * pb_z * fz * fgb

                     + 3.75 * pa_xz * fx * fx * fx * fz

                     + 15.0 * pa_x * fx * fx * fx * fz * pb_z

                     + 9.0 * pa_xyyz * fz * fx * fx

                     + 36.0 * pa_xyy * fz * fx * fx * pb_z

                     - 3.0 * pa_xz * fx * pb_zz * fz * fgb

                     - 3.0 * pa_xz * fz * fga * pb_zz * fx

                     - 2.0 * pa_x * fz * fga * fx * pb_zzz

                     - 6.0 * pa_xyyz * pb_zz * fz * fgb

                     + 18.0 * pa_xz * fx * fx * fz * pb_zz

                     + 12.0 * pa_x * fx * fx * fz * pb_zzz

                     + 42.0 * pa_xyyz * fz * pb_zz * fx

                     + 28.0 * pa_xyy * fz * fx * pb_zzz

                     - pa_xz * fz * fga * pb_zzzz

                     + 7.0 * pa_xz * fx * fz * pb_zzzz

                     + 16.0 * pa_xyyz * fz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xxxx_s_0(double fx,
                                     double pa_xy,
                                     double pa_xyzz,
                                     double pa_y,
                                     double pa_yzz,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxx,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * pa_xy * fx * fx * fx

                     + 1.5 * fx * fx * fx * pa_y * pb_x

                     + 0.75 * pa_xyzz * fx * fx

                     + 3.0 * fx * fx * pa_yzz * pb_x

                     + 1.5 * pa_xy * fx * fx * pb_xx

                     + fx * fx * pa_y * pb_xxx

                     + 3.0 * pa_xyzz * pb_xx * fx

                     + 2.0 * fx * pa_yzz * pb_xxx

                     + 0.5 * pa_xy * fx * pb_xxxx

                     + pa_xyzz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xxxx_r_0(double fga,
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
                                     double pb_xxxx,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xy * fx * fx * fz * fgb

                     - 0.75 * pa_xy * fz * fga * fx * fx

                     - 3.0 * fx * fx * pa_y * pb_x * fz * fgb

                     - 3.0 * fx * fx * pa_y * fz * fga * pb_x

                     - 3.0 * pa_xyzz * fx * fz * fgb

                     + 3.75 * pa_xy * fz * fx * fx * fx

                     - 6.0 * fx * pa_yzz * pb_x * fz * fgb

                     + 15.0 * fx * fx * fx * pa_y * fz * pb_x

                     + 9.0 * pa_xyzz * fz * fx * fx

                     + 36.0 * fx * fx * pa_yzz * fz * pb_x

                     - 3.0 * pa_xy * fx * pb_xx * fz * fgb

                     - 3.0 * pa_xy * fz * fga * pb_xx * fx

                     - 2.0 * fx * pa_y * fz * fga * pb_xxx

                     - 6.0 * pa_xyzz * pb_xx * fz * fgb

                     + 18.0 * pa_xy * fz * fx * fx * pb_xx

                     + 12.0 * fx * fx * pa_y * fz * pb_xxx

                     + 42.0 * pa_xyzz * fz * pb_xx * fx

                     + 28.0 * fx * pa_yzz * fz * pb_xxx

                     - pa_xy * fz * fga * pb_xxxx

                     + 7.0 * pa_xy * fz * fx * pb_xxxx

                     + 16.0 * pa_xyzz * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xxxy_s_0(double fx,
                                     double pa_x,
                                     double pa_xy,
                                     double pa_xyzz,
                                     double pa_xzz,
                                     double pa_y,
                                     double pa_yzz,
                                     double pa_zz,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxy,
                                     double pb_xxy,
                                     double pb_xy,
                                     double pb_y,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * fx * fx * fx * pa_zz

                     + 0.375 * pa_x * fx * fx * fx * pb_x

                     + 0.375 * fx * fx * fx * pa_y * pb_y

                     + 0.375 * fx * fx * fx * pb_xx

                     + 0.75 * pa_xzz * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_yzz * pb_y

                     + 0.75 * fx * fx * pa_zz * pb_xx

                     + 0.75 * pa_xy * fx * fx * pb_xy

                     + 0.25 * pa_x * fx * fx * pb_xxx

                     + 0.75 * fx * fx * pa_y * pb_xxy

                     + 1.5 * pa_xyzz * pb_xy * fx

                     + 0.5 * pa_xzz * fx * pb_xxx

                     + 1.5 * fx * pa_yzz * pb_xxy

                     + 0.5 * pa_xy * fx * pb_xxxy

                     + pa_xyzz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xxxy_r_0(double fga,
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
                                     double pb_xxx,
                                     double pb_xxxy,
                                     double pb_xxy,
                                     double pb_xy,
                                     double pb_y,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.375 * fx * fx * fx * fz * fgb

                     - 0.375 * fx * fx * fx * fz * fga

                     - 0.75 * fx * fx * pa_zz * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     + 3.75 * fx * fx * fx * pa_zz * fz

                     - 0.75 * pa_x * fx * fx * pb_x * fz * fgb

                     - 0.75 * pa_x * fx * fx * fz * fga * pb_x

                     - 0.75 * fx * fx * pa_y * fz * fgb * pb_y

                     - 0.75 * fx * fx * pa_y * fz * fga * pb_y

                     - 0.75 * fx * fx * fz * fga * pb_xx

                     - 1.5 * pa_xzz * fx * pb_x * fz * fgb

                     + 3.75 * pa_x * fx * fx * fx * fz * pb_x

                     - 1.5 * fx * pa_yzz * fz * fgb * pb_y

                     + 3.75 * fx * fx * fx * pa_y * fz * pb_y

                     + 3.75 * fx * fx * fx * fz * pb_xx

                     + 9.0 * pa_xzz * fx * fx * fz * pb_x

                     + 9.0 * fx * fx * pa_yzz * fz * pb_y

                     + 9.0 * fx * fx * pa_zz * fz * pb_xx

                     - 1.5 * pa_xy * fx * pb_xy * fz * fgb

                     - 1.5 * pa_xy * fz * fga * pb_xy * fx

                     - 0.5 * pa_x * fx * fz * fga * pb_xxx

                     - 1.5 * fx * pa_y * fz * fga * pb_xxy

                     - 3.0 * pa_xyzz * pb_xy * fz * fgb

                     + 9.0 * pa_xy * fz * fx * fx * pb_xy

                     + 3.0 * pa_x * fx * fx * fz * pb_xxx

                     + 9.0 * fx * fx * pa_y * fz * pb_xxy

                     + 21.0 * pa_xyzz * fz * pb_xy * fx

                     + 7.0 * pa_xzz * fx * fz * pb_xxx

                     + 21.0 * fx * pa_yzz * fz * pb_xxy

                     - pa_xy * fz * fga * pb_xxxy

                     + 7.0 * pa_xy * fz * fx * pb_xxxy

                     + 16.0 * pa_xyzz * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xxxz_s_0(double fx,
                                     double pa_xy,
                                     double pa_xyz,
                                     double pa_xyzz,
                                     double pa_y,
                                     double pa_yz,
                                     double pa_yzz,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxz,
                                     double pb_xxz,
                                     double pb_xz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * fx * pa_yz

                     + 0.375 * fx * fx * fx * pa_y * pb_z

                     + 1.5 * pa_xyz * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_yzz * pb_z

                     + 1.5 * fx * fx * pa_yz * pb_xx

                     + 0.75 * pa_xy * fx * fx * pb_xz

                     + 0.75 * fx * fx * pa_y * pb_xxz

                     + 1.5 * pa_xyzz * pb_xz * fx

                     + pa_xyz * fx * pb_xxx

                     + 1.5 * fx * pa_yzz * pb_xxz

                     + 0.5 * pa_xy * fx * pb_xxxz

                     + pa_xyzz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xxxz_r_0(double fga,
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
                                     double pb_xxx,
                                     double pb_xxxz,
                                     double pb_xxz,
                                     double pb_xz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_yz * fz * fgb

                     + 7.5 * fx * fx * fx * pa_yz * fz

                     - 0.75 * fx * fx * pa_y * fz * fgb * pb_z

                     - 0.75 * fx * fx * pa_y * fz * fga * pb_z

                     - 3.0 * pa_xyz * fx * pb_x * fz * fgb

                     - 1.5 * fx * pa_yzz * fz * fgb * pb_z

                     + 3.75 * fx * fx * fx * pa_y * fz * pb_z

                     + 18.0 * pa_xyz * fz * fx * fx * pb_x

                     + 9.0 * fx * fx * pa_yzz * fz * pb_z

                     + 18.0 * fx * fx * pa_yz * fz * pb_xx

                     - 1.5 * pa_xy * fx * pb_xz * fz * fgb

                     - 1.5 * pa_xy * fz * fga * pb_xz * fx

                     - 1.5 * fx * pa_y * fz * fga * pb_xxz

                     - 3.0 * pa_xyzz * pb_xz * fz * fgb

                     + 9.0 * pa_xy * fz * fx * fx * pb_xz

                     + 9.0 * fx * fx * pa_y * fz * pb_xxz

                     + 21.0 * pa_xyzz * fz * pb_xz * fx

                     + 14.0 * pa_xyz * fz * fx * pb_xxx

                     + 21.0 * fx * pa_yzz * fz * pb_xxz

                     - pa_xy * fz * fga * pb_xxxz

                     + 7.0 * pa_xy * fz * fx * pb_xxxz

                     + 16.0 * pa_xyzz * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xxyy_s_0(double fx,
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
                                     double pb_xxyy,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_y,
                                     double pb_yy,
                                     double s_0_0)
    {
        return s_0_0 * (0.125 * pa_xy * fx * fx * fx

                     + 0.25 * pa_x * fx * fx * fx * pb_y

                     + 0.25 * fx * fx * fx * pa_y * pb_x

                     + 0.5 * fx * fx * fx * pb_xy

                     + 0.25 * pa_xyzz * fx * fx

                     + 0.5 * pa_xzz * fx * fx * pb_y

                     + 0.5 * fx * fx * pa_yzz * pb_x

                     + fx * fx * pa_zz * pb_xy

                     + 0.25 * pa_xy * fx * fx * pb_xx

                     + 0.25 * pa_xy * fx * fx * pb_yy

                     + 0.5 * pa_x * fx * fx * pb_xxy

                     + 0.5 * fx * fx * pa_y * pb_xyy

                     + 0.5 * pa_xyzz * pb_xx * fx

                     + 0.5 * pa_xyzz * fx * pb_yy

                     + pa_xzz * fx * pb_xxy

                     + fx * pa_yzz * pb_xyy

                     + 0.5 * pa_xy * fx * pb_xxyy

                     + pa_xyzz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xxyy_r_0(double fga,
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
                                     double pb_xxyy,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_y,
                                     double pb_yy,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_xy * fx * fx * fz * fgb

                     - 0.25 * pa_xy * fz * fga * fx * fx

                     - 0.5 * pa_x * fx * fx * fz * fgb * pb_y

                     - 0.5 * pa_x * fx * fx * fz * fga * pb_y

                     - 0.5 * fx * fx * pa_y * pb_x * fz * fgb

                     - 0.5 * fx * fx * pa_y * fz * fga * pb_x

                     - fx * fx * fz * fga * pb_xy

                     - pa_xyzz * fx * fz * fgb

                     + 1.25 * pa_xy * fz * fx * fx * fx

                     - pa_xzz * fx * fz * fgb * pb_y

                     + 2.5 * pa_x * fx * fx * fx * fz * pb_y

                     - fx * pa_yzz * pb_x * fz * fgb

                     + 2.5 * fx * fx * fx * pa_y * fz * pb_x

                     + 5.0 * fx * fx * fx * fz * pb_xy

                     + 3.0 * pa_xyzz * fz * fx * fx

                     + 6.0 * pa_xzz * fx * fx * fz * pb_y

                     + 6.0 * fx * fx * pa_yzz * fz * pb_x

                     + 12.0 * fx * fx * pa_zz * fz * pb_xy

                     - 0.5 * pa_xy * fx * pb_xx * fz * fgb

                     - 0.5 * pa_xy * fx * fz * fgb * pb_yy

                     - 0.5 * pa_xy * fz * fga * pb_xx * fx

                     - 0.5 * pa_xy * fz * fga * fx * pb_yy

                     - pa_x * fx * fz * fga * pb_xxy

                     - fx * pa_y * fz * fga * pb_xyy

                     - pa_xyzz * pb_xx * fz * fgb

                     - pa_xyzz * fz * fgb * pb_yy

                     + 3.0 * pa_xy * fz * fx * fx * pb_xx

                     + 3.0 * pa_xy * fz * fx * fx * pb_yy

                     + 6.0 * pa_x * fx * fx * fz * pb_xxy

                     + 6.0 * fx * fx * pa_y * fz * pb_xyy

                     + 7.0 * pa_xyzz * fz * pb_xx * fx

                     + 7.0 * pa_xyzz * fz * fx * pb_yy

                     + 14.0 * pa_xzz * fx * fz * pb_xxy

                     + 14.0 * fx * pa_yzz * fz * pb_xyy

                     - pa_xy * fz * fga * pb_xxyy

                     + 7.0 * pa_xy * fz * fx * pb_xxyy

                     + 16.0 * pa_xyzz * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xxyz_s_0(double fx,
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
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyz,
                                     double pb_xxz,
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.25 * pa_xz * fx * fx * fx

                     + 0.5 * fx * fx * fx * pa_z * pb_x

                     + 0.125 * pa_x * fx * fx * fx * pb_z

                     + 0.25 * fx * fx * fx * pb_xz

                     + 0.5 * pa_xyz * fx * fx * pb_y

                     + 0.25 * pa_xzz * fx * fx * pb_z

                     + 0.5 * pa_xz * fx * fx * pb_xx

                     + fx * fx * pa_yz * pb_xy

                     + 0.5 * fx * fx * pa_zz * pb_xz

                     + 0.25 * pa_xy * fx * fx * pb_yz

                     + 0.25 * pa_x * fx * fx * pb_xxz

                     + 0.5 * fx * fx * pa_y * pb_xyz

                     + 0.5 * pa_xyzz * fx * pb_yz

                     + pa_xyz * fx * pb_xxy

                     + 0.5 * pa_xzz * fx * pb_xxz

                     + fx * pa_yzz * pb_xyz

                     + 0.5 * pa_xy * fx * pb_xxyz

                     + pa_xyzz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xxyz_r_0(double fga,
                                     double fgb,
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
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyz,
                                     double pb_xxz,
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_xz * fx * fx * fz * fgb

                     + 2.5 * pa_xz * fx * fx * fx * fz

                     + 5.0 * fx * fx * fx * pa_z * fz * pb_x

                     - 0.25 * pa_x * fx * fx * fz * fgb * pb_z

                     - 0.25 * pa_x * fx * fx * fz * fga * pb_z

                     - 0.5 * fx * fx * fz * fga * pb_xz

                     - pa_xyz * fx * fz * fgb * pb_y

                     - 0.5 * pa_xzz * fx * fz * fgb * pb_z

                     + 1.25 * pa_x * fx * fx * fx * fz * pb_z

                     + 2.5 * fx * fx * fx * fz * pb_xz

                     + 6.0 * pa_xyz * fz * fx * fx * pb_y

                     + 3.0 * pa_xzz * fx * fx * fz * pb_z

                     + 6.0 * pa_xz * fx * fx * fz * pb_xx

                     + 12.0 * fx * fx * pa_yz * fz * pb_xy

                     + 6.0 * fx * fx * pa_zz * fz * pb_xz

                     - 0.5 * pa_xy * fx * fz * fgb * pb_yz

                     - 0.5 * pa_xy * fz * fga * fx * pb_yz

                     - 0.5 * pa_x * fx * fz * fga * pb_xxz

                     - fx * pa_y * fz * fga * pb_xyz

                     - pa_xyzz * fz * fgb * pb_yz

                     + 3.0 * pa_xy * fz * fx * fx * pb_yz

                     + 3.0 * pa_x * fx * fx * fz * pb_xxz

                     + 6.0 * fx * fx * pa_y * fz * pb_xyz

                     + 7.0 * pa_xyzz * fz * fx * pb_yz

                     + 14.0 * pa_xyz * fz * fx * pb_xxy

                     + 7.0 * pa_xzz * fx * fz * pb_xxz

                     + 14.0 * fx * pa_yzz * fz * pb_xyz

                     - pa_xy * fz * fga * pb_xxyz

                     + 7.0 * pa_xy * fz * fx * pb_xxyz

                     + 16.0 * pa_xyzz * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xxzz_s_0(double fx,
                                     double pa_xy,
                                     double pa_xyz,
                                     double pa_xyzz,
                                     double pa_y,
                                     double pa_yz,
                                     double pa_yzz,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxz,
                                     double pb_xxzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * pa_xy * fx * fx * fx

                     + 0.75 * fx * fx * fx * pa_y * pb_x

                     + 0.25 * pa_xyzz * fx * fx

                     + pa_xyz * fx * fx * pb_z

                     + 0.75 * pa_xy * fx * fx * pb_xx

                     + 0.5 * fx * fx * pa_yzz * pb_x

                     + 2.0 * fx * fx * pa_yz * pb_xz

                     + 0.25 * pa_xy * fx * fx * pb_zz

                     + 0.5 * fx * fx * pa_y * pb_xzz

                     + 0.5 * pa_xyzz * pb_xx * fx

                     + 0.5 * pa_xyzz * fx * pb_zz

                     + 2.0 * pa_xyz * fx * pb_xxz

                     + fx * pa_yzz * pb_xzz

                     + 0.5 * pa_xy * fx * pb_xxzz

                     + pa_xyzz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xxzz_r_0(double fga,
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
                                     double pb_xxzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- pa_xy * fx * fx * fz * fgb

                     + 3.75 * pa_xy * fx * fx * fx * fz

                     + 7.5 * fx * fx * fx * pa_y * fz * pb_x

                     - 0.25 * pa_xy * fz * fga * fx * fx

                     - 0.5 * fx * fx * pa_y * pb_x * fz * fgb

                     - 0.5 * fx * fx * pa_y * fz * fga * pb_x

                     - pa_xyzz * fx * fz * fgb

                     - 2.0 * pa_xyz * fx * fz * fgb * pb_z

                     - fx * pa_yzz * pb_x * fz * fgb

                     + 3.0 * pa_xyzz * fz * fx * fx

                     + 12.0 * pa_xyz * fz * fx * fx * pb_z

                     + 9.0 * pa_xy * fx * fx * fz * pb_xx

                     + 6.0 * fx * fx * pa_yzz * fz * pb_x

                     + 24.0 * fx * fx * pa_yz * fz * pb_xz

                     - 0.5 * pa_xy * fx * pb_xx * fz * fgb

                     - 0.5 * pa_xy * fx * fz * fgb * pb_zz

                     - 0.5 * pa_xy * fz * fga * pb_xx * fx

                     - 0.5 * pa_xy * fz * fga * fx * pb_zz

                     - fx * pa_y * fz * fga * pb_xzz

                     - pa_xyzz * pb_xx * fz * fgb

                     - pa_xyzz * fz * fgb * pb_zz

                     + 3.0 * pa_xy * fz * fx * fx * pb_zz

                     + 6.0 * fx * fx * pa_y * fz * pb_xzz

                     + 7.0 * pa_xyzz * fz * pb_xx * fx

                     + 7.0 * pa_xyzz * fz * fx * pb_zz

                     + 28.0 * pa_xyz * fz * fx * pb_xxz

                     + 14.0 * fx * pa_yzz * fz * pb_xzz

                     - pa_xy * fz * fga * pb_xxzz

                     + 7.0 * pa_xy * fz * fx * pb_xxzz

                     + 16.0 * pa_xyzz * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xyyy_s_0(double fx,
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
                                     double pb_xyyy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * fx * fx * fx * pa_zz

                     + 0.375 * pa_x * fx * fx * fx * pb_x

                     + 0.375 * fx * fx * fx * pa_y * pb_y

                     + 0.375 * fx * fx * fx * pb_yy

                     + 0.75 * pa_xzz * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_yzz * pb_y

                     + 0.75 * fx * fx * pa_zz * pb_yy

                     + 0.75 * pa_xy * fx * fx * pb_xy

                     + 0.75 * pa_x * fx * fx * pb_xyy

                     + 0.25 * fx * fx * pa_y * pb_yyy

                     + 1.5 * pa_xyzz * pb_xy * fx

                     + 1.5 * pa_xzz * fx * pb_xyy

                     + 0.5 * fx * pa_yzz * pb_yyy

                     + 0.5 * pa_xy * fx * pb_xyyy

                     + pa_xyzz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xyyy_r_0(double fga,
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
                                     double pb_xyyy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.375 * fx * fx * fx * fz * fgb

                     - 0.375 * fx * fx * fx * fz * fga

                     - 0.75 * fx * fx * pa_zz * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     + 3.75 * fx * fx * fx * pa_zz * fz

                     - 0.75 * pa_x * fx * fx * pb_x * fz * fgb

                     - 0.75 * pa_x * fx * fx * fz * fga * pb_x

                     - 0.75 * fx * fx * pa_y * pb_y * fz * fgb

                     - 0.75 * fx * fx * pa_y * fz * fga * pb_y

                     - 0.75 * fx * fx * fz * fga * pb_yy

                     - 1.5 * pa_xzz * fx * pb_x * fz * fgb

                     + 3.75 * pa_x * fx * fx * fx * fz * pb_x

                     - 1.5 * fx * pa_yzz * pb_y * fz * fgb

                     + 3.75 * fx * fx * fx * pa_y * fz * pb_y

                     + 3.75 * fx * fx * fx * fz * pb_yy

                     + 9.0 * pa_xzz * fx * fx * fz * pb_x

                     + 9.0 * fx * fx * pa_yzz * fz * pb_y

                     + 9.0 * fx * fx * pa_zz * fz * pb_yy

                     - 1.5 * pa_xy * fx * pb_xy * fz * fgb

                     - 1.5 * pa_xy * fz * fga * pb_xy * fx

                     - 1.5 * pa_x * fx * fz * fga * pb_xyy

                     - 0.5 * fx * pa_y * fz * fga * pb_yyy

                     - 3.0 * pa_xyzz * pb_xy * fz * fgb

                     + 9.0 * pa_xy * fz * fx * fx * pb_xy

                     + 9.0 * pa_x * fx * fx * fz * pb_xyy

                     + 3.0 * fx * fx * pa_y * fz * pb_yyy

                     + 21.0 * pa_xyzz * fz * pb_xy * fx

                     + 21.0 * pa_xzz * fx * fz * pb_xyy

                     + 7.0 * fx * pa_yzz * fz * pb_yyy

                     - pa_xy * fz * fga * pb_xyyy

                     + 7.0 * pa_xy * fz * fx * pb_xyyy

                     + 16.0 * pa_xyzz * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xyyz_s_0(double fx,
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
                                     double pb_xyy,
                                     double pb_xyyz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * fx * pa_yz

                     + 0.5 * fx * fx * fx * pa_z * pb_y

                     + 0.125 * fx * fx * fx * pa_y * pb_z

                     + 0.25 * fx * fx * fx * pb_yz

                     + 0.5 * pa_xyz * fx * fx * pb_x

                     + pa_xz * fx * fx * pb_xy

                     + 0.25 * fx * fx * pa_yzz * pb_z

                     + 0.5 * fx * fx * pa_yz * pb_yy

                     + 0.5 * fx * fx * pa_zz * pb_yz

                     + 0.25 * pa_xy * fx * fx * pb_xz

                     + 0.5 * pa_x * fx * fx * pb_xyz

                     + 0.25 * fx * fx * pa_y * pb_yyz

                     + 0.5 * pa_xyzz * pb_xz * fx

                     + pa_xyz * fx * pb_xyy

                     + pa_xzz * fx * pb_xyz

                     + 0.5 * fx * pa_yzz * pb_yyz

                     + 0.5 * pa_xy * fx * pb_xyyz

                     + pa_xyzz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xyyz_r_0(double fga,
                                     double fgb,
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
                                     double pb_xyy,
                                     double pb_xyyz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * pa_yz * fz * fgb

                     + 2.5 * fx * fx * fx * pa_yz * fz

                     + 5.0 * fx * fx * fx * pa_z * fz * pb_y

                     - 0.25 * fx * fx * pa_y * fz * fgb * pb_z

                     - 0.25 * fx * fx * pa_y * fz * fga * pb_z

                     - 0.5 * fx * fx * fz * fga * pb_yz

                     - pa_xyz * fx * pb_x * fz * fgb

                     - 0.5 * fx * pa_yzz * fz * fgb * pb_z

                     + 1.25 * fx * fx * fx * pa_y * fz * pb_z

                     + 2.5 * fx * fx * fx * fz * pb_yz

                     + 6.0 * pa_xyz * fz * fx * fx * pb_x

                     + 12.0 * pa_xz * fx * fx * fz * pb_xy

                     + 3.0 * fx * fx * pa_yzz * fz * pb_z

                     + 6.0 * fx * fx * pa_yz * fz * pb_yy

                     + 6.0 * fx * fx * pa_zz * fz * pb_yz

                     - 0.5 * pa_xy * fx * pb_xz * fz * fgb

                     - 0.5 * pa_xy * fz * fga * pb_xz * fx

                     - pa_x * fx * fz * fga * pb_xyz

                     - 0.5 * fx * pa_y * fz * fga * pb_yyz

                     - pa_xyzz * pb_xz * fz * fgb

                     + 3.0 * pa_xy * fz * fx * fx * pb_xz

                     + 6.0 * pa_x * fx * fx * fz * pb_xyz

                     + 3.0 * fx * fx * pa_y * fz * pb_yyz

                     + 7.0 * pa_xyzz * fz * pb_xz * fx

                     + 14.0 * pa_xyz * fz * fx * pb_xyy

                     + 14.0 * pa_xzz * fx * fz * pb_xyz

                     + 7.0 * fx * pa_yzz * fz * pb_yyz

                     - pa_xy * fz * fga * pb_xyyz

                     + 7.0 * pa_xy * fz * fx * pb_xyyz

                     + 16.0 * pa_xyzz * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xyzz_s_0(double fx,
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
                                     double pb_xyzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * pa_x * fx * fx * fx * pb_x

                     + 0.375 * fx * fx * fx * pa_y * pb_y

                     + 0.125 * fx * fx * fx * pa_zz

                     + 0.5 * fx * fx * fx * pa_z * pb_z

                     + 0.125 * fx * fx * fx * pb_zz

                     + 0.75 * pa_xy * fx * fx * pb_xy

                     + 0.25 * pa_xzz * fx * fx * pb_x

                     + pa_xz * fx * fx * pb_xz

                     + 0.25 * fx * fx * pa_yzz * pb_y

                     + fx * fx * pa_yz * pb_yz

                     + 0.25 * fx * fx * pa_zz * pb_zz

                     + 0.25 * pa_x * fx * fx * pb_xzz

                     + 0.25 * fx * fx * pa_y * pb_yzz

                     + 0.5 * pa_xyzz * pb_xy * fx

                     + 2.0 * pa_xyz * fx * pb_xyz

                     + 0.5 * pa_xzz * fx * pb_xzz

                     + 0.5 * fx * pa_yzz * pb_yzz

                     + 0.5 * pa_xy * fx * pb_xyzz

                     + pa_xyzz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xyzz_r_0(double fga,
                                     double fgb,
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
                                     double pb_xyzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (1.5 * fx * fx * fx * fx * fz

                     - 0.125 * fx * fx * fx * fz * fgb

                     - 0.125 * fx * fx * fx * fz * fga

                     - 0.25 * fx * fx * pa_zz * fz * fgb

                     + 3.75 * pa_x * fx * fx * fx * fz * pb_x

                     + 3.75 * fx * fx * fx * pa_y * fz * pb_y

                     + 1.25 * fx * fx * fx * pa_zz * fz

                     + 5.0 * fx * fx * fx * pa_z * fz * pb_z

                     - 0.25 * pa_x * fx * fx * pb_x * fz * fgb

                     - 0.25 * pa_x * fx * fx * fz * fga * pb_x

                     - 0.25 * fx * fx * pa_y * pb_y * fz * fgb

                     - 0.25 * fx * fx * pa_y * fz * fga * pb_y

                     - 0.25 * fx * fx * fz * fga * pb_zz

                     - 0.5 * pa_xzz * fx * pb_x * fz * fgb

                     - 0.5 * fx * pa_yzz * pb_y * fz * fgb

                     + 1.25 * fx * fx * fx * fz * pb_zz

                     + 9.0 * pa_xy * fx * fx * fz * pb_xy

                     + 3.0 * pa_xzz * fx * fx * fz * pb_x

                     + 12.0 * pa_xz * fx * fx * fz * pb_xz

                     + 3.0 * fx * fx * pa_yzz * fz * pb_y

                     + 12.0 * fx * fx * pa_yz * fz * pb_yz

                     + 3.0 * fx * fx * pa_zz * fz * pb_zz

                     - 0.5 * pa_xy * fx * pb_xy * fz * fgb

                     - 0.5 * pa_xy * fz * fga * pb_xy * fx

                     - 0.5 * pa_x * fx * fz * fga * pb_xzz

                     - 0.5 * fx * pa_y * fz * fga * pb_yzz

                     - pa_xyzz * pb_xy * fz * fgb

                     + 3.0 * pa_x * fx * fx * fz * pb_xzz

                     + 3.0 * fx * fx * pa_y * fz * pb_yzz

                     + 7.0 * pa_xyzz * fz * pb_xy * fx

                     + 28.0 * pa_xyz * fz * fx * pb_xyz

                     + 7.0 * pa_xzz * fx * fz * pb_xzz

                     + 7.0 * fx * pa_yzz * fz * pb_yzz

                     - pa_xy * fz * fga * pb_xyzz

                     + 7.0 * pa_xy * fz * fx * pb_xyzz

                     + 16.0 * pa_xyzz * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xzzz_s_0(double fx,
                                     double pa_xy,
                                     double pa_xyz,
                                     double pa_xyzz,
                                     double pa_y,
                                     double pa_yz,
                                     double pa_yzz,
                                     double pb_x,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_xzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * fx * pa_yz

                     + 1.125 * fx * fx * fx * pa_y * pb_z

                     + 1.5 * pa_xyz * fx * fx * pb_x

                     + 2.25 * pa_xy * fx * fx * pb_xz

                     + 0.75 * fx * fx * pa_yzz * pb_z

                     + 1.5 * fx * fx * pa_yz * pb_zz

                     + 0.25 * fx * fx * pa_y * pb_zzz

                     + 1.5 * pa_xyzz * pb_xz * fx

                     + 3.0 * pa_xyz * fx * pb_xzz

                     + 0.5 * fx * pa_yzz * pb_zzz

                     + 0.5 * pa_xy * fx * pb_xzzz

                     + pa_xyzz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xzzz_r_0(double fga,
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
                                     double pb_xzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_yz * fz * fgb

                     + 7.5 * fx * fx * fx * pa_yz * fz

                     + 11.25 * fx * fx * fx * pa_y * fz * pb_z

                     - 0.75 * fx * fx * pa_y * pb_z * fz * fgb

                     - 0.75 * fx * fx * pa_y * fz * fga * pb_z

                     - 3.0 * pa_xyz * fx * pb_x * fz * fgb

                     - 1.5 * fx * pa_yzz * pb_z * fz * fgb

                     + 18.0 * pa_xyz * fz * fx * fx * pb_x

                     + 27.0 * pa_xy * fx * fx * fz * pb_xz

                     + 9.0 * fx * fx * pa_yzz * fz * pb_z

                     + 18.0 * fx * fx * pa_yz * fz * pb_zz

                     - 1.5 * pa_xy * fx * pb_xz * fz * fgb

                     - 1.5 * pa_xy * fz * fga * pb_xz * fx

                     - 0.5 * fx * pa_y * fz * fga * pb_zzz

                     - 3.0 * pa_xyzz * pb_xz * fz * fgb

                     + 3.0 * fx * fx * pa_y * fz * pb_zzz

                     + 21.0 * pa_xyzz * fz * pb_xz * fx

                     + 42.0 * pa_xyz * fz * fx * pb_xzz

                     + 7.0 * fx * pa_yzz * fz * pb_zzz

                     - pa_xy * fz * fga * pb_xzzz

                     + 7.0 * pa_xy * fz * fx * pb_xzzz

                     + 16.0 * pa_xyzz * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_yyyy_s_0(double fx,
                                     double pa_x,
                                     double pa_xy,
                                     double pa_xyzz,
                                     double pa_xzz,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double pb_yyyy,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * pa_xy * fx * fx * fx

                     + 1.5 * pa_x * fx * fx * fx * pb_y

                     + 0.75 * pa_xyzz * fx * fx

                     + 3.0 * pa_xzz * fx * fx * pb_y

                     + 1.5 * pa_xy * fx * fx * pb_yy

                     + pa_x * fx * fx * pb_yyy

                     + 3.0 * pa_xyzz * pb_yy * fx

                     + 2.0 * pa_xzz * fx * pb_yyy

                     + 0.5 * pa_xy * fx * pb_yyyy

                     + pa_xyzz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_yyyy_r_0(double fga,
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
                                     double pb_yyyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xy * fx * fx * fz * fgb

                     - 0.75 * pa_xy * fz * fga * fx * fx

                     - 3.0 * pa_x * fx * fx * pb_y * fz * fgb

                     - 3.0 * pa_x * fx * fx * fz * fga * pb_y

                     - 3.0 * pa_xyzz * fx * fz * fgb

                     + 3.75 * pa_xy * fz * fx * fx * fx

                     - 6.0 * pa_xzz * fx * pb_y * fz * fgb

                     + 15.0 * pa_x * fx * fx * fx * fz * pb_y

                     + 9.0 * pa_xyzz * fz * fx * fx

                     + 36.0 * pa_xzz * fx * fx * fz * pb_y

                     - 3.0 * pa_xy * fx * pb_yy * fz * fgb

                     - 3.0 * pa_xy * fz * fga * pb_yy * fx

                     - 2.0 * pa_x * fx * fz * fga * pb_yyy

                     - 6.0 * pa_xyzz * pb_yy * fz * fgb

                     + 18.0 * pa_xy * fz * fx * fx * pb_yy

                     + 12.0 * pa_x * fx * fx * fz * pb_yyy

                     + 42.0 * pa_xyzz * fz * pb_yy * fx

                     + 28.0 * pa_xzz * fx * fz * pb_yyy

                     - pa_xy * fz * fga * pb_yyyy

                     + 7.0 * pa_xy * fz * fx * pb_yyyy

                     + 16.0 * pa_xyzz * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_yyyz_s_0(double fx,
                                     double pa_x,
                                     double pa_xy,
                                     double pa_xyz,
                                     double pa_xyzz,
                                     double pa_xz,
                                     double pa_xzz,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double pb_yyyz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xz * fx * fx * fx

                     + 0.375 * pa_x * fx * fx * fx * pb_z

                     + 1.5 * pa_xyz * fx * fx * pb_y

                     + 0.75 * pa_xzz * fx * fx * pb_z

                     + 1.5 * pa_xz * fx * fx * pb_yy

                     + 0.75 * pa_xy * fx * fx * pb_yz

                     + 0.75 * pa_x * fx * fx * pb_yyz

                     + 1.5 * pa_xyzz * pb_yz * fx

                     + pa_xyz * fx * pb_yyy

                     + 1.5 * pa_xzz * fx * pb_yyz

                     + 0.5 * pa_xy * fx * pb_yyyz

                     + pa_xyzz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_yyyz_r_0(double fga,
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
                                     double pb_yyy,
                                     double pb_yyyz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xz * fx * fx * fz * fgb

                     + 7.5 * pa_xz * fx * fx * fx * fz

                     - 0.75 * pa_x * fx * fx * fz * fgb * pb_z

                     - 0.75 * pa_x * fx * fx * fz * fga * pb_z

                     - 3.0 * pa_xyz * fx * pb_y * fz * fgb

                     - 1.5 * pa_xzz * fx * fz * fgb * pb_z

                     + 3.75 * pa_x * fx * fx * fx * fz * pb_z

                     + 18.0 * pa_xyz * fz * fx * fx * pb_y

                     + 9.0 * pa_xzz * fx * fx * fz * pb_z

                     + 18.0 * pa_xz * fx * fx * fz * pb_yy

                     - 1.5 * pa_xy * fx * pb_yz * fz * fgb

                     - 1.5 * pa_xy * fz * fga * pb_yz * fx

                     - 1.5 * pa_x * fx * fz * fga * pb_yyz

                     - 3.0 * pa_xyzz * pb_yz * fz * fgb

                     + 9.0 * pa_xy * fz * fx * fx * pb_yz

                     + 9.0 * pa_x * fx * fx * fz * pb_yyz

                     + 21.0 * pa_xyzz * fz * pb_yz * fx

                     + 14.0 * pa_xyz * fz * fx * pb_yyy

                     + 21.0 * pa_xzz * fx * fz * pb_yyz

                     - pa_xy * fz * fga * pb_yyyz

                     + 7.0 * pa_xy * fz * fx * pb_yyyz

                     + 16.0 * pa_xyzz * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_yyzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xy,
                                     double pa_xyz,
                                     double pa_xyzz,
                                     double pa_xz,
                                     double pa_xzz,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_yyzz,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * pa_xy * fx * fx * fx

                     + 0.75 * pa_x * fx * fx * fx * pb_y

                     + 0.25 * pa_xyzz * fx * fx

                     + pa_xyz * fx * fx * pb_z

                     + 0.75 * pa_xy * fx * fx * pb_yy

                     + 0.5 * pa_xzz * fx * fx * pb_y

                     + 2.0 * pa_xz * fx * fx * pb_yz

                     + 0.25 * pa_xy * fx * fx * pb_zz

                     + 0.5 * pa_x * fx * fx * pb_yzz

                     + 0.5 * pa_xyzz * pb_yy * fx

                     + 0.5 * pa_xyzz * fx * pb_zz

                     + 2.0 * pa_xyz * fx * pb_yyz

                     + pa_xzz * fx * pb_yzz

                     + 0.5 * pa_xy * fx * pb_yyzz

                     + pa_xyzz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_yyzz_r_0(double fga,
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
                                     double pb_yyzz,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- pa_xy * fx * fx * fz * fgb

                     + 3.75 * pa_xy * fx * fx * fx * fz

                     + 7.5 * pa_x * fx * fx * fx * fz * pb_y

                     - 0.25 * pa_xy * fz * fga * fx * fx

                     - 0.5 * pa_x * fx * fx * pb_y * fz * fgb

                     - 0.5 * pa_x * fx * fx * fz * fga * pb_y

                     - pa_xyzz * fx * fz * fgb

                     - 2.0 * pa_xyz * fx * fz * fgb * pb_z

                     - pa_xzz * fx * pb_y * fz * fgb

                     + 3.0 * pa_xyzz * fz * fx * fx

                     + 12.0 * pa_xyz * fz * fx * fx * pb_z

                     + 9.0 * pa_xy * fx * fx * fz * pb_yy

                     + 6.0 * pa_xzz * fx * fx * fz * pb_y

                     + 24.0 * pa_xz * fx * fx * fz * pb_yz

                     - 0.5 * pa_xy * fx * pb_yy * fz * fgb

                     - 0.5 * pa_xy * fx * fz * fgb * pb_zz

                     - 0.5 * pa_xy * fz * fga * pb_yy * fx

                     - 0.5 * pa_xy * fz * fga * fx * pb_zz

                     - pa_x * fx * fz * fga * pb_yzz

                     - pa_xyzz * pb_yy * fz * fgb

                     - pa_xyzz * fz * fgb * pb_zz

                     + 3.0 * pa_xy * fz * fx * fx * pb_zz

                     + 6.0 * pa_x * fx * fx * fz * pb_yzz

                     + 7.0 * pa_xyzz * fz * pb_yy * fx

                     + 7.0 * pa_xyzz * fz * fx * pb_zz

                     + 28.0 * pa_xyz * fz * fx * pb_yyz

                     + 14.0 * pa_xzz * fx * fz * pb_yzz

                     - pa_xy * fz * fga * pb_yyzz

                     + 7.0 * pa_xy * fz * fx * pb_yyzz

                     + 16.0 * pa_xyzz * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_yzzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xy,
                                     double pa_xyz,
                                     double pa_xyzz,
                                     double pa_xz,
                                     double pa_xzz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_yzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xz * fx * fx * fx

                     + 1.125 * pa_x * fx * fx * fx * pb_z

                     + 1.5 * pa_xyz * fx * fx * pb_y

                     + 2.25 * pa_xy * fx * fx * pb_yz

                     + 0.75 * pa_xzz * fx * fx * pb_z

                     + 1.5 * pa_xz * fx * fx * pb_zz

                     + 0.25 * pa_x * fx * fx * pb_zzz

                     + 1.5 * pa_xyzz * pb_yz * fx

                     + 3.0 * pa_xyz * fx * pb_yzz

                     + 0.5 * pa_xzz * fx * pb_zzz

                     + 0.5 * pa_xy * fx * pb_yzzz

                     + pa_xyzz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_yzzz_r_0(double fga,
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
                                     double pb_yzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xz * fx * fx * fz * fgb

                     + 7.5 * pa_xz * fx * fx * fx * fz

                     + 11.25 * pa_x * fx * fx * fx * fz * pb_z

                     - 0.75 * pa_x * fx * fx * pb_z * fz * fgb

                     - 0.75 * pa_x * fx * fx * fz * fga * pb_z

                     - 3.0 * pa_xyz * fx * pb_y * fz * fgb

                     - 1.5 * pa_xzz * fx * pb_z * fz * fgb

                     + 18.0 * pa_xyz * fz * fx * fx * pb_y

                     + 27.0 * pa_xy * fx * fx * fz * pb_yz

                     + 9.0 * pa_xzz * fx * fx * fz * pb_z

                     + 18.0 * pa_xz * fx * fx * fz * pb_zz

                     - 1.5 * pa_xy * fx * pb_yz * fz * fgb

                     - 1.5 * pa_xy * fz * fga * pb_yz * fx

                     - 0.5 * pa_x * fx * fz * fga * pb_zzz

                     - 3.0 * pa_xyzz * pb_yz * fz * fgb

                     + 3.0 * pa_x * fx * fx * fz * pb_zzz

                     + 21.0 * pa_xyzz * fz * pb_yz * fx

                     + 42.0 * pa_xyz * fz * fx * pb_yzz

                     + 7.0 * pa_xzz * fx * fz * pb_zzz

                     - pa_xy * fz * fga * pb_yzzz

                     + 7.0 * pa_xy * fz * fx * pb_yzzz

                     + 16.0 * pa_xyzz * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_zzzz_s_0(double fx,
                                     double pa_xy,
                                     double pa_xyz,
                                     double pa_xyzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double pb_zzzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.875 * pa_xy * fx * fx * fx

                     + 0.75 * pa_xyzz * fx * fx

                     + 6.0 * pa_xyz * fx * fx * pb_z

                     + 4.5 * pa_xy * fx * fx * pb_zz

                     + 3.0 * pa_xyzz * pb_zz * fx

                     + 4.0 * pa_xyz * fx * pb_zzz

                     + 0.5 * pa_xy * fx * pb_zzzz

                     + pa_xyzz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_zzzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xy,
                                     double pa_xyz,
                                     double pa_xyzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double pb_zzzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_xy * fx * fx * fz * fgb

                     + 18.75 * pa_xy * fx * fx * fx * fz

                     - 0.75 * pa_xy * fz * fga * fx * fx

                     - 3.0 * pa_xyzz * fx * fz * fgb

                     - 12.0 * pa_xyz * fx * pb_z * fz * fgb

                     + 9.0 * pa_xyzz * fz * fx * fx

                     + 72.0 * pa_xyz * fz * fx * fx * pb_z

                     + 54.0 * pa_xy * fx * fx * fz * pb_zz

                     - 3.0 * pa_xy * fx * pb_zz * fz * fgb

                     - 3.0 * pa_xy * fz * fga * pb_zz * fx

                     - 6.0 * pa_xyzz * pb_zz * fz * fgb

                     + 42.0 * pa_xyzz * fz * pb_zz * fx

                     + 56.0 * pa_xyz * fz * fx * pb_zzz

                     - pa_xy * fz * fga * pb_zzzz

                     + 7.0 * pa_xy * fz * fx * pb_zzzz

                     + 16.0 * pa_xyzz * fz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xxxx_s_0(double fx,
                                     double pa_xz,
                                     double pa_xzzz,
                                     double pa_z,
                                     double pa_zzz,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxx,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_xz * fx * fx * fx

                     + 4.5 * fx * fx * fx * pa_z * pb_x

                     + 0.75 * pa_xzzz * fx * fx

                     + 3.0 * fx * fx * pa_zzz * pb_x

                     + 4.5 * pa_xz * fx * fx * pb_xx

                     + 3.0 * fx * fx * pa_z * pb_xxx

                     + 3.0 * pa_xzzz * pb_xx * fx

                     + 2.0 * fx * pa_zzz * pb_xxx

                     + 1.5 * pa_xz * fx * pb_xxxx

                     + pa_xzzz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xxxx_r_0(double fga,
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
                                     double pb_xxxx,
                                     double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_xz * fx * fx * fz * fgb

                     - 2.25 * pa_xz * fz * fga * fx * fx

                     - 9.0 * fx * fx * pa_z * pb_x * fz * fgb

                     - 9.0 * fx * fx * pa_z * fz * fga * pb_x

                     - 3.0 * pa_xzzz * fx * fz * fgb

                     + 11.25 * pa_xz * fz * fx * fx * fx

                     - 6.0 * fx * pa_zzz * pb_x * fz * fgb

                     + 45.0 * fx * fx * fx * pa_z * fz * pb_x

                     + 9.0 * pa_xzzz * fz * fx * fx

                     + 36.0 * fx * fx * pa_zzz * fz * pb_x

                     - 9.0 * pa_xz * fx * pb_xx * fz * fgb

                     - 9.0 * pa_xz * fz * fga * pb_xx * fx

                     - 6.0 * fx * pa_z * fz * fga * pb_xxx

                     - 6.0 * pa_xzzz * pb_xx * fz * fgb

                     + 54.0 * pa_xz * fz * fx * fx * pb_xx

                     + 36.0 * fx * fx * pa_z * fz * pb_xxx

                     + 42.0 * pa_xzzz * fz * pb_xx * fx

                     + 28.0 * fx * pa_zzz * fz * pb_xxx

                     - 3.0 * pa_xz * fz * fga * pb_xxxx

                     + 21.0 * pa_xz * fz * fx * pb_xxxx

                     + 16.0 * pa_xzzz * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xxxy_s_0(double fx,
                                     double pa_xz,
                                     double pa_xzzz,
                                     double pa_z,
                                     double pa_zzz,
                                     double pb_xxxy,
                                     double pb_xxy,
                                     double pb_xy,
                                     double pb_y,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pa_z * pb_y

                     + 0.75 * fx * fx * pa_zzz * pb_y

                     + 2.25 * pa_xz * fx * fx * pb_xy

                     + 2.25 * fx * fx * pa_z * pb_xxy

                     + 1.5 * pa_xzzz * pb_xy * fx

                     + 1.5 * fx * pa_zzz * pb_xxy

                     + 1.5 * pa_xz * fx * pb_xxxy

                     + pa_xzzz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xxxy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xz,
                                     double pa_xzzz,
                                     double pa_z,
                                     double pa_zzz,
                                     double pb_xxxy,
                                     double pb_xxy,
                                     double pb_xy,
                                     double pb_y,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pa_z * fz * fgb * pb_y

                     - 2.25 * fx * fx * pa_z * fz * fga * pb_y

                     - 1.5 * fx * pa_zzz * fz * fgb * pb_y

                     + 11.25 * fx * fx * fx * pa_z * fz * pb_y

                     + 9.0 * fx * fx * pa_zzz * fz * pb_y

                     - 4.5 * pa_xz * fx * pb_xy * fz * fgb

                     - 4.5 * pa_xz * fz * fga * pb_xy * fx

                     - 4.5 * fx * pa_z * fz * fga * pb_xxy

                     - 3.0 * pa_xzzz * pb_xy * fz * fgb

                     + 27.0 * pa_xz * fz * fx * fx * pb_xy

                     + 27.0 * fx * fx * pa_z * fz * pb_xxy

                     + 21.0 * pa_xzzz * fz * pb_xy * fx

                     + 21.0 * fx * pa_zzz * fz * pb_xxy

                     - 3.0 * pa_xz * fz * fga * pb_xxxy

                     + 21.0 * pa_xz * fz * fx * pb_xxxy

                     + 16.0 * pa_xzzz * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xxxz_s_0(double fx,
                                     double pa_x,
                                     double pa_xz,
                                     double pa_xzz,
                                     double pa_xzzz,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxx,
                                     double pb_xxxz,
                                     double pb_xxz,
                                     double pb_xz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.5625 * fx * fx * fx * fx

                     + 1.125 * fx * fx * fx * pa_zz

                     + 1.125 * pa_x * fx * fx * fx * pb_x

                     + 1.125 * fx * fx * fx * pa_z * pb_z

                     + 1.125 * fx * fx * fx * pb_xx

                     + 2.25 * pa_xzz * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_zzz * pb_z

                     + 2.25 * fx * fx * pa_zz * pb_xx

                     + 2.25 * pa_xz * fx * fx * pb_xz

                     + 0.75 * pa_x * fx * fx * pb_xxx

                     + 2.25 * fx * fx * pa_z * pb_xxz

                     + 1.5 * pa_xzzz * pb_xz * fx

                     + 1.5 * pa_xzz * fx * pb_xxx

                     + 1.5 * fx * pa_zzz * pb_xxz

                     + 1.5 * pa_xz * fx * pb_xxxz

                     + pa_xzzz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xxxz_r_0(double fga,
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
                                     double pb_xxx,
                                     double pb_xxxz,
                                     double pb_xxz,
                                     double pb_xz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.125 * fx * fx * fx * fz * fgb

                     - 1.125 * fx * fx * fx * fz * fga

                     - 2.25 * fx * fx * pa_zz * fz * fgb

                     + 4.5 * fx * fx * fx * fx * fz

                     + 11.25 * fx * fx * fx * pa_zz * fz

                     - 2.25 * pa_x * fx * fx * pb_x * fz * fgb

                     - 2.25 * pa_x * fx * fx * fz * fga * pb_x

                     - 2.25 * fx * fx * pa_z * fz * fgb * pb_z

                     - 2.25 * fx * fx * pa_z * fz * fga * pb_z

                     - 2.25 * fx * fx * fz * fga * pb_xx

                     - 4.5 * pa_xzz * fx * pb_x * fz * fgb

                     + 11.25 * pa_x * fx * fx * fx * fz * pb_x

                     - 1.5 * fx * pa_zzz * fz * fgb * pb_z

                     + 11.25 * fx * fx * fx * pa_z * fz * pb_z

                     + 11.25 * fx * fx * fx * fz * pb_xx

                     + 27.0 * pa_xzz * fz * fx * fx * pb_x

                     + 9.0 * fx * fx * pa_zzz * fz * pb_z

                     + 27.0 * fx * fx * pa_zz * fz * pb_xx

                     - 4.5 * pa_xz * fx * pb_xz * fz * fgb

                     - 4.5 * pa_xz * fz * fga * pb_xz * fx

                     - 1.5 * pa_x * fx * fz * fga * pb_xxx

                     - 4.5 * fx * pa_z * fz * fga * pb_xxz

                     - 3.0 * pa_xzzz * pb_xz * fz * fgb

                     + 27.0 * pa_xz * fz * fx * fx * pb_xz

                     + 9.0 * pa_x * fx * fx * fz * pb_xxx

                     + 27.0 * fx * fx * pa_z * fz * pb_xxz

                     + 21.0 * pa_xzzz * fz * pb_xz * fx

                     + 21.0 * pa_xzz * fz * fx * pb_xxx

                     + 21.0 * fx * pa_zzz * fz * pb_xxz

                     - 3.0 * pa_xz * fz * fga * pb_xxxz

                     + 21.0 * pa_xz * fz * fx * pb_xxxz

                     + 16.0 * pa_xzzz * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xxyy_s_0(double fx,
                                     double pa_xz,
                                     double pa_xzzz,
                                     double pa_z,
                                     double pa_zzz,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxyy,
                                     double pb_xyy,
                                     double pb_yy,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * pa_xz * fx * fx * fx

                     + 0.75 * fx * fx * fx * pa_z * pb_x

                     + 0.25 * pa_xzzz * fx * fx

                     + 0.5 * fx * fx * pa_zzz * pb_x

                     + 0.75 * pa_xz * fx * fx * pb_xx

                     + 0.75 * pa_xz * fx * fx * pb_yy

                     + 1.5 * fx * fx * pa_z * pb_xyy

                     + 0.5 * pa_xzzz * pb_xx * fx

                     + 0.5 * pa_xzzz * fx * pb_yy

                     + fx * pa_zzz * pb_xyy

                     + 1.5 * pa_xz * fx * pb_xxyy

                     + pa_xzzz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xxyy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xz,
                                     double pa_xzzz,
                                     double pa_z,
                                     double pa_zzz,
                                     double pb_x,
                                     double pb_xx,
                                     double pb_xxyy,
                                     double pb_xyy,
                                     double pb_yy,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xz * fx * fx * fz * fgb

                     - 0.75 * pa_xz * fz * fga * fx * fx

                     - 1.5 * fx * fx * pa_z * pb_x * fz * fgb

                     - 1.5 * fx * fx * pa_z * fz * fga * pb_x

                     - pa_xzzz * fx * fz * fgb

                     + 3.75 * pa_xz * fz * fx * fx * fx

                     - fx * pa_zzz * pb_x * fz * fgb

                     + 7.5 * fx * fx * fx * pa_z * fz * pb_x

                     + 3.0 * pa_xzzz * fz * fx * fx

                     + 6.0 * fx * fx * pa_zzz * fz * pb_x

                     - 1.5 * pa_xz * fx * pb_xx * fz * fgb

                     - 1.5 * pa_xz * fx * fz * fgb * pb_yy

                     - 1.5 * pa_xz * fz * fga * pb_xx * fx

                     - 1.5 * pa_xz * fz * fga * fx * pb_yy

                     - 3.0 * fx * pa_z * fz * fga * pb_xyy

                     - pa_xzzz * pb_xx * fz * fgb

                     - pa_xzzz * fz * fgb * pb_yy

                     + 9.0 * pa_xz * fz * fx * fx * pb_xx

                     + 9.0 * pa_xz * fz * fx * fx * pb_yy

                     + 18.0 * fx * fx * pa_z * fz * pb_xyy

                     + 7.0 * pa_xzzz * fz * pb_xx * fx

                     + 7.0 * pa_xzzz * fz * fx * pb_yy

                     + 14.0 * fx * pa_zzz * fz * pb_xyy

                     - 3.0 * pa_xz * fz * fga * pb_xxyy

                     + 21.0 * pa_xz * fz * fx * pb_xxyy

                     + 16.0 * pa_xzzz * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xxyz_s_0(double fx,
                                     double pa_x,
                                     double pa_xz,
                                     double pa_xzz,
                                     double pa_xzzz,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pb_xxy,
                                     double pb_xxyz,
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_y,
                                     double pb_yz,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * pa_x * fx * fx * fx * pb_y

                     + 0.75 * fx * fx * fx * pb_xy

                     + 0.75 * pa_xzz * fx * fx * pb_y

                     + 1.5 * fx * fx * pa_zz * pb_xy

                     + 0.75 * pa_xz * fx * fx * pb_yz

                     + 0.75 * pa_x * fx * fx * pb_xxy

                     + 1.5 * fx * fx * pa_z * pb_xyz

                     + 0.5 * pa_xzzz * fx * pb_yz

                     + 1.5 * pa_xzz * fx * pb_xxy

                     + fx * pa_zzz * pb_xyz

                     + 1.5 * pa_xz * fx * pb_xxyz

                     + pa_xzzz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xxyz_r_0(double fga,
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
                                     double pb_xxy,
                                     double pb_xxyz,
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_y,
                                     double pb_yz,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.75 * pa_x * fx * fx * fz * fgb * pb_y

                     - 0.75 * pa_x * fx * fx * fz * fga * pb_y

                     - 1.5 * fx * fx * fz * fga * pb_xy

                     - 1.5 * pa_xzz * fx * fz * fgb * pb_y

                     + 3.75 * pa_x * fx * fx * fx * fz * pb_y

                     + 7.5 * fx * fx * fx * fz * pb_xy

                     + 9.0 * pa_xzz * fz * fx * fx * pb_y

                     + 18.0 * fx * fx * pa_zz * fz * pb_xy

                     - 1.5 * pa_xz * fx * fz * fgb * pb_yz

                     - 1.5 * pa_xz * fz * fga * fx * pb_yz

                     - 1.5 * pa_x * fx * fz * fga * pb_xxy

                     - 3.0 * fx * pa_z * fz * fga * pb_xyz

                     - pa_xzzz * fz * fgb * pb_yz

                     + 9.0 * pa_xz * fz * fx * fx * pb_yz

                     + 9.0 * pa_x * fx * fx * fz * pb_xxy

                     + 18.0 * fx * fx * pa_z * fz * pb_xyz

                     + 7.0 * pa_xzzz * fz * fx * pb_yz

                     + 21.0 * pa_xzz * fz * fx * pb_xxy

                     + 14.0 * fx * pa_zzz * fz * pb_xyz

                     - 3.0 * pa_xz * fz * fga * pb_xxyz

                     + 21.0 * pa_xz * fz * fx * pb_xxyz

                     + 16.0 * pa_xzzz * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xxzz_s_0(double fx,
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
                                     double pb_xxzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_xz * fx * fx * fx

                     + 2.25 * fx * fx * fx * pa_z * pb_x

                     + 0.75 * pa_x * fx * fx * fx * pb_z

                     + 1.5 * fx * fx * fx * pb_xz

                     + 0.25 * pa_xzzz * fx * fx

                     + 1.5 * pa_xzz * fx * fx * pb_z

                     + 2.25 * pa_xz * fx * fx * pb_xx

                     + 0.5 * fx * fx * pa_zzz * pb_x

                     + 3.0 * fx * fx * pa_zz * pb_xz

                     + 0.75 * pa_xz * fx * fx * pb_zz

                     + 1.5 * pa_x * fx * fx * pb_xxz

                     + 1.5 * fx * fx * pa_z * pb_xzz

                     + 0.5 * pa_xzzz * pb_xx * fx

                     + 0.5 * pa_xzzz * fx * pb_zz

                     + 3.0 * pa_xzz * fx * pb_xxz

                     + fx * pa_zzz * pb_xzz

                     + 1.5 * pa_xz * fx * pb_xxzz

                     + pa_xzzz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xxzz_r_0(double fga,
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
                                     double pb_xxzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xz * fx * fx * fz * fgb

                     + 11.25 * pa_xz * fx * fx * fx * fz

                     + 22.5 * fx * fx * fx * pa_z * fz * pb_x

                     - 0.75 * pa_xz * fz * fga * fx * fx

                     - 1.5 * pa_x * fx * fx * fz * fgb * pb_z

                     - 1.5 * pa_x * fx * fx * fz * fga * pb_z

                     - 1.5 * fx * fx * pa_z * pb_x * fz * fgb

                     - 1.5 * fx * fx * pa_z * fz * fga * pb_x

                     - 3.0 * fx * fx * fz * fga * pb_xz

                     - pa_xzzz * fx * fz * fgb

                     - 3.0 * pa_xzz * fx * fz * fgb * pb_z

                     + 7.5 * pa_x * fx * fx * fx * fz * pb_z

                     - fx * pa_zzz * pb_x * fz * fgb

                     + 15.0 * fx * fx * fx * fz * pb_xz

                     + 3.0 * pa_xzzz * fz * fx * fx

                     + 18.0 * pa_xzz * fz * fx * fx * pb_z

                     + 27.0 * pa_xz * fx * fx * fz * pb_xx

                     + 6.0 * fx * fx * pa_zzz * fz * pb_x

                     + 36.0 * fx * fx * pa_zz * fz * pb_xz

                     - 1.5 * pa_xz * fx * pb_xx * fz * fgb

                     - 1.5 * pa_xz * fx * fz * fgb * pb_zz

                     - 1.5 * pa_xz * fz * fga * pb_xx * fx

                     - 1.5 * pa_xz * fz * fga * fx * pb_zz

                     - 3.0 * pa_x * fx * fz * fga * pb_xxz

                     - 3.0 * fx * pa_z * fz * fga * pb_xzz

                     - pa_xzzz * pb_xx * fz * fgb

                     - pa_xzzz * fz * fgb * pb_zz

                     + 9.0 * pa_xz * fz * fx * fx * pb_zz

                     + 18.0 * pa_x * fx * fx * fz * pb_xxz

                     + 18.0 * fx * fx * pa_z * fz * pb_xzz

                     + 7.0 * pa_xzzz * fz * pb_xx * fx

                     + 7.0 * pa_xzzz * fz * fx * pb_zz

                     + 42.0 * pa_xzz * fz * fx * pb_xxz

                     + 14.0 * fx * pa_zzz * fz * pb_xzz

                     - 3.0 * pa_xz * fz * fga * pb_xxzz

                     + 21.0 * pa_xz * fz * fx * pb_xxzz

                     + 16.0 * pa_xzzz * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xyyy_s_0(double fx,
                                     double pa_xz,
                                     double pa_xzzz,
                                     double pa_z,
                                     double pa_zzz,
                                     double pb_xy,
                                     double pb_xyyy,
                                     double pb_y,
                                     double pb_yyy,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pa_z * pb_y

                     + 0.75 * fx * fx * pa_zzz * pb_y

                     + 2.25 * pa_xz * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_z * pb_yyy

                     + 1.5 * pa_xzzz * pb_xy * fx

                     + 0.5 * fx * pa_zzz * pb_yyy

                     + 1.5 * pa_xz * fx * pb_xyyy

                     + pa_xzzz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xyyy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xz,
                                     double pa_xzzz,
                                     double pa_z,
                                     double pa_zzz,
                                     double pb_xy,
                                     double pb_xyyy,
                                     double pb_y,
                                     double pb_yyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pa_z * pb_y * fz * fgb

                     - 2.25 * fx * fx * pa_z * fz * fga * pb_y

                     - 1.5 * fx * pa_zzz * pb_y * fz * fgb

                     + 11.25 * fx * fx * fx * pa_z * fz * pb_y

                     + 9.0 * fx * fx * pa_zzz * fz * pb_y

                     - 4.5 * pa_xz * fx * pb_xy * fz * fgb

                     - 4.5 * pa_xz * fz * fga * pb_xy * fx

                     - 1.5 * fx * pa_z * fz * fga * pb_yyy

                     - 3.0 * pa_xzzz * pb_xy * fz * fgb

                     + 27.0 * pa_xz * fz * fx * fx * pb_xy

                     + 9.0 * fx * fx * pa_z * fz * pb_yyy

                     + 21.0 * pa_xzzz * fz * pb_xy * fx

                     + 7.0 * fx * pa_zzz * fz * pb_yyy

                     - 3.0 * pa_xz * fz * fga * pb_xyyy

                     + 21.0 * pa_xz * fz * fx * pb_xyyy

                     + 16.0 * pa_xzzz * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xyyz_s_0(double fx,
                                     double pa_x,
                                     double pa_xz,
                                     double pa_xzz,
                                     double pa_xzzz,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pb_x,
                                     double pb_xyy,
                                     double pb_xyyz,
                                     double pb_xz,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * fx * fx * fx * pa_zz

                     + 0.375 * pa_x * fx * fx * fx * pb_x

                     + 0.375 * fx * fx * fx * pa_z * pb_z

                     + 0.375 * fx * fx * fx * pb_yy

                     + 0.75 * pa_xzz * fx * fx * pb_x

                     + 0.25 * fx * fx * pa_zzz * pb_z

                     + 0.75 * fx * fx * pa_zz * pb_yy

                     + 0.75 * pa_xz * fx * fx * pb_xz

                     + 0.75 * pa_x * fx * fx * pb_xyy

                     + 0.75 * fx * fx * pa_z * pb_yyz

                     + 0.5 * pa_xzzz * pb_xz * fx

                     + 1.5 * pa_xzz * fx * pb_xyy

                     + 0.5 * fx * pa_zzz * pb_yyz

                     + 1.5 * pa_xz * fx * pb_xyyz

                     + pa_xzzz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xyyz_r_0(double fga,
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
                                     double pb_xyy,
                                     double pb_xyyz,
                                     double pb_xz,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.375 * fx * fx * fx * fz * fgb

                     - 0.375 * fx * fx * fx * fz * fga

                     - 0.75 * fx * fx * pa_zz * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     + 3.75 * fx * fx * fx * pa_zz * fz

                     - 0.75 * pa_x * fx * fx * pb_x * fz * fgb

                     - 0.75 * pa_x * fx * fx * fz * fga * pb_x

                     - 0.75 * fx * fx * pa_z * fz * fgb * pb_z

                     - 0.75 * fx * fx * pa_z * fz * fga * pb_z

                     - 0.75 * fx * fx * fz * fga * pb_yy

                     - 1.5 * pa_xzz * fx * pb_x * fz * fgb

                     + 3.75 * pa_x * fx * fx * fx * fz * pb_x

                     - 0.5 * fx * pa_zzz * fz * fgb * pb_z

                     + 3.75 * fx * fx * fx * pa_z * fz * pb_z

                     + 3.75 * fx * fx * fx * fz * pb_yy

                     + 9.0 * pa_xzz * fz * fx * fx * pb_x

                     + 3.0 * fx * fx * pa_zzz * fz * pb_z

                     + 9.0 * fx * fx * pa_zz * fz * pb_yy

                     - 1.5 * pa_xz * fx * pb_xz * fz * fgb

                     - 1.5 * pa_xz * fz * fga * pb_xz * fx

                     - 1.5 * pa_x * fx * fz * fga * pb_xyy

                     - 1.5 * fx * pa_z * fz * fga * pb_yyz

                     - pa_xzzz * pb_xz * fz * fgb

                     + 9.0 * pa_xz * fz * fx * fx * pb_xz

                     + 9.0 * pa_x * fx * fx * fz * pb_xyy

                     + 9.0 * fx * fx * pa_z * fz * pb_yyz

                     + 7.0 * pa_xzzz * fz * pb_xz * fx

                     + 21.0 * pa_xzz * fz * fx * pb_xyy

                     + 7.0 * fx * pa_zzz * fz * pb_yyz

                     - 3.0 * pa_xz * fz * fga * pb_xyyz

                     + 21.0 * pa_xz * fz * fx * pb_xyyz

                     + 16.0 * pa_xzzz * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xyzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xz,
                                     double pa_xzz,
                                     double pa_xzzz,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_xyzz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_yzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pa_z * pb_y

                     + 0.75 * fx * fx * fx * pb_yz

                     + 2.25 * pa_xz * fx * fx * pb_xy

                     + 0.25 * fx * fx * pa_zzz * pb_y

                     + 1.5 * fx * fx * pa_zz * pb_yz

                     + 1.5 * pa_x * fx * fx * pb_xyz

                     + 0.75 * fx * fx * pa_z * pb_yzz

                     + 0.5 * pa_xzzz * pb_xy * fx

                     + 3.0 * pa_xzz * fx * pb_xyz

                     + 0.5 * fx * pa_zzz * pb_yzz

                     + 1.5 * pa_xz * fx * pb_xyzz

                     + pa_xzzz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xyzz_r_0(double fga,
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
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_xyzz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_yzz,
                                     double r_0_0)
    {
        return r_0_0 * (11.25 * fx * fx * fx * pa_z * fz * pb_y

                     - 0.75 * fx * fx * pa_z * pb_y * fz * fgb

                     - 0.75 * fx * fx * pa_z * fz * fga * pb_y

                     - 1.5 * fx * fx * fz * fga * pb_yz

                     - 0.5 * fx * pa_zzz * pb_y * fz * fgb

                     + 7.5 * fx * fx * fx * fz * pb_yz

                     + 27.0 * pa_xz * fx * fx * fz * pb_xy

                     + 3.0 * fx * fx * pa_zzz * fz * pb_y

                     + 18.0 * fx * fx * pa_zz * fz * pb_yz

                     - 1.5 * pa_xz * fx * pb_xy * fz * fgb

                     - 1.5 * pa_xz * fz * fga * pb_xy * fx

                     - 3.0 * pa_x * fx * fz * fga * pb_xyz

                     - 1.5 * fx * pa_z * fz * fga * pb_yzz

                     - pa_xzzz * pb_xy * fz * fgb

                     + 18.0 * pa_x * fx * fx * fz * pb_xyz

                     + 9.0 * fx * fx * pa_z * fz * pb_yzz

                     + 7.0 * pa_xzzz * fz * pb_xy * fx

                     + 42.0 * pa_xzz * fz * fx * pb_xyz

                     + 7.0 * fx * pa_zzz * fz * pb_yzz

                     - 3.0 * pa_xz * fz * fga * pb_xyzz

                     + 21.0 * pa_xz * fz * fx * pb_xyzz

                     + 16.0 * pa_xzzz * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xzzz_s_0(double fx,
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
                                     double pb_xzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.9375 * fx * fx * fx * fx

                     + 1.875 * pa_x * fx * fx * fx * pb_x

                     + 1.125 * fx * fx * fx * pa_zz

                     + 3.375 * fx * fx * fx * pa_z * pb_z

                     + 1.125 * fx * fx * fx * pb_zz

                     + 2.25 * pa_xzz * fx * fx * pb_x

                     + 6.75 * pa_xz * fx * fx * pb_xz

                     + 0.75 * fx * fx * pa_zzz * pb_z

                     + 2.25 * fx * fx * pa_zz * pb_zz

                     + 2.25 * pa_x * fx * fx * pb_xzz

                     + 0.75 * fx * fx * pa_z * pb_zzz

                     + 1.5 * pa_xzzz * pb_xz * fx

                     + 4.5 * pa_xzz * fx * pb_xzz

                     + 0.5 * fx * pa_zzz * pb_zzz

                     + 1.5 * pa_xz * fx * pb_xzzz

                     + pa_xzzz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xzzz_r_0(double fga,
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
                                     double pb_xzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double r_0_0)
    {
        return r_0_0 * (7.5 * fx * fx * fx * fx * fz

                     - 1.125 * fx * fx * fx * fz * fgb

                     - 1.125 * fx * fx * fx * fz * fga

                     - 2.25 * fx * fx * pa_zz * fz * fgb

                     + 18.75 * pa_x * fx * fx * fx * fz * pb_x

                     + 11.25 * fx * fx * fx * pa_zz * fz

                     + 33.75 * fx * fx * fx * pa_z * fz * pb_z

                     - 2.25 * pa_x * fx * fx * pb_x * fz * fgb

                     - 2.25 * pa_x * fx * fx * fz * fga * pb_x

                     - 2.25 * fx * fx * pa_z * pb_z * fz * fgb

                     - 2.25 * fx * fx * pa_z * fz * fga * pb_z

                     - 2.25 * fx * fx * fz * fga * pb_zz

                     - 4.5 * pa_xzz * fx * pb_x * fz * fgb

                     - 1.5 * fx * pa_zzz * pb_z * fz * fgb

                     + 11.25 * fx * fx * fx * fz * pb_zz

                     + 27.0 * pa_xzz * fz * fx * fx * pb_x

                     + 81.0 * pa_xz * fx * fx * fz * pb_xz

                     + 9.0 * fx * fx * pa_zzz * fz * pb_z

                     + 27.0 * fx * fx * pa_zz * fz * pb_zz

                     - 4.5 * pa_xz * fx * pb_xz * fz * fgb

                     - 4.5 * pa_xz * fz * fga * pb_xz * fx

                     - 4.5 * pa_x * fx * fz * fga * pb_xzz

                     - 1.5 * fx * pa_z * fz * fga * pb_zzz

                     - 3.0 * pa_xzzz * pb_xz * fz * fgb

                     + 27.0 * pa_x * fx * fx * fz * pb_xzz

                     + 9.0 * fx * fx * pa_z * fz * pb_zzz

                     + 21.0 * pa_xzzz * fz * pb_xz * fx

                     + 63.0 * pa_xzz * fz * fx * pb_xzz

                     + 7.0 * fx * pa_zzz * fz * pb_zzz

                     - 3.0 * pa_xz * fz * fga * pb_xzzz

                     + 21.0 * pa_xz * fz * fx * pb_xzzz

                     + 16.0 * pa_xzzz * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_yyyy_s_0(double fx,
                                     double pa_xz,
                                     double pa_xzzz,
                                     double pb_yy,
                                     double pb_yyyy,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_xz * fx * fx * fx

                     + 0.75 * pa_xzzz * fx * fx

                     + 4.5 * pa_xz * fx * fx * pb_yy

                     + 3.0 * pa_xzzz * pb_yy * fx

                     + 1.5 * pa_xz * fx * pb_yyyy

                     + pa_xzzz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_yyyy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_xz,
                                     double pa_xzzz,
                                     double pb_yy,
                                     double pb_yyyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_xz * fx * fx * fz * fgb

                     - 2.25 * pa_xz * fz * fga * fx * fx

                     - 3.0 * pa_xzzz * fx * fz * fgb

                     + 11.25 * pa_xz * fz * fx * fx * fx

                     + 9.0 * pa_xzzz * fz * fx * fx

                     - 9.0 * pa_xz * fx * pb_yy * fz * fgb

                     - 9.0 * pa_xz * fz * fga * pb_yy * fx

                     - 6.0 * pa_xzzz * pb_yy * fz * fgb

                     + 54.0 * pa_xz * fz * fx * fx * pb_yy

                     + 42.0 * pa_xzzz * fz * pb_yy * fx

                     - 3.0 * pa_xz * fz * fga * pb_yyyy

                     + 21.0 * pa_xz * fz * fx * pb_yyyy

                     + 16.0 * pa_xzzz * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_yyyz_s_0(double fx,
                                     double pa_x,
                                     double pa_xz,
                                     double pa_xzz,
                                     double pa_xzzz,
                                     double pb_y,
                                     double pb_yyy,
                                     double pb_yyyz,
                                     double pb_yz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_x * fx * fx * fx * pb_y

                     + 2.25 * pa_xzz * fx * fx * pb_y

                     + 2.25 * pa_xz * fx * fx * pb_yz

                     + 0.75 * pa_x * fx * fx * pb_yyy

                     + 1.5 * pa_xzzz * pb_yz * fx

                     + 1.5 * pa_xzz * fx * pb_yyy

                     + 1.5 * pa_xz * fx * pb_yyyz

                     + pa_xzzz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_yyyz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xz,
                                     double pa_xzz,
                                     double pa_xzzz,
                                     double pb_y,
                                     double pb_yyy,
                                     double pb_yyyz,
                                     double pb_yz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * pa_x * fx * fx * pb_y * fz * fgb

                     - 2.25 * pa_x * fx * fx * fz * fga * pb_y

                     - 4.5 * pa_xzz * fx * pb_y * fz * fgb

                     + 11.25 * pa_x * fx * fx * fx * fz * pb_y

                     + 27.0 * pa_xzz * fz * fx * fx * pb_y

                     - 4.5 * pa_xz * fx * pb_yz * fz * fgb

                     - 4.5 * pa_xz * fz * fga * pb_yz * fx

                     - 1.5 * pa_x * fx * fz * fga * pb_yyy

                     - 3.0 * pa_xzzz * pb_yz * fz * fgb

                     + 27.0 * pa_xz * fz * fx * fx * pb_yz

                     + 9.0 * pa_x * fx * fx * fz * pb_yyy

                     + 21.0 * pa_xzzz * fz * pb_yz * fx

                     + 21.0 * pa_xzz * fz * fx * pb_yyy

                     - 3.0 * pa_xz * fz * fga * pb_yyyz

                     + 21.0 * pa_xz * fz * fx * pb_yyyz

                     + 16.0 * pa_xzzz * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_yyzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xz,
                                     double pa_xzz,
                                     double pa_xzzz,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_yyzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_xz * fx * fx * fx

                     + 0.75 * pa_x * fx * fx * fx * pb_z

                     + 0.25 * pa_xzzz * fx * fx

                     + 1.5 * pa_xzz * fx * fx * pb_z

                     + 2.25 * pa_xz * fx * fx * pb_yy

                     + 0.75 * pa_xz * fx * fx * pb_zz

                     + 1.5 * pa_x * fx * fx * pb_yyz

                     + 0.5 * pa_xzzz * pb_yy * fx

                     + 0.5 * pa_xzzz * fx * pb_zz

                     + 3.0 * pa_xzz * fx * pb_yyz

                     + 1.5 * pa_xz * fx * pb_yyzz

                     + pa_xzzz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_yyzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_x,
                                     double pa_xz,
                                     double pa_xzz,
                                     double pa_xzzz,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_yyzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xz * fx * fx * fz * fgb

                     + 11.25 * pa_xz * fx * fx * fx * fz

                     - 0.75 * pa_xz * fz * fga * fx * fx

                     - 1.5 * pa_x * fx * fx * fz * fgb * pb_z

                     - 1.5 * pa_x * fx * fx * fz * fga * pb_z

                     - pa_xzzz * fx * fz * fgb

                     - 3.0 * pa_xzz * fx * fz * fgb * pb_z

                     + 7.5 * pa_x * fx * fx * fx * fz * pb_z

                     + 3.0 * pa_xzzz * fz * fx * fx

                     + 18.0 * pa_xzz * fz * fx * fx * pb_z

                     + 27.0 * pa_xz * fx * fx * fz * pb_yy

                     - 1.5 * pa_xz * fx * pb_yy * fz * fgb

                     - 1.5 * pa_xz * fx * fz * fgb * pb_zz

                     - 1.5 * pa_xz * fz * fga * pb_yy * fx

                     - 1.5 * pa_xz * fz * fga * fx * pb_zz

                     - 3.0 * pa_x * fx * fz * fga * pb_yyz

                     - pa_xzzz * pb_yy * fz * fgb

                     - pa_xzzz * fz * fgb * pb_zz

                     + 9.0 * pa_xz * fz * fx * fx * pb_zz

                     + 18.0 * pa_x * fx * fx * fz * pb_yyz

                     + 7.0 * pa_xzzz * fz * pb_yy * fx

                     + 7.0 * pa_xzzz * fz * fx * pb_zz

                     + 42.0 * pa_xzz * fz * fx * pb_yyz

                     - 3.0 * pa_xz * fz * fga * pb_yyzz

                     + 21.0 * pa_xz * fz * fx * pb_yyzz

                     + 16.0 * pa_xzzz * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_yzzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xz,
                                     double pa_xzz,
                                     double pa_xzzz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_yzzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.875 * pa_x * fx * fx * fx * pb_y

                     + 2.25 * pa_xzz * fx * fx * pb_y

                     + 6.75 * pa_xz * fx * fx * pb_yz

                     + 2.25 * pa_x * fx * fx * pb_yzz

                     + 1.5 * pa_xzzz * pb_yz * fx

                     + 4.5 * pa_xzz * fx * pb_yzz

                     + 1.5 * pa_xz * fx * pb_yzzz

                     + pa_xzzz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_yzzz_r_0(double fga,
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
                                     double pb_yzzz,
                                     double r_0_0)
    {
        return r_0_0 * (18.75 * pa_x * fx * fx * fx * fz * pb_y

                     - 2.25 * pa_x * fx * fx * pb_y * fz * fgb

                     - 2.25 * pa_x * fx * fx * fz * fga * pb_y

                     - 4.5 * pa_xzz * fx * pb_y * fz * fgb

                     + 27.0 * pa_xzz * fz * fx * fx * pb_y

                     + 81.0 * pa_xz * fx * fx * fz * pb_yz

                     - 4.5 * pa_xz * fx * pb_yz * fz * fgb

                     - 4.5 * pa_xz * fz * fga * pb_yz * fx

                     - 4.5 * pa_x * fx * fz * fga * pb_yzz

                     - 3.0 * pa_xzzz * pb_yz * fz * fgb

                     + 27.0 * pa_x * fx * fx * fz * pb_yzz

                     + 21.0 * pa_xzzz * fz * pb_yz * fx

                     + 63.0 * pa_xzz * fz * fx * pb_yzz

                     - 3.0 * pa_xz * fz * fga * pb_yzzz

                     + 21.0 * pa_xz * fz * fx * pb_yzzz

                     + 16.0 * pa_xzzz * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_zzzz_s_0(double fx,
                                     double pa_x,
                                     double pa_xz,
                                     double pa_xzz,
                                     double pa_xzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double pb_zzzz,
                                     double s_0_0)
    {
        return s_0_0 * (5.625 * pa_xz * fx * fx * fx

                     + 7.5 * pa_x * fx * fx * fx * pb_z

                     + 0.75 * pa_xzzz * fx * fx

                     + 9.0 * pa_xzz * fx * fx * pb_z

                     + 13.5 * pa_xz * fx * fx * pb_zz

                     + 3.0 * pa_x * fx * fx * pb_zzz

                     + 3.0 * pa_xzzz * pb_zz * fx

                     + 6.0 * pa_xzz * fx * pb_zzz

                     + 1.5 * pa_xz * fx * pb_zzzz

                     + pa_xzzz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_zzzz_r_0(double fga,
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
                                     double pb_zzzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 13.5 * pa_xz * fx * fx * fz * fgb

                     + 56.25 * pa_xz * fx * fx * fx * fz

                     + 75.0 * pa_x * fx * fx * fx * fz * pb_z

                     - 2.25 * pa_xz * fz * fga * fx * fx

                     - 9.0 * pa_x * fx * fx * pb_z * fz * fgb

                     - 9.0 * pa_x * fx * fx * fz * fga * pb_z

                     - 3.0 * pa_xzzz * fx * fz * fgb

                     - 18.0 * pa_xzz * fx * pb_z * fz * fgb

                     + 9.0 * pa_xzzz * fz * fx * fx

                     + 108.0 * pa_xzz * fz * fx * fx * pb_z

                     + 162.0 * pa_xz * fx * fx * fz * pb_zz

                     - 9.0 * pa_xz * fx * pb_zz * fz * fgb

                     - 9.0 * pa_xz * fz * fga * pb_zz * fx

                     - 6.0 * pa_x * fx * fz * fga * pb_zzz

                     - 6.0 * pa_xzzz * pb_zz * fz * fgb

                     + 36.0 * pa_x * fx * fx * fz * pb_zzz

                     + 42.0 * pa_xzzz * fz * pb_zz * fx

                     + 84.0 * pa_xzz * fz * fx * pb_zzz

                     - 3.0 * pa_xz * fz * fga * pb_zzzz

                     + 21.0 * pa_xz * fz * fx * pb_zzzz

                     + 16.0 * pa_xzzz * fz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xxxx_s_0(double fx,
                                     double pa_yy,
                                     double pa_yyyy,
                                     double pb_xx,
                                     double pb_xxxx,
                                     double s_0_0)
    {
        return s_0_0 * (0.5625 * fx * fx * fx * fx

                     + 2.25 * pa_yy * fx * fx * fx

                     + 0.75 * pa_yyyy * fx * fx

                     + 2.25 * fx * fx * fx * pb_xx

                     + 9.0 * pa_yy * fx * fx * pb_xx

                     + 3.0 * pa_yyyy * pb_xx * fx

                     + 0.75 * fx * fx * pb_xxxx

                     + 3.0 * pa_yy * fx * pb_xxxx

                     + pa_yyyy * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xxxx_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_yy,
                                     double pa_yyyy,
                                     double pb_xx,
                                     double pb_xxxx,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * fx * fz * fgb

                     - 2.25 * fx * fx * fx * fz * fga

                     - 9.0 * pa_yy * fx * fx * fz * fgb

                     - 4.5 * pa_yy * fz * fga * fx * fx

                     + 4.5 * fx * fx * fx * fx * fz

                     - 3.0 * pa_yyyy * fx * fz * fgb

                     + 22.5 * pa_yy * fz * fx * fx * fx

                     + 9.0 * pa_yyyy * fz * fx * fx

                     - 4.5 * fx * fx * pb_xx * fz * fgb

                     - 9.0 * fx * fx * fz * fga * pb_xx

                     - 18.0 * pa_yy * fx * pb_xx * fz * fgb

                     - 18.0 * pa_yy * fz * fga * pb_xx * fx

                     + 22.5 * fx * fx * fx * fz * pb_xx

                     - 6.0 * pa_yyyy * pb_xx * fz * fgb

                     + 108.0 * pa_yy * fz * fx * fx * pb_xx

                     + 42.0 * pa_yyyy * fz * pb_xx * fx

                     - 3.0 * fx * fz * fga * pb_xxxx

                     - 6.0 * pa_yy * fz * fga * pb_xxxx

                     + 9.0 * fx * fx * fz * pb_xxxx

                     + 42.0 * pa_yy * fz * fx * pb_xxxx

                     + 16.0 * pa_yyyy * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xxxy_s_0(double fx,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyy,
                                     double pb_x,
                                     double pb_xxx,
                                     double pb_xxxy,
                                     double pb_xy,
                                     double s_0_0)
    {
        return s_0_0 * (4.5 * pa_y * fx * fx * fx * pb_x

                     + 3.0 * pa_yyy * fx * fx * pb_x

                     + 1.125 * fx * fx * fx * pb_xy

                     + 4.5 * pa_yy * fx * fx * pb_xy

                     + 3.0 * pa_y * fx * fx * pb_xxx

                     + 1.5 * pa_yyyy * pb_xy * fx

                     + 2.0 * pa_yyy * fx * pb_xxx

                     + 0.75 * fx * fx * pb_xxxy

                     + 3.0 * pa_yy * fx * pb_xxxy

                     + pa_yyyy * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xxxy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyy,
                                     double pb_x,
                                     double pb_xxx,
                                     double pb_xxxy,
                                     double pb_xy,
                                     double r_0_0)
    {
        return r_0_0 * (- 9.0 * pa_y * fx * fx * pb_x * fz * fgb

                     - 9.0 * pa_y * fx * fx * fz * fga * pb_x

                     - 6.0 * pa_yyy * fx * pb_x * fz * fgb

                     + 45.0 * pa_y * fx * fx * fx * fz * pb_x

                     + 36.0 * pa_yyy * fz * fx * fx * pb_x

                     - 2.25 * fx * fx * pb_xy * fz * fgb

                     - 4.5 * fx * fx * fz * fga * pb_xy

                     - 9.0 * pa_yy * fx * pb_xy * fz * fgb

                     - 9.0 * pa_yy * fz * fga * pb_xy * fx

                     - 6.0 * pa_y * fx * fz * fga * pb_xxx

                     + 11.25 * fx * fx * fx * fz * pb_xy

                     - 3.0 * pa_yyyy * pb_xy * fz * fgb

                     + 54.0 * pa_yy * fz * fx * fx * pb_xy

                     + 36.0 * pa_y * fx * fx * fz * pb_xxx

                     + 21.0 * pa_yyyy * fz * pb_xy * fx

                     + 28.0 * pa_yyy * fz * fx * pb_xxx

                     - 3.0 * fx * fz * fga * pb_xxxy

                     - 6.0 * pa_yy * fz * fga * pb_xxxy

                     + 9.0 * fx * fx * fz * pb_xxxy

                     + 42.0 * pa_yy * fz * fx * pb_xxxy

                     + 16.0 * pa_yyyy * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xxxz_s_0(double fx,
                                     double pa_yy,
                                     double pa_yyyy,
                                     double pb_xxxz,
                                     double pb_xz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pb_xz

                     + 4.5 * pa_yy * fx * fx * pb_xz

                     + 1.5 * pa_yyyy * pb_xz * fx

                     + 0.75 * fx * fx * pb_xxxz

                     + 3.0 * pa_yy * fx * pb_xxxz

                     + pa_yyyy * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xxxz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_yy,
                                     double pa_yyyy,
                                     double pb_xxxz,
                                     double pb_xz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pb_xz * fz * fgb

                     - 4.5 * fx * fx * fz * fga * pb_xz

                     - 9.0 * pa_yy * fx * pb_xz * fz * fgb

                     - 9.0 * pa_yy * fz * fga * pb_xz * fx

                     + 11.25 * fx * fx * fx * fz * pb_xz

                     - 3.0 * pa_yyyy * pb_xz * fz * fgb

                     + 54.0 * pa_yy * fz * fx * fx * pb_xz

                     + 21.0 * pa_yyyy * fz * pb_xz * fx

                     - 3.0 * fx * fz * fga * pb_xxxz

                     - 6.0 * pa_yy * fz * fga * pb_xxxz

                     + 9.0 * fx * fx * fz * pb_xxxz

                     + 42.0 * pa_yy * fz * fx * pb_xxxz

                     + 16.0 * pa_yyyy * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xxyy_s_0(double fx,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyy,
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyy,
                                     double pb_y,
                                     double pb_yy,
                                     double s_0_0)
    {
        return s_0_0 * (0.9375 * fx * fx * fx * fx

                     + 2.25 * pa_yy * fx * fx * fx

                     + 3.0 * pa_y * fx * fx * fx * pb_y

                     + 1.875 * fx * fx * fx * pb_xx

                     + 0.25 * pa_yyyy * fx * fx

                     + 2.0 * pa_yyy * fx * fx * pb_y

                     + 4.5 * pa_yy * fx * fx * pb_xx

                     + 0.375 * fx * fx * fx * pb_yy

                     + 1.5 * pa_yy * fx * fx * pb_yy

                     + 6.0 * pa_y * fx * fx * pb_xxy

                     + 0.5 * pa_yyyy * pb_xx * fx

                     + 0.5 * pa_yyyy * fx * pb_yy

                     + 4.0 * pa_yyy * fx * pb_xxy

                     + 0.75 * fx * fx * pb_xxyy

                     + 3.0 * pa_yy * fx * pb_xxyy

                     + pa_yyyy * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xxyy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyy,
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyy,
                                     double pb_y,
                                     double pb_yy,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * fx * fz * fgb

                     - 2.25 * fx * fx * fx * fz * fga

                     - 6.0 * pa_yy * fx * fx * fz * fgb

                     + 7.5 * fx * fx * fx * fx * fz

                     + 22.5 * pa_yy * fx * fx * fx * fz

                     - 1.5 * pa_yy * fz * fga * fx * fx

                     - 6.0 * pa_y * fx * fx * fz * fgb * pb_y

                     - 6.0 * pa_y * fx * fx * fz * fga * pb_y

                     - 4.5 * fx * fx * fz * fga * pb_xx

                     - pa_yyyy * fx * fz * fgb

                     - 4.0 * pa_yyy * fx * fz * fgb * pb_y

                     + 30.0 * pa_y * fx * fx * fx * fz * pb_y

                     + 18.75 * fx * fx * fx * fz * pb_xx

                     + 3.0 * pa_yyyy * fz * fx * fx

                     + 24.0 * pa_yyy * fz * fx * fx * pb_y

                     + 54.0 * pa_yy * fx * fx * fz * pb_xx

                     - 0.75 * fx * fx * pb_xx * fz * fgb

                     - 0.75 * fx * fx * fz * fgb * pb_yy

                     - 1.5 * fx * fx * fz * fga * pb_yy

                     - 3.0 * pa_yy * fx * pb_xx * fz * fgb

                     - 3.0 * pa_yy * fx * fz * fgb * pb_yy

                     - 3.0 * pa_yy * fz * fga * pb_xx * fx

                     - 3.0 * pa_yy * fz * fga * fx * pb_yy

                     - 12.0 * pa_y * fx * fz * fga * pb_xxy

                     + 3.75 * fx * fx * fx * fz * pb_yy

                     - pa_yyyy * pb_xx * fz * fgb

                     - pa_yyyy * fz * fgb * pb_yy

                     + 18.0 * pa_yy * fz * fx * fx * pb_yy

                     + 72.0 * pa_y * fx * fx * fz * pb_xxy

                     + 7.0 * pa_yyyy * fz * pb_xx * fx

                     + 7.0 * pa_yyyy * fz * fx * pb_yy

                     + 56.0 * pa_yyy * fz * fx * pb_xxy

                     - 3.0 * fx * fz * fga * pb_xxyy

                     - 6.0 * pa_yy * fz * fga * pb_xxyy

                     + 9.0 * fx * fx * fz * pb_xxyy

                     + 42.0 * pa_yy * fz * fx * pb_xxyy

                     + 16.0 * pa_yyyy * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xxyz_s_0(double fx,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyy,
                                     double pb_xxyz,
                                     double pb_xxz,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (1.5 * pa_y * fx * fx * fx * pb_z

                     + pa_yyy * fx * fx * pb_z

                     + 0.375 * fx * fx * fx * pb_yz

                     + 1.5 * pa_yy * fx * fx * pb_yz

                     + 3.0 * pa_y * fx * fx * pb_xxz

                     + 0.5 * pa_yyyy * fx * pb_yz

                     + 2.0 * pa_yyy * fx * pb_xxz

                     + 0.75 * fx * fx * pb_xxyz

                     + 3.0 * pa_yy * fx * pb_xxyz

                     + pa_yyyy * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xxyz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyy,
                                     double pb_xxyz,
                                     double pb_xxz,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_y * fx * fx * fz * fgb * pb_z

                     - 3.0 * pa_y * fx * fx * fz * fga * pb_z

                     - 2.0 * pa_yyy * fx * fz * fgb * pb_z

                     + 15.0 * pa_y * fx * fx * fx * fz * pb_z

                     + 12.0 * pa_yyy * fz * fx * fx * pb_z

                     - 0.75 * fx * fx * fz * fgb * pb_yz

                     - 1.5 * fx * fx * fz * fga * pb_yz

                     - 3.0 * pa_yy * fx * fz * fgb * pb_yz

                     - 3.0 * pa_yy * fz * fga * fx * pb_yz

                     - 6.0 * pa_y * fx * fz * fga * pb_xxz

                     + 3.75 * fx * fx * fx * fz * pb_yz

                     - pa_yyyy * fz * fgb * pb_yz

                     + 18.0 * pa_yy * fz * fx * fx * pb_yz

                     + 36.0 * pa_y * fx * fx * fz * pb_xxz

                     + 7.0 * pa_yyyy * fz * fx * pb_yz

                     + 28.0 * pa_yyy * fz * fx * pb_xxz

                     - 3.0 * fx * fz * fga * pb_xxyz

                     - 6.0 * pa_yy * fz * fga * pb_xxyz

                     + 9.0 * fx * fx * fz * pb_xxyz

                     + 42.0 * pa_yy * fz * fx * pb_xxyz

                     + 16.0 * pa_yyyy * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xxzz_s_0(double fx,
                                     double pa_yy,
                                     double pa_yyyy,
                                     double pb_xx,
                                     double pb_xxzz,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.75 * pa_yy * fx * fx * fx

                     + 0.25 * pa_yyyy * fx * fx

                     + 0.375 * fx * fx * fx * pb_xx

                     + 0.375 * fx * fx * fx * pb_zz

                     + 1.5 * pa_yy * fx * fx * pb_xx

                     + 1.5 * pa_yy * fx * fx * pb_zz

                     + 0.5 * pa_yyyy * pb_xx * fx

                     + 0.5 * pa_yyyy * fx * pb_zz

                     + 0.75 * fx * fx * pb_xxzz

                     + 3.0 * pa_yy * fx * pb_xxzz

                     + pa_yyyy * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xxzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_yy,
                                     double pa_yyyy,
                                     double pb_xx,
                                     double pb_xxzz,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fx * fz * fgb

                     - 0.75 * fx * fx * fx * fz * fga

                     - 3.0 * pa_yy * fx * fx * fz * fgb

                     - 1.5 * pa_yy * fz * fga * fx * fx

                     + 1.5 * fx * fx * fx * fx * fz

                     - pa_yyyy * fx * fz * fgb

                     + 7.5 * pa_yy * fz * fx * fx * fx

                     + 3.0 * pa_yyyy * fz * fx * fx

                     - 0.75 * fx * fx * pb_xx * fz * fgb

                     - 0.75 * fx * fx * fz * fgb * pb_zz

                     - 1.5 * fx * fx * fz * fga * pb_xx

                     - 1.5 * fx * fx * fz * fga * pb_zz

                     - 3.0 * pa_yy * fx * pb_xx * fz * fgb

                     - 3.0 * pa_yy * fx * fz * fgb * pb_zz

                     - 3.0 * pa_yy * fz * fga * pb_xx * fx

                     - 3.0 * pa_yy * fz * fga * fx * pb_zz

                     + 3.75 * fx * fx * fx * fz * pb_xx

                     + 3.75 * fx * fx * fx * fz * pb_zz

                     - pa_yyyy * pb_xx * fz * fgb

                     - pa_yyyy * fz * fgb * pb_zz

                     + 18.0 * pa_yy * fz * fx * fx * pb_xx

                     + 18.0 * pa_yy * fz * fx * fx * pb_zz

                     + 7.0 * pa_yyyy * fz * pb_xx * fx

                     + 7.0 * pa_yyyy * fz * fx * pb_zz

                     - 3.0 * fx * fz * fga * pb_xxzz

                     - 6.0 * pa_yy * fz * fga * pb_xxzz

                     + 9.0 * fx * fx * fz * pb_xxzz

                     + 42.0 * pa_yy * fz * fx * pb_xxzz

                     + 16.0 * pa_yyyy * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xyyy_s_0(double fx,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyy,
                                     double pb_x,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_xyyy,
                                     double s_0_0)
    {
        return s_0_0 * (7.5 * pa_y * fx * fx * fx * pb_x

                     + 5.625 * fx * fx * fx * pb_xy

                     + 3.0 * pa_yyy * fx * fx * pb_x

                     + 13.5 * pa_yy * fx * fx * pb_xy

                     + 9.0 * pa_y * fx * fx * pb_xyy

                     + 1.5 * pa_yyyy * pb_xy * fx

                     + 6.0 * pa_yyy * fx * pb_xyy

                     + 0.75 * fx * fx * pb_xyyy

                     + 3.0 * pa_yy * fx * pb_xyyy

                     + pa_yyyy * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xyyy_r_0(double fga,
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
                                     double pb_xyyy,
                                     double r_0_0)
    {
        return r_0_0 * (75.0 * pa_y * fx * fx * fx * fz * pb_x

                     - 9.0 * pa_y * fx * fx * pb_x * fz * fgb

                     - 9.0 * pa_y * fx * fx * fz * fga * pb_x

                     - 13.5 * fx * fx * fz * fga * pb_xy

                     - 6.0 * pa_yyy * fx * pb_x * fz * fgb

                     + 56.25 * fx * fx * fx * fz * pb_xy

                     + 36.0 * pa_yyy * fz * fx * fx * pb_x

                     + 162.0 * pa_yy * fx * fx * fz * pb_xy

                     - 2.25 * fx * fx * pb_xy * fz * fgb

                     - 9.0 * pa_yy * fx * pb_xy * fz * fgb

                     - 9.0 * pa_yy * fz * fga * pb_xy * fx

                     - 18.0 * pa_y * fx * fz * fga * pb_xyy

                     - 3.0 * pa_yyyy * pb_xy * fz * fgb

                     + 108.0 * pa_y * fx * fx * fz * pb_xyy

                     + 21.0 * pa_yyyy * fz * pb_xy * fx

                     + 84.0 * pa_yyy * fz * fx * pb_xyy

                     - 3.0 * fx * fz * fga * pb_xyyy

                     - 6.0 * pa_yy * fz * fga * pb_xyyy

                     + 9.0 * fx * fx * fz * pb_xyyy

                     + 42.0 * pa_yy * fz * fx * pb_xyyy

                     + 16.0 * pa_yyyy * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xyyz_s_0(double fx,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyy,
                                     double pb_xyyz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pb_xz

                     + 4.5 * pa_yy * fx * fx * pb_xz

                     + 6.0 * pa_y * fx * fx * pb_xyz

                     + 0.5 * pa_yyyy * pb_xz * fx

                     + 4.0 * pa_yyy * fx * pb_xyz

                     + 0.75 * fx * fx * pb_xyyz

                     + 3.0 * pa_yy * fx * pb_xyyz

                     + pa_yyyy * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xyyz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyy,
                                     double pb_xyyz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double r_0_0)
    {
        return r_0_0 * (- 4.5 * fx * fx * fz * fga * pb_xz

                     + 18.75 * fx * fx * fx * fz * pb_xz

                     + 54.0 * pa_yy * fx * fx * fz * pb_xz

                     - 0.75 * fx * fx * pb_xz * fz * fgb

                     - 3.0 * pa_yy * fx * pb_xz * fz * fgb

                     - 3.0 * pa_yy * fz * fga * pb_xz * fx

                     - 12.0 * pa_y * fx * fz * fga * pb_xyz

                     - pa_yyyy * pb_xz * fz * fgb

                     + 72.0 * pa_y * fx * fx * fz * pb_xyz

                     + 7.0 * pa_yyyy * fz * pb_xz * fx

                     + 56.0 * pa_yyy * fz * fx * pb_xyz

                     - 3.0 * fx * fz * fga * pb_xyyz

                     - 6.0 * pa_yy * fz * fga * pb_xyyz

                     + 9.0 * fx * fx * fz * pb_xyyz

                     + 42.0 * pa_yy * fz * fx * pb_xyyz

                     + 16.0 * pa_yyyy * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xyzz_s_0(double fx,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyy,
                                     double pb_x,
                                     double pb_xy,
                                     double pb_xyzz,
                                     double pb_xzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.5 * pa_y * fx * fx * fx * pb_x

                     + pa_yyy * fx * fx * pb_x

                     + 0.375 * fx * fx * fx * pb_xy

                     + 1.5 * pa_yy * fx * fx * pb_xy

                     + 3.0 * pa_y * fx * fx * pb_xzz

                     + 0.5 * pa_yyyy * pb_xy * fx

                     + 2.0 * pa_yyy * fx * pb_xzz

                     + 0.75 * fx * fx * pb_xyzz

                     + 3.0 * pa_yy * fx * pb_xyzz

                     + pa_yyyy * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xyzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyy,
                                     double pb_x,
                                     double pb_xy,
                                     double pb_xyzz,
                                     double pb_xzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_y * fx * fx * pb_x * fz * fgb

                     - 3.0 * pa_y * fx * fx * fz * fga * pb_x

                     - 2.0 * pa_yyy * fx * pb_x * fz * fgb

                     + 15.0 * pa_y * fx * fx * fx * fz * pb_x

                     + 12.0 * pa_yyy * fz * fx * fx * pb_x

                     - 0.75 * fx * fx * pb_xy * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_xy

                     - 3.0 * pa_yy * fx * pb_xy * fz * fgb

                     - 3.0 * pa_yy * fz * fga * pb_xy * fx

                     - 6.0 * pa_y * fx * fz * fga * pb_xzz

                     + 3.75 * fx * fx * fx * fz * pb_xy

                     - pa_yyyy * pb_xy * fz * fgb

                     + 18.0 * pa_yy * fz * fx * fx * pb_xy

                     + 36.0 * pa_y * fx * fx * fz * pb_xzz

                     + 7.0 * pa_yyyy * fz * pb_xy * fx

                     + 28.0 * pa_yyy * fz * fx * pb_xzz

                     - 3.0 * fx * fz * fga * pb_xyzz

                     - 6.0 * pa_yy * fz * fga * pb_xyzz

                     + 9.0 * fx * fx * fz * pb_xyzz

                     + 42.0 * pa_yy * fz * fx * pb_xyzz

                     + 16.0 * pa_yyyy * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xzzz_s_0(double fx,
                                     double pa_yy,
                                     double pa_yyyy,
                                     double pb_xz,
                                     double pb_xzzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pb_xz

                     + 4.5 * pa_yy * fx * fx * pb_xz

                     + 1.5 * pa_yyyy * pb_xz * fx

                     + 0.75 * fx * fx * pb_xzzz

                     + 3.0 * pa_yy * fx * pb_xzzz

                     + pa_yyyy * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xzzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_yy,
                                     double pa_yyyy,
                                     double pb_xz,
                                     double pb_xzzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pb_xz * fz * fgb

                     - 4.5 * fx * fx * fz * fga * pb_xz

                     - 9.0 * pa_yy * fx * pb_xz * fz * fgb

                     - 9.0 * pa_yy * fz * fga * pb_xz * fx

                     + 11.25 * fx * fx * fx * fz * pb_xz

                     - 3.0 * pa_yyyy * pb_xz * fz * fgb

                     + 54.0 * pa_yy * fz * fx * fx * pb_xz

                     + 21.0 * pa_yyyy * fz * pb_xz * fx

                     - 3.0 * fx * fz * fga * pb_xzzz

                     - 6.0 * pa_yy * fz * fga * pb_xzzz

                     + 9.0 * fx * fx * fz * pb_xzzz

                     + 42.0 * pa_yy * fz * fx * pb_xzzz

                     + 16.0 * pa_yyyy * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_yyyy_s_0(double fx,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double pb_yyyy,
                                     double s_0_0)
    {
        return s_0_0 * (6.5625 * fx * fx * fx * fx

                     + 11.25 * pa_yy * fx * fx * fx

                     + 30.0 * pa_y * fx * fx * fx * pb_y

                     + 11.25 * fx * fx * fx * pb_yy

                     + 0.75 * pa_yyyy * fx * fx

                     + 12.0 * pa_yyy * fx * fx * pb_y

                     + 27.0 * pa_yy * fx * fx * pb_yy

                     + 12.0 * pa_y * fx * fx * pb_yyy

                     + 3.0 * pa_yyyy * pb_yy * fx

                     + 8.0 * pa_yyy * fx * pb_yyy

                     + 0.75 * fx * fx * pb_yyyy

                     + 3.0 * pa_yy * fx * pb_yyyy

                     + pa_yyyy * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_yyyy_r_0(double fga,
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
                                     double pb_yyyy,
                                     double r_0_0)
    {
        return r_0_0 * (52.5 * fx * fx * fx * fx * fz

                     - 11.25 * fx * fx * fx * fz * fgb

                     - 11.25 * fx * fx * fx * fz * fga

                     - 27.0 * pa_yy * fx * fx * fz * fgb

                     + 112.5 * pa_yy * fx * fx * fx * fz

                     + 300.0 * pa_y * fx * fx * fx * fz * pb_y

                     - 4.5 * pa_yy * fz * fga * fx * fx

                     - 36.0 * pa_y * fx * fx * pb_y * fz * fgb

                     - 36.0 * pa_y * fx * fx * fz * fga * pb_y

                     - 27.0 * fx * fx * fz * fga * pb_yy

                     - 3.0 * pa_yyyy * fx * fz * fgb

                     - 24.0 * pa_yyy * fx * pb_y * fz * fgb

                     + 112.5 * fx * fx * fx * fz * pb_yy

                     + 9.0 * pa_yyyy * fz * fx * fx

                     + 144.0 * pa_yyy * fz * fx * fx * pb_y

                     + 324.0 * pa_yy * fx * fx * fz * pb_yy

                     - 4.5 * fx * fx * pb_yy * fz * fgb

                     - 18.0 * pa_yy * fx * pb_yy * fz * fgb

                     - 18.0 * pa_yy * fz * fga * pb_yy * fx

                     - 24.0 * pa_y * fx * fz * fga * pb_yyy

                     - 6.0 * pa_yyyy * pb_yy * fz * fgb

                     + 144.0 * pa_y * fx * fx * fz * pb_yyy

                     + 42.0 * pa_yyyy * fz * pb_yy * fx

                     + 112.0 * pa_yyy * fz * fx * pb_yyy

                     - 3.0 * fx * fz * fga * pb_yyyy

                     - 6.0 * pa_yy * fz * fga * pb_yyyy

                     + 9.0 * fx * fx * fz * pb_yyyy

                     + 42.0 * pa_yy * fz * fx * pb_yyyy

                     + 16.0 * pa_yyyy * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_yyyz_s_0(double fx,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyy,
                                     double pb_yyyz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (7.5 * pa_y * fx * fx * fx * pb_z

                     + 5.625 * fx * fx * fx * pb_yz

                     + 3.0 * pa_yyy * fx * fx * pb_z

                     + 13.5 * pa_yy * fx * fx * pb_yz

                     + 9.0 * pa_y * fx * fx * pb_yyz

                     + 1.5 * pa_yyyy * pb_yz * fx

                     + 6.0 * pa_yyy * fx * pb_yyz

                     + 0.75 * fx * fx * pb_yyyz

                     + 3.0 * pa_yy * fx * pb_yyyz

                     + pa_yyyy * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_yyyz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyy,
                                     double pb_yyyz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (75.0 * pa_y * fx * fx * fx * fz * pb_z

                     - 9.0 * pa_y * fx * fx * fz * fgb * pb_z

                     - 9.0 * pa_y * fx * fx * fz * fga * pb_z

                     - 13.5 * fx * fx * fz * fga * pb_yz

                     - 6.0 * pa_yyy * fx * fz * fgb * pb_z

                     + 56.25 * fx * fx * fx * fz * pb_yz

                     + 36.0 * pa_yyy * fz * fx * fx * pb_z

                     + 162.0 * pa_yy * fx * fx * fz * pb_yz

                     - 2.25 * fx * fx * pb_yz * fz * fgb

                     - 9.0 * pa_yy * fx * pb_yz * fz * fgb

                     - 9.0 * pa_yy * fz * fga * pb_yz * fx

                     - 18.0 * pa_y * fx * fz * fga * pb_yyz

                     - 3.0 * pa_yyyy * pb_yz * fz * fgb

                     + 108.0 * pa_y * fx * fx * fz * pb_yyz

                     + 21.0 * pa_yyyy * fz * pb_yz * fx

                     + 84.0 * pa_yyy * fz * fx * pb_yyz

                     - 3.0 * fx * fz * fga * pb_yyyz

                     - 6.0 * pa_yy * fz * fga * pb_yyyz

                     + 9.0 * fx * fx * fz * pb_yyyz

                     + 42.0 * pa_yy * fz * fx * pb_yyyz

                     + 16.0 * pa_yyyy * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_yyzz_s_0(double fx,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyzz,
                                     double pb_yzz,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.9375 * fx * fx * fx * fx

                     + 2.25 * pa_yy * fx * fx * fx

                     + 3.0 * pa_y * fx * fx * fx * pb_y

                     + 1.875 * fx * fx * fx * pb_zz

                     + 0.25 * pa_yyyy * fx * fx

                     + 2.0 * pa_yyy * fx * fx * pb_y

                     + 4.5 * pa_yy * fx * fx * pb_zz

                     + 0.375 * fx * fx * fx * pb_yy

                     + 1.5 * pa_yy * fx * fx * pb_yy

                     + 6.0 * pa_y * fx * fx * pb_yzz

                     + 0.5 * pa_yyyy * pb_yy * fx

                     + 0.5 * pa_yyyy * fx * pb_zz

                     + 4.0 * pa_yyy * fx * pb_yzz

                     + 0.75 * fx * fx * pb_yyzz

                     + 3.0 * pa_yy * fx * pb_yyzz

                     + pa_yyyy * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_yyzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyy,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyzz,
                                     double pb_yzz,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * fx * fz * fgb

                     - 2.25 * fx * fx * fx * fz * fga

                     - 6.0 * pa_yy * fx * fx * fz * fgb

                     + 7.5 * fx * fx * fx * fx * fz

                     + 22.5 * pa_yy * fx * fx * fx * fz

                     - 1.5 * pa_yy * fz * fga * fx * fx

                     - 6.0 * pa_y * fx * fx * pb_y * fz * fgb

                     - 6.0 * pa_y * fx * fx * fz * fga * pb_y

                     - 4.5 * fx * fx * fz * fga * pb_zz

                     - pa_yyyy * fx * fz * fgb

                     - 4.0 * pa_yyy * fx * pb_y * fz * fgb

                     + 30.0 * pa_y * fx * fx * fx * fz * pb_y

                     + 18.75 * fx * fx * fx * fz * pb_zz

                     + 3.0 * pa_yyyy * fz * fx * fx

                     + 24.0 * pa_yyy * fz * fx * fx * pb_y

                     + 54.0 * pa_yy * fx * fx * fz * pb_zz

                     - 0.75 * fx * fx * pb_yy * fz * fgb

                     - 0.75 * fx * fx * fz * fgb * pb_zz

                     - 1.5 * fx * fx * fz * fga * pb_yy

                     - 3.0 * pa_yy * fx * pb_yy * fz * fgb

                     - 3.0 * pa_yy * fx * fz * fgb * pb_zz

                     - 3.0 * pa_yy * fz * fga * pb_yy * fx

                     - 3.0 * pa_yy * fz * fga * fx * pb_zz

                     - 12.0 * pa_y * fx * fz * fga * pb_yzz

                     + 3.75 * fx * fx * fx * fz * pb_yy

                     - pa_yyyy * pb_yy * fz * fgb

                     - pa_yyyy * fz * fgb * pb_zz

                     + 18.0 * pa_yy * fz * fx * fx * pb_yy

                     + 72.0 * pa_y * fx * fx * fz * pb_yzz

                     + 7.0 * pa_yyyy * fz * pb_yy * fx

                     + 7.0 * pa_yyyy * fz * fx * pb_zz

                     + 56.0 * pa_yyy * fz * fx * pb_yzz

                     - 3.0 * fx * fz * fga * pb_yyzz

                     - 6.0 * pa_yy * fz * fga * pb_yyzz

                     + 9.0 * fx * fx * fz * pb_yyzz

                     + 42.0 * pa_yy * fz * fx * pb_yyzz

                     + 16.0 * pa_yyyy * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_yzzz_s_0(double fx,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyy,
                                     double pb_yz,
                                     double pb_yzzz,
                                     double pb_z,
                                     double pb_zzz,
                                     double s_0_0)
    {
        return s_0_0 * (4.5 * pa_y * fx * fx * fx * pb_z

                     + 3.0 * pa_yyy * fx * fx * pb_z

                     + 1.125 * fx * fx * fx * pb_yz

                     + 4.5 * pa_yy * fx * fx * pb_yz

                     + 3.0 * pa_y * fx * fx * pb_zzz

                     + 1.5 * pa_yyyy * pb_yz * fx

                     + 2.0 * pa_yyy * fx * pb_zzz

                     + 0.75 * fx * fx * pb_yzzz

                     + 3.0 * pa_yy * fx * pb_yzzz

                     + pa_yyyy * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_yzzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyy,
                                     double pb_yz,
                                     double pb_yzzz,
                                     double pb_z,
                                     double pb_zzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 9.0 * pa_y * fx * fx * pb_z * fz * fgb

                     - 9.0 * pa_y * fx * fx * fz * fga * pb_z

                     - 6.0 * pa_yyy * fx * pb_z * fz * fgb

                     + 45.0 * pa_y * fx * fx * fx * fz * pb_z

                     + 36.0 * pa_yyy * fz * fx * fx * pb_z

                     - 2.25 * fx * fx * pb_yz * fz * fgb

                     - 4.5 * fx * fx * fz * fga * pb_yz

                     - 9.0 * pa_yy * fx * pb_yz * fz * fgb

                     - 9.0 * pa_yy * fz * fga * pb_yz * fx

                     - 6.0 * pa_y * fx * fz * fga * pb_zzz

                     + 11.25 * fx * fx * fx * fz * pb_yz

                     - 3.0 * pa_yyyy * pb_yz * fz * fgb

                     + 54.0 * pa_yy * fz * fx * fx * pb_yz

                     + 36.0 * pa_y * fx * fx * fz * pb_zzz

                     + 21.0 * pa_yyyy * fz * pb_yz * fx

                     + 28.0 * pa_yyy * fz * fx * pb_zzz

                     - 3.0 * fx * fz * fga * pb_yzzz

                     - 6.0 * pa_yy * fz * fga * pb_yzzz

                     + 9.0 * fx * fx * fz * pb_yzzz

                     + 42.0 * pa_yy * fz * fx * pb_yzzz

                     + 16.0 * pa_yyyy * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_zzzz_s_0(double fx,
                                     double pa_yy,
                                     double pa_yyyy,
                                     double pb_zz,
                                     double pb_zzzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.5625 * fx * fx * fx * fx

                     + 2.25 * pa_yy * fx * fx * fx

                     + 0.75 * pa_yyyy * fx * fx

                     + 2.25 * fx * fx * fx * pb_zz

                     + 9.0 * pa_yy * fx * fx * pb_zz

                     + 3.0 * pa_yyyy * pb_zz * fx

                     + 0.75 * fx * fx * pb_zzzz

                     + 3.0 * pa_yy * fx * pb_zzzz

                     + pa_yyyy * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_zzzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_yy,
                                     double pa_yyyy,
                                     double pb_zz,
                                     double pb_zzzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * fx * fz * fgb

                     - 2.25 * fx * fx * fx * fz * fga

                     - 9.0 * pa_yy * fx * fx * fz * fgb

                     - 4.5 * pa_yy * fz * fga * fx * fx

                     + 4.5 * fx * fx * fx * fx * fz

                     - 3.0 * pa_yyyy * fx * fz * fgb

                     + 22.5 * pa_yy * fz * fx * fx * fx

                     + 9.0 * pa_yyyy * fz * fx * fx

                     - 4.5 * fx * fx * pb_zz * fz * fgb

                     - 9.0 * fx * fx * fz * fga * pb_zz

                     - 18.0 * pa_yy * fx * pb_zz * fz * fgb

                     - 18.0 * pa_yy * fz * fga * pb_zz * fx

                     + 22.5 * fx * fx * fx * fz * pb_zz

                     - 6.0 * pa_yyyy * pb_zz * fz * fgb

                     + 108.0 * pa_yy * fz * fx * fx * pb_zz

                     + 42.0 * pa_yyyy * fz * pb_zz * fx

                     - 3.0 * fx * fz * fga * pb_zzzz

                     - 6.0 * pa_yy * fz * fga * pb_zzzz

                     + 9.0 * fx * fx * fz * pb_zzzz

                     + 42.0 * pa_yy * fz * fx * pb_zzzz

                     + 16.0 * pa_yyyy * fz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xxxx_s_0(double fx,
                                     double pa_yyyz,
                                     double pa_yz,
                                     double pb_xx,
                                     double pb_xxxx,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_yz * fx * fx * fx

                     + 0.75 * pa_yyyz * fx * fx

                     + 4.5 * pa_yz * fx * fx * pb_xx

                     + 3.0 * pa_yyyz * pb_xx * fx

                     + 1.5 * pa_yz * fx * pb_xxxx

                     + pa_yyyz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xxxx_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_yyyz,
                                     double pa_yz,
                                     double pb_xx,
                                     double pb_xxxx,
                                     double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_yz * fx * fx * fz * fgb

                     - 2.25 * pa_yz * fz * fga * fx * fx

                     - 3.0 * pa_yyyz * fx * fz * fgb

                     + 11.25 * pa_yz * fx * fx * fx * fz

                     + 9.0 * pa_yyyz * fz * fx * fx

                     - 9.0 * pa_yz * fx * pb_xx * fz * fgb

                     - 9.0 * pa_yz * fz * fga * pb_xx * fx

                     - 6.0 * pa_yyyz * pb_xx * fz * fgb

                     + 54.0 * pa_yz * fx * fx * fz * pb_xx

                     + 42.0 * pa_yyyz * fz * pb_xx * fx

                     - 3.0 * pa_yz * fz * fga * pb_xxxx

                     + 21.0 * pa_yz * fx * fz * pb_xxxx

                     + 16.0 * pa_yyyz * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xxxy_s_0(double fx,
                                     double pa_yyyz,
                                     double pa_yyz,
                                     double pa_yz,
                                     double pa_z,
                                     double pb_x,
                                     double pb_xxx,
                                     double pb_xxxy,
                                     double pb_xy,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pa_z * pb_x

                     + 2.25 * pa_yyz * fx * fx * pb_x

                     + 2.25 * pa_yz * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_z * pb_xxx

                     + 1.5 * pa_yyyz * pb_xy * fx

                     + 1.5 * pa_yyz * fx * pb_xxx

                     + 1.5 * pa_yz * fx * pb_xxxy

                     + pa_yyyz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xxxy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_yyyz,
                                     double pa_yyz,
                                     double pa_yz,
                                     double pa_z,
                                     double pb_x,
                                     double pb_xxx,
                                     double pb_xxxy,
                                     double pb_xy,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pa_z * pb_x * fz * fgb

                     - 2.25 * fx * fx * fz * fga * pa_z * pb_x

                     - 4.5 * pa_yyz * fx * pb_x * fz * fgb

                     + 11.25 * fx * fx * fx * fz * pa_z * pb_x

                     + 27.0 * pa_yyz * fx * fx * fz * pb_x

                     - 4.5 * pa_yz * fx * pb_xy * fz * fgb

                     - 4.5 * pa_yz * fz * fga * pb_xy * fx

                     - 1.5 * fx * fz * fga * pa_z * pb_xxx

                     - 3.0 * pa_yyyz * pb_xy * fz * fgb

                     + 27.0 * pa_yz * fx * fx * fz * pb_xy

                     + 9.0 * fx * fx * fz * pa_z * pb_xxx

                     + 21.0 * pa_yyyz * fz * pb_xy * fx

                     + 21.0 * pa_yyz * fx * fz * pb_xxx

                     - 3.0 * pa_yz * fz * fga * pb_xxxy

                     + 21.0 * pa_yz * fx * fz * pb_xxxy

                     + 16.0 * pa_yyyz * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xxxz_s_0(double fx,
                                     double pa_y,
                                     double pa_yyy,
                                     double pa_yyyz,
                                     double pa_yz,
                                     double pb_x,
                                     double pb_xxx,
                                     double pb_xxxz,
                                     double pb_xz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_y * fx * fx * fx * pb_x

                     + 0.75 * pa_yyy * fx * fx * pb_x

                     + 2.25 * pa_yz * fx * fx * pb_xz

                     + 0.75 * pa_y * fx * fx * pb_xxx

                     + 1.5 * pa_yyyz * pb_xz * fx

                     + 0.5 * pa_yyy * fx * pb_xxx

                     + 1.5 * pa_yz * fx * pb_xxxz

                     + pa_yyyz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xxxz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_y,
                                     double pa_yyy,
                                     double pa_yyyz,
                                     double pa_yz,
                                     double pb_x,
                                     double pb_xxx,
                                     double pb_xxxz,
                                     double pb_xz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * pa_y * fx * fx * pb_x * fz * fgb

                     - 2.25 * pa_y * fz * fga * fx * fx * pb_x

                     - 1.5 * pa_yyy * fx * pb_x * fz * fgb

                     + 11.25 * pa_y * fx * fx * fx * fz * pb_x

                     + 9.0 * pa_yyy * fz * fx * fx * pb_x

                     - 4.5 * pa_yz * fx * pb_xz * fz * fgb

                     - 4.5 * pa_yz * fz * fga * pb_xz * fx

                     - 1.5 * pa_y * fz * fga * fx * pb_xxx

                     - 3.0 * pa_yyyz * pb_xz * fz * fgb

                     + 27.0 * pa_yz * fx * fx * fz * pb_xz

                     + 9.0 * pa_y * fx * fx * fz * pb_xxx

                     + 21.0 * pa_yyyz * fz * pb_xz * fx

                     + 7.0 * pa_yyy * fz * fx * pb_xxx

                     - 3.0 * pa_yz * fz * fga * pb_xxxz

                     + 21.0 * pa_yz * fx * fz * pb_xxxz

                     + 16.0 * pa_yyyz * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xxyy_s_0(double fx,
                                     double pa_yyyz,
                                     double pa_yyz,
                                     double pa_yz,
                                     double pa_z,
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyy,
                                     double pb_y,
                                     double pb_yy,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_yz * fx * fx * fx

                     + 0.75 * fx * fx * fx * pa_z * pb_y

                     + 0.25 * pa_yyyz * fx * fx

                     + 1.5 * pa_yyz * fx * fx * pb_y

                     + 2.25 * pa_yz * fx * fx * pb_xx

                     + 0.75 * pa_yz * fx * fx * pb_yy

                     + 1.5 * fx * fx * pa_z * pb_xxy

                     + 0.5 * pa_yyyz * pb_xx * fx

                     + 0.5 * pa_yyyz * fx * pb_yy

                     + 3.0 * pa_yyz * fx * pb_xxy

                     + 1.5 * pa_yz * fx * pb_xxyy

                     + pa_yyyz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xxyy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_yyyz,
                                     double pa_yyz,
                                     double pa_yz,
                                     double pa_z,
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyy,
                                     double pb_y,
                                     double pb_yy,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_yz * fx * fx * fz * fgb

                     + 11.25 * pa_yz * fx * fx * fx * fz

                     - 0.75 * pa_yz * fz * fga * fx * fx

                     - 1.5 * fx * fx * pa_z * fz * fgb * pb_y

                     - 1.5 * fx * fx * fz * fga * pa_z * pb_y

                     - pa_yyyz * fx * fz * fgb

                     - 3.0 * pa_yyz * fx * fz * fgb * pb_y

                     + 7.5 * fx * fx * fx * fz * pa_z * pb_y

                     + 3.0 * pa_yyyz * fz * fx * fx

                     + 18.0 * pa_yyz * fx * fx * fz * pb_y

                     + 27.0 * pa_yz * fx * fx * fz * pb_xx

                     - 1.5 * pa_yz * fx * pb_xx * fz * fgb

                     - 1.5 * pa_yz * fx * fz * fgb * pb_yy

                     - 1.5 * pa_yz * fz * fga * pb_xx * fx

                     - 1.5 * pa_yz * fz * fga * fx * pb_yy

                     - 3.0 * fx * fz * fga * pa_z * pb_xxy

                     - pa_yyyz * pb_xx * fz * fgb

                     - pa_yyyz * fz * fgb * pb_yy

                     + 9.0 * pa_yz * fx * fx * fz * pb_yy

                     + 18.0 * fx * fx * fz * pa_z * pb_xxy

                     + 7.0 * pa_yyyz * fz * pb_xx * fx

                     + 7.0 * pa_yyyz * fz * fx * pb_yy

                     + 42.0 * pa_yyz * fx * fz * pb_xxy

                     - 3.0 * pa_yz * fz * fga * pb_xxyy

                     + 21.0 * pa_yz * fx * fz * pb_xxyy

                     + 16.0 * pa_yyyz * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xxyz_s_0(double fx,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyz,
                                     double pa_yyz,
                                     double pa_yz,
                                     double pa_z,
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyz,
                                     double pb_xxz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * pa_yy * fx * fx * fx

                     + 0.375 * pa_y * fx * fx * fx * pb_y

                     + 0.375 * fx * fx * fx * pa_z * pb_z

                     + 0.375 * fx * fx * fx * pb_xx

                     + 0.25 * pa_yyy * fx * fx * pb_y

                     + 0.75 * pa_yyz * fx * fx * pb_z

                     + 0.75 * pa_yy * fx * fx * pb_xx

                     + 0.75 * pa_yz * fx * fx * pb_yz

                     + 0.75 * pa_y * fx * fx * pb_xxy

                     + 0.75 * fx * fx * pa_z * pb_xxz

                     + 0.5 * pa_yyyz * fx * pb_yz

                     + 0.5 * pa_yyy * fx * pb_xxy

                     + 1.5 * pa_yyz * fx * pb_xxz

                     + 1.5 * pa_yz * fx * pb_xxyz

                     + pa_yyyz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xxyz_r_0(double fga,
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
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyz,
                                     double pb_xxz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.375 * fx * fx * fx * fz * fgb

                     - 0.375 * fx * fx * fx * fz * fga

                     - 0.75 * pa_yy * fx * fx * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     + 3.75 * pa_yy * fx * fx * fx * fz

                     - 0.75 * pa_y * fx * fx * fz * fgb * pb_y

                     - 0.75 * pa_y * fz * fga * fx * fx * pb_y

                     - 0.75 * fx * fx * pa_z * fz * fgb * pb_z

                     - 0.75 * fx * fx * fz * fga * pa_z * pb_z

                     - 0.75 * fx * fx * fz * fga * pb_xx

                     - 0.5 * pa_yyy * fx * fz * fgb * pb_y

                     - 1.5 * pa_yyz * fx * fz * fgb * pb_z

                     + 3.75 * pa_y * fx * fx * fx * fz * pb_y

                     + 3.75 * fx * fx * fx * fz * pa_z * pb_z

                     + 3.75 * fx * fx * fx * fz * pb_xx

                     + 3.0 * pa_yyy * fz * fx * fx * pb_y

                     + 9.0 * pa_yyz * fx * fx * fz * pb_z

                     + 9.0 * pa_yy * fx * fx * fz * pb_xx

                     - 1.5 * pa_yz * fx * fz * fgb * pb_yz

                     - 1.5 * pa_yz * fz * fga * fx * pb_yz

                     - 1.5 * pa_y * fz * fga * fx * pb_xxy

                     - 1.5 * fx * fz * fga * pa_z * pb_xxz

                     - pa_yyyz * fz * fgb * pb_yz

                     + 9.0 * pa_yz * fx * fx * fz * pb_yz

                     + 9.0 * pa_y * fx * fx * fz * pb_xxy

                     + 9.0 * fx * fx * fz * pa_z * pb_xxz

                     + 7.0 * pa_yyyz * fz * fx * pb_yz

                     + 7.0 * pa_yyy * fz * fx * pb_xxy

                     + 21.0 * pa_yyz * fx * fz * pb_xxz

                     - 3.0 * pa_yz * fz * fga * pb_xxyz

                     + 21.0 * pa_yz * fx * fz * pb_xxyz

                     + 16.0 * pa_yyyz * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xxzz_s_0(double fx,
                                     double pa_y,
                                     double pa_yyy,
                                     double pa_yyyz,
                                     double pa_yz,
                                     double pb_xx,
                                     double pb_xxz,
                                     double pb_xxzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * pa_yz * fx * fx * fx

                     + 0.75 * pa_y * fx * fx * fx * pb_z

                     + 0.25 * pa_yyyz * fx * fx

                     + 0.5 * pa_yyy * fx * fx * pb_z

                     + 0.75 * pa_yz * fx * fx * pb_xx

                     + 0.75 * pa_yz * fx * fx * pb_zz

                     + 1.5 * pa_y * fx * fx * pb_xxz

                     + 0.5 * pa_yyyz * pb_xx * fx

                     + 0.5 * pa_yyyz * fx * pb_zz

                     + pa_yyy * fx * pb_xxz

                     + 1.5 * pa_yz * fx * pb_xxzz

                     + pa_yyyz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xxzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_y,
                                     double pa_yyy,
                                     double pa_yyyz,
                                     double pa_yz,
                                     double pb_xx,
                                     double pb_xxz,
                                     double pb_xxzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_yz * fx * fx * fz * fgb

                     - 1.5 * pa_y * fx * fx * fz * fgb * pb_z

                     - 0.75 * pa_yz * fz * fga * fx * fx

                     - 1.5 * pa_y * fz * fga * fx * fx * pb_z

                     - pa_yyyz * fx * fz * fgb

                     - pa_yyy * fx * fz * fgb * pb_z

                     + 3.75 * pa_yz * fx * fx * fx * fz

                     + 7.5 * pa_y * fx * fx * fx * fz * pb_z

                     + 3.0 * pa_yyyz * fz * fx * fx

                     + 6.0 * pa_yyy * fz * fx * fx * pb_z

                     - 1.5 * pa_yz * fx * pb_xx * fz * fgb

                     - 1.5 * pa_yz * fx * fz * fgb * pb_zz

                     - 1.5 * pa_yz * fz * fga * pb_xx * fx

                     - 1.5 * pa_yz * fz * fga * fx * pb_zz

                     - 3.0 * pa_y * fz * fga * fx * pb_xxz

                     - pa_yyyz * pb_xx * fz * fgb

                     - pa_yyyz * fz * fgb * pb_zz

                     + 9.0 * pa_yz * fx * fx * fz * pb_xx

                     + 9.0 * pa_yz * fx * fx * fz * pb_zz

                     + 18.0 * pa_y * fx * fx * fz * pb_xxz

                     + 7.0 * pa_yyyz * fz * pb_xx * fx

                     + 7.0 * pa_yyyz * fz * fx * pb_zz

                     + 14.0 * pa_yyy * fz * fx * pb_xxz

                     - 3.0 * pa_yz * fz * fga * pb_xxzz

                     + 21.0 * pa_yz * fx * fz * pb_xxzz

                     + 16.0 * pa_yyyz * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xyyy_s_0(double fx,
                                     double pa_yyyz,
                                     double pa_yyz,
                                     double pa_yz,
                                     double pa_z,
                                     double pb_x,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_xyyy,
                                     double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pa_z * pb_x

                     + 2.25 * pa_yyz * fx * fx * pb_x

                     + 6.75 * pa_yz * fx * fx * pb_xy

                     + 2.25 * fx * fx * pa_z * pb_xyy

                     + 1.5 * pa_yyyz * pb_xy * fx

                     + 4.5 * pa_yyz * fx * pb_xyy

                     + 1.5 * pa_yz * fx * pb_xyyy

                     + pa_yyyz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xyyy_r_0(double fga,
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
                                     double pb_xyyy,
                                     double r_0_0)
    {
        return r_0_0 * (18.75 * fx * fx * fx * fz * pa_z * pb_x

                     - 2.25 * fx * fx * pa_z * pb_x * fz * fgb

                     - 2.25 * fx * fx * fz * fga * pa_z * pb_x

                     - 4.5 * pa_yyz * fx * pb_x * fz * fgb

                     + 27.0 * pa_yyz * fx * fx * fz * pb_x

                     + 81.0 * pa_yz * fx * fx * fz * pb_xy

                     - 4.5 * pa_yz * fx * pb_xy * fz * fgb

                     - 4.5 * pa_yz * fz * fga * pb_xy * fx

                     - 4.5 * fx * fz * fga * pa_z * pb_xyy

                     - 3.0 * pa_yyyz * pb_xy * fz * fgb

                     + 27.0 * fx * fx * fz * pa_z * pb_xyy

                     + 21.0 * pa_yyyz * fz * pb_xy * fx

                     + 63.0 * pa_yyz * fx * fz * pb_xyy

                     - 3.0 * pa_yz * fz * fga * pb_xyyy

                     + 21.0 * pa_yz * fx * fz * pb_xyyy

                     + 16.0 * pa_yyyz * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xyyz_s_0(double fx,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyz,
                                     double pa_yyz,
                                     double pa_yz,
                                     double pa_z,
                                     double pb_x,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_xyyz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_y * fx * fx * fx * pb_x

                     + 0.75 * fx * fx * fx * pb_xy

                     + 0.25 * pa_yyy * fx * fx * pb_x

                     + 1.5 * pa_yy * fx * fx * pb_xy

                     + 2.25 * pa_yz * fx * fx * pb_xz

                     + 0.75 * pa_y * fx * fx * pb_xyy

                     + 1.5 * fx * fx * pa_z * pb_xyz

                     + 0.5 * pa_yyyz * pb_xz * fx

                     + 0.5 * pa_yyy * fx * pb_xyy

                     + 3.0 * pa_yyz * fx * pb_xyz

                     + 1.5 * pa_yz * fx * pb_xyyz

                     + pa_yyyz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xyyz_r_0(double fga,
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
                                     double pb_x,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_xyyz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double r_0_0)
    {
        return r_0_0 * (11.25 * pa_y * fx * fx * fx * fz * pb_x

                     - 0.75 * pa_y * fx * fx * pb_x * fz * fgb

                     - 0.75 * pa_y * fz * fga * fx * fx * pb_x

                     - 1.5 * fx * fx * fz * fga * pb_xy

                     - 0.5 * pa_yyy * fx * pb_x * fz * fgb

                     + 7.5 * fx * fx * fx * fz * pb_xy

                     + 3.0 * pa_yyy * fz * fx * fx * pb_x

                     + 18.0 * pa_yy * fx * fx * fz * pb_xy

                     + 27.0 * pa_yz * fx * fx * fz * pb_xz

                     - 1.5 * pa_yz * fx * pb_xz * fz * fgb

                     - 1.5 * pa_yz * fz * fga * pb_xz * fx

                     - 1.5 * pa_y * fz * fga * fx * pb_xyy

                     - 3.0 * fx * fz * fga * pa_z * pb_xyz

                     - pa_yyyz * pb_xz * fz * fgb

                     + 9.0 * pa_y * fx * fx * fz * pb_xyy

                     + 18.0 * fx * fx * fz * pa_z * pb_xyz

                     + 7.0 * pa_yyyz * fz * pb_xz * fx

                     + 7.0 * pa_yyy * fz * fx * pb_xyy

                     + 42.0 * pa_yyz * fx * fz * pb_xyz

                     - 3.0 * pa_yz * fz * fga * pb_xyyz

                     + 21.0 * pa_yz * fx * fz * pb_xyyz

                     + 16.0 * pa_yyyz * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xyzz_s_0(double fx,
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
                                     double pb_xyzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_z * pb_x

                     + 0.75 * fx * fx * fx * pb_xz

                     + 0.75 * pa_yyz * fx * fx * pb_x

                     + 1.5 * pa_yy * fx * fx * pb_xz

                     + 0.75 * pa_yz * fx * fx * pb_xy

                     + 1.5 * pa_y * fx * fx * pb_xyz

                     + 0.75 * fx * fx * pa_z * pb_xzz

                     + 0.5 * pa_yyyz * pb_xy * fx

                     + pa_yyy * fx * pb_xyz

                     + 1.5 * pa_yyz * fx * pb_xzz

                     + 1.5 * pa_yz * fx * pb_xyzz

                     + pa_yyyz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xyzz_r_0(double fga,
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
                                     double pb_x,
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_xyzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pa_z * pb_x * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pa_z * pb_x

                     - 1.5 * fx * fx * fz * fga * pb_xz

                     - 1.5 * pa_yyz * fx * pb_x * fz * fgb

                     + 3.75 * fx * fx * fx * fz * pa_z * pb_x

                     + 7.5 * fx * fx * fx * fz * pb_xz

                     + 9.0 * pa_yyz * fx * fx * fz * pb_x

                     + 18.0 * pa_yy * fx * fx * fz * pb_xz

                     - 1.5 * pa_yz * fx * pb_xy * fz * fgb

                     - 1.5 * pa_yz * fz * fga * pb_xy * fx

                     - 3.0 * pa_y * fz * fga * fx * pb_xyz

                     - 1.5 * fx * fz * fga * pa_z * pb_xzz

                     - pa_yyyz * pb_xy * fz * fgb

                     + 9.0 * pa_yz * fx * fx * fz * pb_xy

                     + 18.0 * pa_y * fx * fx * fz * pb_xyz

                     + 9.0 * fx * fx * fz * pa_z * pb_xzz

                     + 7.0 * pa_yyyz * fz * pb_xy * fx

                     + 14.0 * pa_yyy * fz * fx * pb_xyz

                     + 21.0 * pa_yyz * fx * fz * pb_xzz

                     - 3.0 * pa_yz * fz * fga * pb_xyzz

                     + 21.0 * pa_yz * fx * fz * pb_xyzz

                     + 16.0 * pa_yyyz * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xzzz_s_0(double fx,
                                     double pa_y,
                                     double pa_yyy,
                                     double pa_yyyz,
                                     double pa_yz,
                                     double pb_x,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_xzzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_y * fx * fx * fx * pb_x

                     + 0.75 * pa_yyy * fx * fx * pb_x

                     + 2.25 * pa_yz * fx * fx * pb_xz

                     + 2.25 * pa_y * fx * fx * pb_xzz

                     + 1.5 * pa_yyyz * pb_xz * fx

                     + 1.5 * pa_yyy * fx * pb_xzz

                     + 1.5 * pa_yz * fx * pb_xzzz

                     + pa_yyyz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xzzz_r_0(double fga,
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
                                     double pb_xzzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * pa_y * fx * fx * pb_x * fz * fgb

                     - 2.25 * pa_y * fz * fga * fx * fx * pb_x

                     - 1.5 * pa_yyy * fx * pb_x * fz * fgb

                     + 11.25 * pa_y * fx * fx * fx * fz * pb_x

                     + 9.0 * pa_yyy * fz * fx * fx * pb_x

                     - 4.5 * pa_yz * fx * pb_xz * fz * fgb

                     - 4.5 * pa_yz * fz * fga * pb_xz * fx

                     - 4.5 * pa_y * fz * fga * fx * pb_xzz

                     - 3.0 * pa_yyyz * pb_xz * fz * fgb

                     + 27.0 * pa_yz * fx * fx * fz * pb_xz

                     + 27.0 * pa_y * fx * fx * fz * pb_xzz

                     + 21.0 * pa_yyyz * fz * pb_xz * fx

                     + 21.0 * pa_yyy * fz * fx * pb_xzz

                     - 3.0 * pa_yz * fz * fga * pb_xzzz

                     + 21.0 * pa_yz * fx * fz * pb_xzzz

                     + 16.0 * pa_yyyz * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_yyyy_s_0(double fx,
                                     double pa_yyyz,
                                     double pa_yyz,
                                     double pa_yz,
                                     double pa_z,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double pb_yyyy,
                                     double s_0_0)
    {
        return s_0_0 * (5.625 * pa_yz * fx * fx * fx

                     + 7.5 * fx * fx * fx * pa_z * pb_y

                     + 0.75 * pa_yyyz * fx * fx

                     + 9.0 * pa_yyz * fx * fx * pb_y

                     + 13.5 * pa_yz * fx * fx * pb_yy

                     + 3.0 * fx * fx * pa_z * pb_yyy

                     + 3.0 * pa_yyyz * pb_yy * fx

                     + 6.0 * pa_yyz * fx * pb_yyy

                     + 1.5 * pa_yz * fx * pb_yyyy

                     + pa_yyyz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_yyyy_r_0(double fga,
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
                                     double pb_yyyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 13.5 * pa_yz * fx * fx * fz * fgb

                     + 56.25 * pa_yz * fx * fx * fx * fz

                     + 75.0 * fx * fx * fx * fz * pa_z * pb_y

                     - 2.25 * pa_yz * fz * fga * fx * fx

                     - 9.0 * fx * fx * pa_z * pb_y * fz * fgb

                     - 9.0 * fx * fx * fz * fga * pa_z * pb_y

                     - 3.0 * pa_yyyz * fx * fz * fgb

                     - 18.0 * pa_yyz * fx * pb_y * fz * fgb

                     + 9.0 * pa_yyyz * fz * fx * fx

                     + 108.0 * pa_yyz * fx * fx * fz * pb_y

                     + 162.0 * pa_yz * fx * fx * fz * pb_yy

                     - 9.0 * pa_yz * fx * pb_yy * fz * fgb

                     - 9.0 * pa_yz * fz * fga * pb_yy * fx

                     - 6.0 * fx * fz * fga * pa_z * pb_yyy

                     - 6.0 * pa_yyyz * pb_yy * fz * fgb

                     + 36.0 * fx * fx * fz * pa_z * pb_yyy

                     + 42.0 * pa_yyyz * fz * pb_yy * fx

                     + 84.0 * pa_yyz * fx * fz * pb_yyy

                     - 3.0 * pa_yz * fz * fga * pb_yyyy

                     + 21.0 * pa_yz * fx * fz * pb_yyyy

                     + 16.0 * pa_yyyz * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_yyyz_s_0(double fx,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyy,
                                     double pa_yyyz,
                                     double pa_yyz,
                                     double pa_yz,
                                     double pa_z,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double pb_yyyz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.9375 * fx * fx * fx * fx

                     + 1.125 * pa_yy * fx * fx * fx

                     + 3.375 * pa_y * fx * fx * fx * pb_y

                     + 1.875 * fx * fx * fx * pa_z * pb_z

                     + 1.125 * fx * fx * fx * pb_yy

                     + 0.75 * pa_yyy * fx * fx * pb_y

                     + 2.25 * pa_yyz * fx * fx * pb_z

                     + 2.25 * pa_yy * fx * fx * pb_yy

                     + 6.75 * pa_yz * fx * fx * pb_yz

                     + 0.75 * pa_y * fx * fx * pb_yyy

                     + 2.25 * fx * fx * pa_z * pb_yyz

                     + 1.5 * pa_yyyz * pb_yz * fx

                     + 0.5 * pa_yyy * fx * pb_yyy

                     + 4.5 * pa_yyz * fx * pb_yyz

                     + 1.5 * pa_yz * fx * pb_yyyz

                     + pa_yyyz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_yyyz_r_0(double fga,
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
                                     double pb_yyy,
                                     double pb_yyyz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (7.5 * fx * fx * fx * fx * fz

                     - 1.125 * fx * fx * fx * fz * fgb

                     - 1.125 * fx * fx * fx * fz * fga

                     - 2.25 * pa_yy * fx * fx * fz * fgb

                     + 11.25 * pa_yy * fx * fx * fx * fz

                     + 33.75 * pa_y * fx * fx * fx * fz * pb_y

                     + 18.75 * fx * fx * fx * fz * pa_z * pb_z

                     - 2.25 * pa_y * fx * fx * pb_y * fz * fgb

                     - 2.25 * pa_y * fz * fga * fx * fx * pb_y

                     - 2.25 * fx * fx * pa_z * fz * fgb * pb_z

                     - 2.25 * fx * fx * fz * fga * pa_z * pb_z

                     - 2.25 * fx * fx * fz * fga * pb_yy

                     - 1.5 * pa_yyy * fx * pb_y * fz * fgb

                     - 4.5 * pa_yyz * fx * fz * fgb * pb_z

                     + 11.25 * fx * fx * fx * fz * pb_yy

                     + 9.0 * pa_yyy * fz * fx * fx * pb_y

                     + 27.0 * pa_yyz * fx * fx * fz * pb_z

                     + 27.0 * pa_yy * fx * fx * fz * pb_yy

                     + 81.0 * pa_yz * fx * fx * fz * pb_yz

                     - 4.5 * pa_yz * fx * pb_yz * fz * fgb

                     - 4.5 * pa_yz * fz * fga * pb_yz * fx

                     - 1.5 * pa_y * fz * fga * fx * pb_yyy

                     - 4.5 * fx * fz * fga * pa_z * pb_yyz

                     - 3.0 * pa_yyyz * pb_yz * fz * fgb

                     + 9.0 * pa_y * fx * fx * fz * pb_yyy

                     + 27.0 * fx * fx * fz * pa_z * pb_yyz

                     + 21.0 * pa_yyyz * fz * pb_yz * fx

                     + 7.0 * pa_yyy * fz * fx * pb_yyy

                     + 63.0 * pa_yyz * fx * fz * pb_yyz

                     - 3.0 * pa_yz * fz * fga * pb_yyyz

                     + 21.0 * pa_yz * fx * fz * pb_yyyz

                     + 16.0 * pa_yyyz * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_yyzz_s_0(double fx,
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
                                     double pb_yyzz,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_yz * fx * fx * fx

                     + 2.25 * pa_y * fx * fx * fx * pb_z

                     + 0.75 * fx * fx * fx * pa_z * pb_y

                     + 1.5 * fx * fx * fx * pb_yz

                     + 0.25 * pa_yyyz * fx * fx

                     + 0.5 * pa_yyy * fx * fx * pb_z

                     + 1.5 * pa_yyz * fx * fx * pb_y

                     + 3.0 * pa_yy * fx * fx * pb_yz

                     + 2.25 * pa_yz * fx * fx * pb_zz

                     + 0.75 * pa_yz * fx * fx * pb_yy

                     + 1.5 * pa_y * fx * fx * pb_yyz

                     + 1.5 * fx * fx * pa_z * pb_yzz

                     + 0.5 * pa_yyyz * pb_yy * fx

                     + 0.5 * pa_yyyz * fx * pb_zz

                     + pa_yyy * fx * pb_yyz

                     + 3.0 * pa_yyz * fx * pb_yzz

                     + 1.5 * pa_yz * fx * pb_yyzz

                     + pa_yyyz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_yyzz_r_0(double fga,
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
                                     double pb_yyzz,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_yz * fx * fx * fz * fgb

                     + 11.25 * pa_yz * fx * fx * fx * fz

                     + 22.5 * pa_y * fx * fx * fx * fz * pb_z

                     - 1.5 * pa_y * fx * fx * fz * fgb * pb_z

                     - 0.75 * pa_yz * fz * fga * fx * fx

                     - 1.5 * pa_y * fz * fga * fx * fx * pb_z

                     - 1.5 * fx * fx * pa_z * pb_y * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pa_z * pb_y

                     - 3.0 * fx * fx * fz * fga * pb_yz

                     - pa_yyyz * fx * fz * fgb

                     - pa_yyy * fx * fz * fgb * pb_z

                     - 3.0 * pa_yyz * fx * pb_y * fz * fgb

                     + 7.5 * fx * fx * fx * fz * pa_z * pb_y

                     + 15.0 * fx * fx * fx * fz * pb_yz

                     + 3.0 * pa_yyyz * fz * fx * fx

                     + 6.0 * pa_yyy * fz * fx * fx * pb_z

                     + 18.0 * pa_yyz * fx * fx * fz * pb_y

                     + 36.0 * pa_yy * fx * fx * fz * pb_yz

                     + 27.0 * pa_yz * fx * fx * fz * pb_zz

                     - 1.5 * pa_yz * fx * pb_yy * fz * fgb

                     - 1.5 * pa_yz * fx * fz * fgb * pb_zz

                     - 1.5 * pa_yz * fz * fga * pb_yy * fx

                     - 1.5 * pa_yz * fz * fga * fx * pb_zz

                     - 3.0 * pa_y * fz * fga * fx * pb_yyz

                     - 3.0 * fx * fz * fga * pa_z * pb_yzz

                     - pa_yyyz * pb_yy * fz * fgb

                     - pa_yyyz * fz * fgb * pb_zz

                     + 9.0 * pa_yz * fx * fx * fz * pb_yy

                     + 18.0 * pa_y * fx * fx * fz * pb_yyz

                     + 18.0 * fx * fx * fz * pa_z * pb_yzz

                     + 7.0 * pa_yyyz * fz * pb_yy * fx

                     + 7.0 * pa_yyyz * fz * fx * pb_zz

                     + 14.0 * pa_yyy * fz * fx * pb_yyz

                     + 42.0 * pa_yyz * fx * fz * pb_yzz

                     - 3.0 * pa_yz * fz * fga * pb_yyzz

                     + 21.0 * pa_yz * fx * fz * pb_yyzz

                     + 16.0 * pa_yyyz * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_yzzz_s_0(double fx,
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
                                     double pb_yzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.5625 * fx * fx * fx * fx

                     + 1.125 * pa_yy * fx * fx * fx

                     + 1.125 * pa_y * fx * fx * fx * pb_y

                     + 1.125 * fx * fx * fx * pa_z * pb_z

                     + 1.125 * fx * fx * fx * pb_zz

                     + 0.75 * pa_yyy * fx * fx * pb_y

                     + 2.25 * pa_yyz * fx * fx * pb_z

                     + 2.25 * pa_yy * fx * fx * pb_zz

                     + 2.25 * pa_yz * fx * fx * pb_yz

                     + 2.25 * pa_y * fx * fx * pb_yzz

                     + 0.75 * fx * fx * pa_z * pb_zzz

                     + 1.5 * pa_yyyz * pb_yz * fx

                     + 1.5 * pa_yyy * fx * pb_yzz

                     + 1.5 * pa_yyz * fx * pb_zzz

                     + 1.5 * pa_yz * fx * pb_yzzz

                     + pa_yyyz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_yzzz_r_0(double fga,
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
                                     double pb_yzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.125 * fx * fx * fx * fz * fgb

                     - 1.125 * fx * fx * fx * fz * fga

                     - 2.25 * pa_yy * fx * fx * fz * fgb

                     + 4.5 * fx * fx * fx * fx * fz

                     + 11.25 * pa_yy * fx * fx * fx * fz

                     - 2.25 * pa_y * fx * fx * pb_y * fz * fgb

                     - 2.25 * pa_y * fz * fga * fx * fx * pb_y

                     - 2.25 * fx * fx * pa_z * pb_z * fz * fgb

                     - 2.25 * fx * fx * fz * fga * pa_z * pb_z

                     - 2.25 * fx * fx * fz * fga * pb_zz

                     - 1.5 * pa_yyy * fx * pb_y * fz * fgb

                     - 4.5 * pa_yyz * fx * pb_z * fz * fgb

                     + 11.25 * pa_y * fx * fx * fx * fz * pb_y

                     + 11.25 * fx * fx * fx * fz * pa_z * pb_z

                     + 11.25 * fx * fx * fx * fz * pb_zz

                     + 9.0 * pa_yyy * fz * fx * fx * pb_y

                     + 27.0 * pa_yyz * fx * fx * fz * pb_z

                     + 27.0 * pa_yy * fx * fx * fz * pb_zz

                     - 4.5 * pa_yz * fx * pb_yz * fz * fgb

                     - 4.5 * pa_yz * fz * fga * pb_yz * fx

                     - 4.5 * pa_y * fz * fga * fx * pb_yzz

                     - 1.5 * fx * fz * fga * pa_z * pb_zzz

                     - 3.0 * pa_yyyz * pb_yz * fz * fgb

                     + 27.0 * pa_yz * fx * fx * fz * pb_yz

                     + 27.0 * pa_y * fx * fx * fz * pb_yzz

                     + 9.0 * fx * fx * fz * pa_z * pb_zzz

                     + 21.0 * pa_yyyz * fz * pb_yz * fx

                     + 21.0 * pa_yyy * fz * fx * pb_yzz

                     + 21.0 * pa_yyz * fx * fz * pb_zzz

                     - 3.0 * pa_yz * fz * fga * pb_yzzz

                     + 21.0 * pa_yz * fx * fz * pb_yzzz

                     + 16.0 * pa_yyyz * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_zzzz_s_0(double fx,
                                     double pa_y,
                                     double pa_yyy,
                                     double pa_yyyz,
                                     double pa_yz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double pb_zzzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_yz * fx * fx * fx

                     + 4.5 * pa_y * fx * fx * fx * pb_z

                     + 0.75 * pa_yyyz * fx * fx

                     + 3.0 * pa_yyy * fx * fx * pb_z

                     + 4.5 * pa_yz * fx * fx * pb_zz

                     + 3.0 * pa_y * fx * fx * pb_zzz

                     + 3.0 * pa_yyyz * pb_zz * fx

                     + 2.0 * pa_yyy * fx * pb_zzz

                     + 1.5 * pa_yz * fx * pb_zzzz

                     + pa_yyyz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_zzzz_r_0(double fga,
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
                                     double pb_zzzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_yz * fx * fx * fz * fgb

                     - 9.0 * pa_y * fx * fx * pb_z * fz * fgb

                     - 2.25 * pa_yz * fz * fga * fx * fx

                     - 9.0 * pa_y * fz * fga * fx * fx * pb_z

                     - 3.0 * pa_yyyz * fx * fz * fgb

                     - 6.0 * pa_yyy * fx * pb_z * fz * fgb

                     + 11.25 * pa_yz * fx * fx * fx * fz

                     + 45.0 * pa_y * fx * fx * fx * fz * pb_z

                     + 9.0 * pa_yyyz * fz * fx * fx

                     + 36.0 * pa_yyy * fz * fx * fx * pb_z

                     - 9.0 * pa_yz * fx * pb_zz * fz * fgb

                     - 9.0 * pa_yz * fz * fga * pb_zz * fx

                     - 6.0 * pa_y * fz * fga * fx * pb_zzz

                     - 6.0 * pa_yyyz * pb_zz * fz * fgb

                     + 54.0 * pa_yz * fx * fx * fz * pb_zz

                     + 36.0 * pa_y * fx * fx * fz * pb_zzz

                     + 42.0 * pa_yyyz * fz * pb_zz * fx

                     + 28.0 * pa_yyy * fz * fx * pb_zzz

                     - 3.0 * pa_yz * fz * fga * pb_zzzz

                     + 21.0 * pa_yz * fx * fz * pb_zzzz

                     + 16.0 * pa_yyyz * fz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xxxx_s_0(double fx,
                                     double pa_yy,
                                     double pa_yyzz,
                                     double pa_zz,
                                     double pb_xx,
                                     double pb_xxxx,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * pa_yy * fx * fx * fx

                     + 0.375 * fx * fx * fx * pa_zz

                     + 0.75 * pa_yyzz * fx * fx

                     + 0.75 * fx * fx * fx * pb_xx

                     + 1.5 * pa_yy * fx * fx * pb_xx

                     + 1.5 * fx * fx * pa_zz * pb_xx

                     + 3.0 * pa_yyzz * pb_xx * fx

                     + 0.25 * fx * fx * pb_xxxx

                     + 0.5 * pa_yy * fx * pb_xxxx

                     + 0.5 * fx * pa_zz * pb_xxxx

                     + pa_yyzz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xxxx_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_yy,
                                     double pa_yyzz,
                                     double pa_zz,
                                     double pb_xx,
                                     double pb_xxxx,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fx * fz * fgb

                     - 0.75 * fx * fx * fx * fz * fga

                     - 1.5 * pa_yy * fx * fx * fz * fgb

                     - 0.75 * pa_yy * fz * fga * fx * fx

                     - 1.5 * fx * fx * pa_zz * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     - 0.75 * fz * fga * pa_zz * fx * fx

                     - 3.0 * pa_yyzz * fx * fz * fgb

                     + 3.75 * pa_yy * fz * fx * fx * fx

                     + 3.75 * fx * fx * fx * pa_zz * fz

                     + 9.0 * pa_yyzz * fz * fx * fx

                     - 1.5 * fx * fx * pb_xx * fz * fgb

                     - 3.0 * fx * fx * fz * fga * pb_xx

                     - 3.0 * pa_yy * fx * pb_xx * fz * fgb

                     - 3.0 * pa_yy * fz * fga * pb_xx * fx

                     - 3.0 * fx * pa_zz * pb_xx * fz * fgb

                     + 7.5 * fx * fx * fx * fz * pb_xx

                     - 3.0 * fz * fga * pa_zz * pb_xx * fx

                     - 6.0 * pa_yyzz * pb_xx * fz * fgb

                     + 18.0 * pa_yy * fz * fx * fx * pb_xx

                     + 18.0 * fx * fx * pa_zz * fz * pb_xx

                     + 42.0 * pa_yyzz * fz * pb_xx * fx

                     - fx * fz * fga * pb_xxxx

                     - pa_yy * fz * fga * pb_xxxx

                     + 3.0 * fx * fx * fz * pb_xxxx

                     - fz * fga * pa_zz * pb_xxxx

                     + 7.0 * pa_yy * fz * fx * pb_xxxx

                     + 7.0 * fx * pa_zz * fz * pb_xxxx

                     + 16.0 * pa_yyzz * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xxxy_s_0(double fx,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyzz,
                                     double pa_yzz,
                                     double pa_zz,
                                     double pb_x,
                                     double pb_xxx,
                                     double pb_xxxy,
                                     double pb_xy,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx * fx * pb_x

                     + 1.5 * pa_yzz * fx * fx * pb_x

                     + 0.375 * fx * fx * fx * pb_xy

                     + 0.75 * pa_yy * fx * fx * pb_xy

                     + 0.5 * pa_y * fx * fx * pb_xxx

                     + 0.75 * fx * fx * pa_zz * pb_xy

                     + 1.5 * pa_yyzz * pb_xy * fx

                     + pa_yzz * fx * pb_xxx

                     + 0.25 * fx * fx * pb_xxxy

                     + 0.5 * pa_yy * fx * pb_xxxy

                     + 0.5 * fx * pa_zz * pb_xxxy

                     + pa_yyzz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xxxy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyzz,
                                     double pa_yzz,
                                     double pa_zz,
                                     double pb_x,
                                     double pb_xxx,
                                     double pb_xxxy,
                                     double pb_xy,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_y * fx * fx * pb_x * fz * fgb

                     - 1.5 * pa_y * fx * fx * fz * fga * pb_x

                     - 3.0 * pa_yzz * fx * pb_x * fz * fgb

                     + 7.5 * pa_y * fx * fx * fx * fz * pb_x

                     + 18.0 * pa_yzz * fx * fx * fz * pb_x

                     - 0.75 * fx * fx * pb_xy * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_xy

                     - 1.5 * pa_yy * fx * pb_xy * fz * fgb

                     - 1.5 * pa_yy * fz * fga * pb_xy * fx

                     - pa_y * fx * fz * fga * pb_xxx

                     - 1.5 * fx * pa_zz * pb_xy * fz * fgb

                     + 3.75 * fx * fx * fx * fz * pb_xy

                     - 1.5 * fz * fga * pa_zz * pb_xy * fx

                     - 3.0 * pa_yyzz * pb_xy * fz * fgb

                     + 9.0 * pa_yy * fz * fx * fx * pb_xy

                     + 6.0 * pa_y * fx * fx * fz * pb_xxx

                     + 9.0 * fx * fx * pa_zz * fz * pb_xy

                     + 21.0 * pa_yyzz * fz * pb_xy * fx

                     + 14.0 * pa_yzz * fx * fz * pb_xxx

                     - fx * fz * fga * pb_xxxy

                     - pa_yy * fz * fga * pb_xxxy

                     + 3.0 * fx * fx * fz * pb_xxxy

                     - fz * fga * pa_zz * pb_xxxy

                     + 7.0 * pa_yy * fz * fx * pb_xxxy

                     + 7.0 * fx * pa_zz * fz * pb_xxxy

                     + 16.0 * pa_yyzz * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xxxz_s_0(double fx,
                                     double pa_yy,
                                     double pa_yyz,
                                     double pa_yyzz,
                                     double pa_z,
                                     double pa_zz,
                                     double pb_x,
                                     double pb_xxx,
                                     double pb_xxxz,
                                     double pb_xz,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * fx * pa_z * pb_x

                     + 1.5 * pa_yyz * fx * fx * pb_x

                     + 0.375 * fx * fx * fx * pb_xz

                     + 0.75 * pa_yy * fx * fx * pb_xz

                     + 0.75 * fx * fx * pa_zz * pb_xz

                     + 0.5 * fx * fx * pa_z * pb_xxx

                     + 1.5 * pa_yyzz * pb_xz * fx

                     + pa_yyz * fx * pb_xxx

                     + 0.25 * fx * fx * pb_xxxz

                     + 0.5 * pa_yy * fx * pb_xxxz

                     + 0.5 * fx * pa_zz * pb_xxxz

                     + pa_yyzz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xxxz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_yy,
                                     double pa_yyz,
                                     double pa_yyzz,
                                     double pa_z,
                                     double pa_zz,
                                     double pb_x,
                                     double pb_xxx,
                                     double pb_xxxz,
                                     double pb_xz,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_z * pb_x * fz * fgb

                     - 1.5 * fz * fga * pa_z * fx * fx * pb_x

                     - 3.0 * pa_yyz * fx * pb_x * fz * fgb

                     + 7.5 * fx * fx * fx * pa_z * fz * pb_x

                     + 18.0 * pa_yyz * fz * fx * fx * pb_x

                     - 0.75 * fx * fx * pb_xz * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_xz

                     - 1.5 * pa_yy * fx * pb_xz * fz * fgb

                     - 1.5 * pa_yy * fz * fga * pb_xz * fx

                     - 1.5 * fx * pa_zz * pb_xz * fz * fgb

                     + 3.75 * fx * fx * fx * fz * pb_xz

                     - 1.5 * fz * fga * pa_zz * pb_xz * fx

                     - fz * fga * pa_z * fx * pb_xxx

                     - 3.0 * pa_yyzz * pb_xz * fz * fgb

                     + 9.0 * pa_yy * fz * fx * fx * pb_xz

                     + 9.0 * fx * fx * pa_zz * fz * pb_xz

                     + 6.0 * fx * fx * pa_z * fz * pb_xxx

                     + 21.0 * pa_yyzz * fz * pb_xz * fx

                     + 14.0 * pa_yyz * fz * fx * pb_xxx

                     - fx * fz * fga * pb_xxxz

                     - pa_yy * fz * fga * pb_xxxz

                     + 3.0 * fx * fx * fz * pb_xxxz

                     - fz * fga * pa_zz * pb_xxxz

                     + 7.0 * pa_yy * fz * fx * pb_xxxz

                     + 7.0 * fx * pa_zz * fz * pb_xxxz

                     + 16.0 * pa_yyzz * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xxyy_s_0(double fx,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyzz,
                                     double pa_yzz,
                                     double pa_zz,
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyy,
                                     double pb_y,
                                     double pb_yy,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * fx * fx * fx * pa_zz

                     + 0.125 * pa_yy * fx * fx * fx

                     + 0.5 * pa_y * fx * fx * fx * pb_y

                     + 0.375 * fx * fx * fx * pb_xx

                     + 0.25 * pa_yyzz * fx * fx

                     + pa_yzz * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_zz * pb_xx

                     + 0.125 * fx * fx * fx * pb_yy

                     + 0.25 * pa_yy * fx * fx * pb_xx

                     + 0.25 * pa_yy * fx * fx * pb_yy

                     + pa_y * fx * fx * pb_xxy

                     + 0.25 * fx * fx * pa_zz * pb_yy

                     + 0.5 * pa_yyzz * pb_xx * fx

                     + 0.5 * pa_yyzz * fx * pb_yy

                     + 2.0 * pa_yzz * fx * pb_xxy

                     + 0.25 * fx * fx * pb_xxyy

                     + 0.5 * pa_yy * fx * pb_xxyy

                     + 0.5 * fx * pa_zz * pb_xxyy

                     + pa_yyzz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xxyy_r_0(double fga,
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
                                     double pb_xxyy,
                                     double pb_y,
                                     double pb_yy,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * fx * fz * fgb

                     - 0.5 * fx * fx * fx * fz * fga

                     - fx * fx * pa_zz * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     + 3.75 * fx * fx * fx * pa_zz * fz

                     - 0.5 * pa_yy * fx * fx * fz * fgb

                     - 0.25 * pa_yy * fz * fga * fx * fx

                     - pa_y * fx * fx * fz * fgb * pb_y

                     - pa_y * fx * fx * fz * fga * pb_y

                     - fx * fx * fz * fga * pb_xx

                     - 0.25 * fz * fga * pa_zz * fx * fx

                     - pa_yyzz * fx * fz * fgb

                     + 1.25 * pa_yy * fz * fx * fx * fx

                     - 2.0 * pa_yzz * fx * fz * fgb * pb_y

                     + 5.0 * pa_y * fx * fx * fx * fz * pb_y

                     + 3.75 * fx * fx * fx * fz * pb_xx

                     + 3.0 * pa_yyzz * fz * fx * fx

                     + 12.0 * pa_yzz * fx * fx * fz * pb_y

                     + 9.0 * fx * fx * pa_zz * fz * pb_xx

                     - 0.25 * fx * fx * pb_xx * fz * fgb

                     - 0.25 * fx * fx * fz * fgb * pb_yy

                     - 0.5 * fx * fx * fz * fga * pb_yy

                     - 0.5 * pa_yy * fx * pb_xx * fz * fgb

                     - 0.5 * pa_yy * fx * fz * fgb * pb_yy

                     - 0.5 * pa_yy * fz * fga * pb_xx * fx

                     - 0.5 * pa_yy * fz * fga * fx * pb_yy

                     - 2.0 * pa_y * fx * fz * fga * pb_xxy

                     - 0.5 * fx * pa_zz * pb_xx * fz * fgb

                     - 0.5 * fx * pa_zz * fz * fgb * pb_yy

                     + 1.25 * fx * fx * fx * fz * pb_yy

                     - 0.5 * fz * fga * pa_zz * pb_xx * fx

                     - 0.5 * fz * fga * pa_zz * fx * pb_yy

                     - pa_yyzz * pb_xx * fz * fgb

                     - pa_yyzz * fz * fgb * pb_yy

                     + 3.0 * pa_yy * fz * fx * fx * pb_xx

                     + 3.0 * pa_yy * fz * fx * fx * pb_yy

                     + 12.0 * pa_y * fx * fx * fz * pb_xxy

                     + 3.0 * fx * fx * pa_zz * fz * pb_yy

                     + 7.0 * pa_yyzz * fz * pb_xx * fx

                     + 7.0 * pa_yyzz * fz * fx * pb_yy

                     + 28.0 * pa_yzz * fx * fz * pb_xxy

                     - fx * fz * fga * pb_xxyy

                     - pa_yy * fz * fga * pb_xxyy

                     + 3.0 * fx * fx * fz * pb_xxyy

                     - fz * fga * pa_zz * pb_xxyy

                     + 7.0 * pa_yy * fz * fx * pb_xxyy

                     + 7.0 * fx * pa_zz * fz * pb_xxyy

                     + 16.0 * pa_yyzz * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xxyz_s_0(double fx,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyz,
                                     double pa_yyzz,
                                     double pa_yz,
                                     double pa_yzz,
                                     double pa_z,
                                     double pa_zz,
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyz,
                                     double pb_xxz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.5 * pa_yz * fx * fx * fx

                     + 0.25 * pa_y * fx * fx * fx * pb_z

                     + 0.25 * fx * fx * fx * pa_z * pb_y

                     + 0.5 * pa_yyz * fx * fx * pb_y

                     + 0.5 * pa_yzz * fx * fx * pb_z

                     + pa_yz * fx * fx * pb_xx

                     + 0.125 * fx * fx * fx * pb_yz

                     + 0.25 * pa_yy * fx * fx * pb_yz

                     + 0.5 * pa_y * fx * fx * pb_xxz

                     + 0.25 * fx * fx * pa_zz * pb_yz

                     + 0.5 * fx * fx * pa_z * pb_xxy

                     + 0.5 * pa_yyzz * fx * pb_yz

                     + pa_yyz * fx * pb_xxy

                     + pa_yzz * fx * pb_xxz

                     + 0.25 * fx * fx * pb_xxyz

                     + 0.5 * pa_yy * fx * pb_xxyz

                     + 0.5 * fx * pa_zz * pb_xxyz

                     + pa_yyzz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xxyz_r_0(double fga,
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
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyz,
                                     double pb_xxz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- pa_yz * fx * fx * fz * fgb

                     + 5.0 * pa_yz * fx * fx * fx * fz

                     - 0.5 * pa_y * fx * fx * fz * fgb * pb_z

                     - 0.5 * pa_y * fx * fx * fz * fga * pb_z

                     - 0.5 * fx * fx * pa_z * fz * fgb * pb_y

                     - 0.5 * fz * fga * pa_z * fx * fx * pb_y

                     - pa_yyz * fx * fz * fgb * pb_y

                     - pa_yzz * fx * fz * fgb * pb_z

                     + 2.5 * pa_y * fx * fx * fx * fz * pb_z

                     + 2.5 * fx * fx * fx * pa_z * fz * pb_y

                     + 6.0 * pa_yyz * fz * fx * fx * pb_y

                     + 6.0 * pa_yzz * fx * fx * fz * pb_z

                     + 12.0 * pa_yz * fx * fx * fz * pb_xx

                     - 0.25 * fx * fx * fz * fgb * pb_yz

                     - 0.5 * fx * fx * fz * fga * pb_yz

                     - 0.5 * pa_yy * fx * fz * fgb * pb_yz

                     - 0.5 * pa_yy * fz * fga * fx * pb_yz

                     - pa_y * fx * fz * fga * pb_xxz

                     - 0.5 * fx * pa_zz * fz * fgb * pb_yz

                     + 1.25 * fx * fx * fx * fz * pb_yz

                     - 0.5 * fz * fga * pa_zz * fx * pb_yz

                     - fz * fga * pa_z * fx * pb_xxy

                     - pa_yyzz * fz * fgb * pb_yz

                     + 3.0 * pa_yy * fz * fx * fx * pb_yz

                     + 6.0 * pa_y * fx * fx * fz * pb_xxz

                     + 3.0 * fx * fx * pa_zz * fz * pb_yz

                     + 6.0 * fx * fx * pa_z * fz * pb_xxy

                     + 7.0 * pa_yyzz * fz * fx * pb_yz

                     + 14.0 * pa_yyz * fz * fx * pb_xxy

                     + 14.0 * pa_yzz * fx * fz * pb_xxz

                     - fx * fz * fga * pb_xxyz

                     - pa_yy * fz * fga * pb_xxyz

                     + 3.0 * fx * fx * fz * pb_xxyz

                     - fz * fga * pa_zz * pb_xxyz

                     + 7.0 * pa_yy * fz * fx * pb_xxyz

                     + 7.0 * fx * pa_zz * fz * pb_xxyz

                     + 16.0 * pa_yyzz * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xxzz_s_0(double fx,
                                     double pa_yy,
                                     double pa_yyz,
                                     double pa_yyzz,
                                     double pa_z,
                                     double pa_zz,
                                     double pb_xx,
                                     double pb_xxz,
                                     double pb_xxzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * pa_yy * fx * fx * fx

                     + 0.125 * fx * fx * fx * pa_zz

                     + 0.5 * fx * fx * fx * pa_z * pb_z

                     + 0.375 * fx * fx * fx * pb_xx

                     + 0.25 * pa_yyzz * fx * fx

                     + pa_yyz * fx * fx * pb_z

                     + 0.75 * pa_yy * fx * fx * pb_xx

                     + 0.125 * fx * fx * fx * pb_zz

                     + 0.25 * pa_yy * fx * fx * pb_zz

                     + 0.25 * fx * fx * pa_zz * pb_xx

                     + 0.25 * fx * fx * pa_zz * pb_zz

                     + fx * fx * pa_z * pb_xxz

                     + 0.5 * pa_yyzz * pb_xx * fx

                     + 0.5 * pa_yyzz * fx * pb_zz

                     + 2.0 * pa_yyz * fx * pb_xxz

                     + 0.25 * fx * fx * pb_xxzz

                     + 0.5 * pa_yy * fx * pb_xxzz

                     + 0.5 * fx * pa_zz * pb_xxzz

                     + pa_yyzz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xxzz_r_0(double fga,
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
                                     double pb_xxzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * fx * fz * fgb

                     - 0.5 * fz * fga * fx * fx * fx

                     - pa_yy * fx * fx * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     + 3.75 * pa_yy * fx * fx * fx * fz

                     - 0.25 * pa_yy * fz * fga * fx * fx

                     - 0.5 * fx * fx * pa_zz * fz * fgb

                     - fx * fx * pa_z * fz * fgb * pb_z

                     - 0.25 * fz * fga * pa_zz * fx * fx

                     - fz * fga * pa_z * fx * fx * pb_z

                     - fz * fga * fx * fx * pb_xx

                     - pa_yyzz * fx * fz * fgb

                     - 2.0 * pa_yyz * fx * fz * fgb * pb_z

                     + 1.25 * fx * fx * fx * pa_zz * fz

                     + 5.0 * fx * fx * fx * pa_z * fz * pb_z

                     + 3.75 * fx * fx * fx * fz * pb_xx

                     + 3.0 * pa_yyzz * fz * fx * fx

                     + 12.0 * pa_yyz * fz * fx * fx * pb_z

                     + 9.0 * pa_yy * fx * fx * fz * pb_xx

                     - 0.25 * fx * fx * pb_xx * fz * fgb

                     - 0.25 * fx * fx * fz * fgb * pb_zz

                     - 0.5 * fx * fx * fz * fga * pb_zz

                     - 0.5 * pa_yy * fx * pb_xx * fz * fgb

                     - 0.5 * pa_yy * fx * fz * fgb * pb_zz

                     - 0.5 * pa_yy * fz * fga * pb_xx * fx

                     - 0.5 * pa_yy * fz * fga * fx * pb_zz

                     - 0.5 * fx * pa_zz * pb_xx * fz * fgb

                     - 0.5 * fx * pa_zz * fz * fgb * pb_zz

                     + 1.25 * fx * fx * fx * fz * pb_zz

                     - 0.5 * fz * fga * pa_zz * pb_xx * fx

                     - 0.5 * fz * fga * pa_zz * fx * pb_zz

                     - 2.0 * fz * fga * pa_z * fx * pb_xxz

                     - pa_yyzz * pb_xx * fz * fgb

                     - pa_yyzz * fz * fgb * pb_zz

                     + 3.0 * pa_yy * fz * fx * fx * pb_zz

                     + 3.0 * fx * fx * pa_zz * fz * pb_xx

                     + 3.0 * fx * fx * pa_zz * fz * pb_zz

                     + 12.0 * fx * fx * pa_z * fz * pb_xxz

                     + 7.0 * pa_yyzz * fz * pb_xx * fx

                     + 7.0 * pa_yyzz * fz * fx * pb_zz

                     + 28.0 * pa_yyz * fz * fx * pb_xxz

                     - fx * fz * fga * pb_xxzz

                     - pa_yy * fz * fga * pb_xxzz

                     + 3.0 * fx * fx * fz * pb_xxzz

                     - fz * fga * pa_zz * pb_xxzz

                     + 7.0 * pa_yy * fz * fx * pb_xxzz

                     + 7.0 * fx * pa_zz * fz * pb_xxzz

                     + 16.0 * pa_yyzz * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xyyy_s_0(double fx,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyzz,
                                     double pa_yzz,
                                     double pa_zz,
                                     double pb_x,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_xyyy,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx * fx * pb_x

                     + 1.125 * fx * fx * fx * pb_xy

                     + 1.5 * pa_yzz * fx * fx * pb_x

                     + 2.25 * fx * fx * pa_zz * pb_xy

                     + 0.75 * pa_yy * fx * fx * pb_xy

                     + 1.5 * pa_y * fx * fx * pb_xyy

                     + 1.5 * pa_yyzz * pb_xy * fx

                     + 3.0 * pa_yzz * fx * pb_xyy

                     + 0.25 * fx * fx * pb_xyyy

                     + 0.5 * pa_yy * fx * pb_xyyy

                     + 0.5 * fx * pa_zz * pb_xyyy

                     + pa_yyzz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xyyy_r_0(double fga,
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
                                     double pb_xyyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_y * fx * fx * pb_x * fz * fgb

                     - 1.5 * pa_y * fx * fx * fz * fga * pb_x

                     - 3.0 * fx * fx * fz * fga * pb_xy

                     - 3.0 * pa_yzz * fx * pb_x * fz * fgb

                     + 7.5 * pa_y * fx * fx * fx * fz * pb_x

                     + 11.25 * fx * fx * fx * fz * pb_xy

                     + 18.0 * pa_yzz * fx * fx * fz * pb_x

                     + 27.0 * fx * fx * pa_zz * fz * pb_xy

                     - 0.75 * fx * fx * pb_xy * fz * fgb

                     - 1.5 * pa_yy * fx * pb_xy * fz * fgb

                     - 1.5 * pa_yy * fz * fga * pb_xy * fx

                     - 3.0 * pa_y * fx * fz * fga * pb_xyy

                     - 1.5 * fx * pa_zz * pb_xy * fz * fgb

                     - 1.5 * fz * fga * pa_zz * pb_xy * fx

                     - 3.0 * pa_yyzz * pb_xy * fz * fgb

                     + 9.0 * pa_yy * fz * fx * fx * pb_xy

                     + 18.0 * pa_y * fx * fx * fz * pb_xyy

                     + 21.0 * pa_yyzz * fz * pb_xy * fx

                     + 42.0 * pa_yzz * fx * fz * pb_xyy

                     - fx * fz * fga * pb_xyyy

                     - pa_yy * fz * fga * pb_xyyy

                     + 3.0 * fx * fx * fz * pb_xyyy

                     - fz * fga * pa_zz * pb_xyyy

                     + 7.0 * pa_yy * fz * fx * pb_xyyy

                     + 7.0 * fx * pa_zz * fz * pb_xyyy

                     + 16.0 * pa_yyzz * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xyyz_s_0(double fx,
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
                                     double pb_xyy,
                                     double pb_xyyz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * fx * pa_z * pb_x

                     + 0.375 * fx * fx * fx * pb_xz

                     + 0.5 * pa_yyz * fx * fx * pb_x

                     + 2.0 * pa_yz * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_zz * pb_xz

                     + 0.25 * pa_yy * fx * fx * pb_xz

                     + pa_y * fx * fx * pb_xyz

                     + 0.5 * fx * fx * pa_z * pb_xyy

                     + 0.5 * pa_yyzz * pb_xz * fx

                     + pa_yyz * fx * pb_xyy

                     + 2.0 * pa_yzz * fx * pb_xyz

                     + 0.25 * fx * fx * pb_xyyz

                     + 0.5 * pa_yy * fx * pb_xyyz

                     + 0.5 * fx * pa_zz * pb_xyyz

                     + pa_yyzz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xyyz_r_0(double fga,
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
                                     double pb_x,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_xyyz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double r_0_0)
    {
        return r_0_0 * (7.5 * fx * fx * fx * pa_z * fz * pb_x

                     - 0.5 * fx * fx * pa_z * pb_x * fz * fgb

                     - fx * fx * fz * fga * pb_xz

                     - 0.5 * fz * fga * pa_z * fx * fx * pb_x

                     - pa_yyz * fx * pb_x * fz * fgb

                     + 3.75 * fx * fx * fx * fz * pb_xz

                     + 6.0 * pa_yyz * fz * fx * fx * pb_x

                     + 24.0 * pa_yz * fx * fx * fz * pb_xy

                     + 9.0 * fx * fx * pa_zz * fz * pb_xz

                     - 0.25 * fx * fx * pb_xz * fz * fgb

                     - 0.5 * pa_yy * fx * pb_xz * fz * fgb

                     - 0.5 * pa_yy * fz * fga * pb_xz * fx

                     - 2.0 * pa_y * fx * fz * fga * pb_xyz

                     - 0.5 * fx * pa_zz * pb_xz * fz * fgb

                     - 0.5 * fz * fga * pa_zz * pb_xz * fx

                     - fz * fga * pa_z * fx * pb_xyy

                     - pa_yyzz * pb_xz * fz * fgb

                     + 3.0 * pa_yy * fz * fx * fx * pb_xz

                     + 12.0 * pa_y * fx * fx * fz * pb_xyz

                     + 6.0 * fx * fx * pa_z * fz * pb_xyy

                     + 7.0 * pa_yyzz * fz * pb_xz * fx

                     + 14.0 * pa_yyz * fz * fx * pb_xyy

                     + 28.0 * pa_yzz * fx * fz * pb_xyz

                     - fx * fz * fga * pb_xyyz

                     - pa_yy * fz * fga * pb_xyyz

                     + 3.0 * fx * fx * fz * pb_xyyz

                     - fz * fga * pa_zz * pb_xyyz

                     + 7.0 * pa_yy * fz * fx * pb_xyyz

                     + 7.0 * fx * pa_zz * fz * pb_xyyz

                     + 16.0 * pa_yyzz * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xyzz_s_0(double fx,
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
                                     double pb_xyzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx * fx * pb_x

                     + 0.375 * fx * fx * fx * pb_xy

                     + 0.75 * pa_yy * fx * fx * pb_xy

                     + 0.5 * pa_yzz * fx * fx * pb_x

                     + 2.0 * pa_yz * fx * fx * pb_xz

                     + 0.5 * pa_y * fx * fx * pb_xzz

                     + 0.25 * fx * fx * pa_zz * pb_xy

                     + fx * fx * pa_z * pb_xyz

                     + 0.5 * pa_yyzz * pb_xy * fx

                     + 2.0 * pa_yyz * fx * pb_xyz

                     + pa_yzz * fx * pb_xzz

                     + 0.25 * fx * fx * pb_xyzz

                     + 0.5 * pa_yy * fx * pb_xyzz

                     + 0.5 * fx * pa_zz * pb_xyzz

                     + pa_yyzz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xyzz_r_0(double fga,
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
                                     double pb_x,
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_xyzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double r_0_0)
    {
        return r_0_0 * (7.5 * pa_y * fx * fx * fx * fz * pb_x

                     - 0.5 * pa_y * fx * fx * pb_x * fz * fgb

                     - 0.5 * pa_y * fx * fx * fz * fga * pb_x

                     - fz * fga * fx * fx * pb_xy

                     - pa_yzz * fx * pb_x * fz * fgb

                     + 3.75 * fx * fx * fx * fz * pb_xy

                     + 9.0 * pa_yy * fx * fx * fz * pb_xy

                     + 6.0 * pa_yzz * fx * fx * fz * pb_x

                     + 24.0 * pa_yz * fx * fx * fz * pb_xz

                     - 0.25 * fx * fx * pb_xy * fz * fgb

                     - 0.5 * pa_yy * fx * pb_xy * fz * fgb

                     - 0.5 * pa_yy * fz * fga * pb_xy * fx

                     - pa_y * fx * fz * fga * pb_xzz

                     - 0.5 * fx * pa_zz * pb_xy * fz * fgb

                     - 0.5 * fz * fga * pa_zz * pb_xy * fx

                     - 2.0 * fz * fga * pa_z * fx * pb_xyz

                     - pa_yyzz * pb_xy * fz * fgb

                     + 6.0 * pa_y * fx * fx * fz * pb_xzz

                     + 3.0 * fx * fx * pa_zz * fz * pb_xy

                     + 12.0 * fx * fx * pa_z * fz * pb_xyz

                     + 7.0 * pa_yyzz * fz * pb_xy * fx

                     + 28.0 * pa_yyz * fz * fx * pb_xyz

                     + 14.0 * pa_yzz * fx * fz * pb_xzz

                     - fx * fz * fga * pb_xyzz

                     - pa_yy * fz * fga * pb_xyzz

                     + 3.0 * fx * fx * fz * pb_xyzz

                     - fz * fga * pa_zz * pb_xyzz

                     + 7.0 * pa_yy * fz * fx * pb_xyzz

                     + 7.0 * fx * pa_zz * fz * pb_xyzz

                     + 16.0 * pa_yyzz * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xzzz_s_0(double fx,
                                     double pa_yy,
                                     double pa_yyz,
                                     double pa_yyzz,
                                     double pa_z,
                                     double pa_zz,
                                     double pb_x,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_xzzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * fx * pa_z * pb_x

                     + 1.125 * fx * fx * fx * pb_xz

                     + 1.5 * pa_yyz * fx * fx * pb_x

                     + 2.25 * pa_yy * fx * fx * pb_xz

                     + 0.75 * fx * fx * pa_zz * pb_xz

                     + 1.5 * fx * fx * pa_z * pb_xzz

                     + 1.5 * pa_yyzz * pb_xz * fx

                     + 3.0 * pa_yyz * fx * pb_xzz

                     + 0.25 * fx * fx * pb_xzzz

                     + 0.5 * pa_yy * fx * pb_xzzz

                     + 0.5 * fx * pa_zz * pb_xzzz

                     + pa_yyzz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xzzz_r_0(double fga,
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
                                     double pb_xzzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_z * pb_x * fz * fgb

                     - 1.5 * fz * fga * pa_z * fx * fx * pb_x

                     - 3.0 * fz * fga * fx * fx * pb_xz

                     - 3.0 * pa_yyz * fx * pb_x * fz * fgb

                     + 7.5 * fx * fx * fx * pa_z * fz * pb_x

                     + 11.25 * fx * fx * fx * fz * pb_xz

                     + 18.0 * pa_yyz * fz * fx * fx * pb_x

                     + 27.0 * pa_yy * fx * fx * fz * pb_xz

                     - 0.75 * fx * fx * pb_xz * fz * fgb

                     - 1.5 * pa_yy * fx * pb_xz * fz * fgb

                     - 1.5 * pa_yy * fz * fga * pb_xz * fx

                     - 1.5 * fx * pa_zz * pb_xz * fz * fgb

                     - 1.5 * fz * fga * pa_zz * pb_xz * fx

                     - 3.0 * fz * fga * pa_z * fx * pb_xzz

                     - 3.0 * pa_yyzz * pb_xz * fz * fgb

                     + 9.0 * fx * fx * pa_zz * fz * pb_xz

                     + 18.0 * fx * fx * pa_z * fz * pb_xzz

                     + 21.0 * pa_yyzz * fz * pb_xz * fx

                     + 42.0 * pa_yyz * fz * fx * pb_xzz

                     - fx * fz * fga * pb_xzzz

                     - pa_yy * fz * fga * pb_xzzz

                     + 3.0 * fx * fx * fz * pb_xzzz

                     - fz * fga * pa_zz * pb_xzzz

                     + 7.0 * pa_yy * fz * fx * pb_xzzz

                     + 7.0 * fx * pa_zz * fz * pb_xzzz

                     + 16.0 * pa_yyzz * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_yyyy_s_0(double fx,
                                     double pa_y,
                                     double pa_yy,
                                     double pa_yyzz,
                                     double pa_yzz,
                                     double pa_zz,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double pb_yyyy,
                                     double s_0_0)
    {
        return s_0_0 * (0.9375 * fx * fx * fx * fx

                     + 1.875 * fx * fx * fx * pa_zz

                     + 0.375 * pa_yy * fx * fx * fx

                     + 3.0 * pa_y * fx * fx * fx * pb_y

                     + 2.25 * fx * fx * fx * pb_yy

                     + 0.75 * pa_yyzz * fx * fx

                     + 6.0 * pa_yzz * fx * fx * pb_y

                     + 4.5 * fx * fx * pa_zz * pb_yy

                     + 1.5 * pa_yy * fx * fx * pb_yy

                     + 2.0 * pa_y * fx * fx * pb_yyy

                     + 3.0 * pa_yyzz * pb_yy * fx

                     + 4.0 * pa_yzz * fx * pb_yyy

                     + 0.25 * fx * fx * pb_yyyy

                     + 0.5 * pa_yy * fx * pb_yyyy

                     + 0.5 * fx * pa_zz * pb_yyyy

                     + pa_yyzz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_yyyy_r_0(double fga,
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
                                     double pb_yyyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * fx * fz * fgb

                     - 2.25 * fx * fx * fx * fz * fga

                     - 4.5 * fx * fx * pa_zz * fz * fgb

                     + 7.5 * fx * fx * fx * fx * fz

                     + 18.75 * fx * fx * fx * pa_zz * fz

                     - 1.5 * pa_yy * fx * fx * fz * fgb

                     - 0.75 * pa_yy * fz * fga * fx * fx

                     - 6.0 * pa_y * fx * fx * pb_y * fz * fgb

                     - 6.0 * pa_y * fx * fx * fz * fga * pb_y

                     - 6.0 * fx * fx * fz * fga * pb_yy

                     - 0.75 * fz * fga * pa_zz * fx * fx

                     - 3.0 * pa_yyzz * fx * fz * fgb

                     + 3.75 * pa_yy * fz * fx * fx * fx

                     - 12.0 * pa_yzz * fx * pb_y * fz * fgb

                     + 30.0 * pa_y * fx * fx * fx * fz * pb_y

                     + 22.5 * fx * fx * fx * fz * pb_yy

                     + 9.0 * pa_yyzz * fz * fx * fx

                     + 72.0 * pa_yzz * fx * fx * fz * pb_y

                     + 54.0 * fx * fx * pa_zz * fz * pb_yy

                     - 1.5 * fx * fx * pb_yy * fz * fgb

                     - 3.0 * pa_yy * fx * pb_yy * fz * fgb

                     - 3.0 * pa_yy * fz * fga * pb_yy * fx

                     - 4.0 * pa_y * fx * fz * fga * pb_yyy

                     - 3.0 * fx * pa_zz * pb_yy * fz * fgb

                     - 3.0 * fz * fga * pa_zz * pb_yy * fx

                     - 6.0 * pa_yyzz * pb_yy * fz * fgb

                     + 18.0 * pa_yy * fz * fx * fx * pb_yy

                     + 24.0 * pa_y * fx * fx * fz * pb_yyy

                     + 42.0 * pa_yyzz * fz * pb_yy * fx

                     + 56.0 * pa_yzz * fx * fz * pb_yyy

                     - fx * fz * fga * pb_yyyy

                     - pa_yy * fz * fga * pb_yyyy

                     + 3.0 * fx * fx * fz * pb_yyyy

                     - fz * fga * pa_zz * pb_yyyy

                     + 7.0 * pa_yy * fz * fx * pb_yyyy

                     + 7.0 * fx * pa_zz * fz * pb_yyyy

                     + 16.0 * pa_yyzz * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_yyyz_s_0(double fx,
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
                                     double pb_yyy,
                                     double pb_yyyz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (1.5 * pa_yz * fx * fx * fx

                     + 2.25 * fx * fx * fx * pa_z * pb_y

                     + 0.75 * pa_y * fx * fx * fx * pb_z

                     + 1.125 * fx * fx * fx * pb_yz

                     + 1.5 * pa_yyz * fx * fx * pb_y

                     + 1.5 * pa_yzz * fx * fx * pb_z

                     + 3.0 * pa_yz * fx * fx * pb_yy

                     + 2.25 * fx * fx * pa_zz * pb_yz

                     + 0.75 * pa_yy * fx * fx * pb_yz

                     + 1.5 * pa_y * fx * fx * pb_yyz

                     + 0.5 * fx * fx * pa_z * pb_yyy

                     + 1.5 * pa_yyzz * pb_yz * fx

                     + pa_yyz * fx * pb_yyy

                     + 3.0 * pa_yzz * fx * pb_yyz

                     + 0.25 * fx * fx * pb_yyyz

                     + 0.5 * pa_yy * fx * pb_yyyz

                     + 0.5 * fx * pa_zz * pb_yyyz

                     + pa_yyzz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_yyyz_r_0(double fga,
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
                                     double pb_yyy,
                                     double pb_yyyz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_yz * fx * fx * fz * fgb

                     + 15.0 * pa_yz * fx * fx * fx * fz

                     + 22.5 * fx * fx * fx * pa_z * fz * pb_y

                     - 1.5 * pa_y * fx * fx * fz * fgb * pb_z

                     - 1.5 * pa_y * fx * fx * fz * fga * pb_z

                     - 1.5 * fx * fx * pa_z * pb_y * fz * fgb

                     - 3.0 * fx * fx * fz * fga * pb_yz

                     - 1.5 * fz * fga * pa_z * fx * fx * pb_y

                     - 3.0 * pa_yyz * fx * pb_y * fz * fgb

                     - 3.0 * pa_yzz * fx * fz * fgb * pb_z

                     + 7.5 * pa_y * fx * fx * fx * fz * pb_z

                     + 11.25 * fx * fx * fx * fz * pb_yz

                     + 18.0 * pa_yyz * fz * fx * fx * pb_y

                     + 18.0 * pa_yzz * fx * fx * fz * pb_z

                     + 36.0 * pa_yz * fx * fx * fz * pb_yy

                     + 27.0 * fx * fx * pa_zz * fz * pb_yz

                     - 0.75 * fx * fx * pb_yz * fz * fgb

                     - 1.5 * pa_yy * fx * pb_yz * fz * fgb

                     - 1.5 * pa_yy * fz * fga * pb_yz * fx

                     - 3.0 * pa_y * fx * fz * fga * pb_yyz

                     - 1.5 * fx * pa_zz * pb_yz * fz * fgb

                     - 1.5 * fz * fga * pa_zz * pb_yz * fx

                     - fz * fga * pa_z * fx * pb_yyy

                     - 3.0 * pa_yyzz * pb_yz * fz * fgb

                     + 9.0 * pa_yy * fz * fx * fx * pb_yz

                     + 18.0 * pa_y * fx * fx * fz * pb_yyz

                     + 6.0 * fx * fx * pa_z * fz * pb_yyy

                     + 21.0 * pa_yyzz * fz * pb_yz * fx

                     + 14.0 * pa_yyz * fz * fx * pb_yyy

                     + 42.0 * pa_yzz * fx * fz * pb_yyz

                     - fx * fz * fga * pb_yyyz

                     - pa_yy * fz * fga * pb_yyyz

                     + 3.0 * fx * fx * fz * pb_yyyz

                     - fz * fga * pa_zz * pb_yyyz

                     + 7.0 * pa_yy * fz * fx * pb_yyyz

                     + 7.0 * fx * pa_zz * fz * pb_yyyz

                     + 16.0 * pa_yyzz * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_yyzz_s_0(double fx,
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
                                     double pb_yyzz,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.5625 * fx * fx * fx * fx

                     + 0.375 * pa_yy * fx * fx * fx

                     + 1.5 * pa_y * fx * fx * fx * pb_y

                     + 0.375 * fx * fx * fx * pa_zz

                     + 1.5 * fx * fx * fx * pa_z * pb_z

                     + 0.375 * fx * fx * fx * pb_yy

                     + 0.375 * fx * fx * fx * pb_zz

                     + 0.25 * pa_yyzz * fx * fx

                     + pa_yyz * fx * fx * pb_z

                     + 0.75 * pa_yy * fx * fx * pb_yy

                     + pa_yzz * fx * fx * pb_y

                     + 4.0 * pa_yz * fx * fx * pb_yz

                     + 0.75 * fx * fx * pa_zz * pb_zz

                     + 0.25 * pa_yy * fx * fx * pb_zz

                     + pa_y * fx * fx * pb_yzz

                     + 0.25 * fx * fx * pa_zz * pb_yy

                     + fx * fx * pa_z * pb_yyz

                     + 0.5 * pa_yyzz * pb_yy * fx

                     + 0.5 * pa_yyzz * fx * pb_zz

                     + 2.0 * pa_yyz * fx * pb_yyz

                     + 2.0 * pa_yzz * fx * pb_yzz

                     + 0.25 * fx * fx * pb_yyzz

                     + 0.5 * pa_yy * fx * pb_yyzz

                     + 0.5 * fx * pa_zz * pb_yyzz

                     + pa_yyzz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_yyzz_r_0(double fga,
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
                                     double pb_yyzz,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (4.5 * fx * fx * fx * fx * fz

                     - 0.75 * fx * fx * fx * fz * fgb

                     - 0.75 * fx * fx * fx * fz * fga

                     - pa_yy * fx * fx * fz * fgb

                     - fx * fx * pa_zz * fz * fgb

                     + 3.75 * pa_yy * fx * fx * fx * fz

                     + 15.0 * pa_y * fx * fx * fx * fz * pb_y

                     + 3.75 * fx * fx * fx * pa_zz * fz

                     + 15.0 * fx * fx * fx * pa_z * fz * pb_z

                     - 0.25 * pa_yy * fz * fga * fx * fx

                     - pa_y * fx * fx * pb_y * fz * fgb

                     - pa_y * fx * fx * fz * fga * pb_y

                     - fx * fx * pa_z * fz * fgb * pb_z

                     - fx * fx * fz * fga * pb_zz

                     - 0.25 * fz * fga * pa_zz * fx * fx

                     - fz * fga * pa_z * fx * fx * pb_z

                     - fz * fga * fx * fx * pb_yy

                     - pa_yyzz * fx * fz * fgb

                     - 2.0 * pa_yyz * fx * fz * fgb * pb_z

                     - 2.0 * pa_yzz * fx * pb_y * fz * fgb

                     + 3.75 * fx * fx * fx * fz * pb_yy

                     + 3.75 * fx * fx * fx * fz * pb_zz

                     + 3.0 * pa_yyzz * fz * fx * fx

                     + 12.0 * pa_yyz * fz * fx * fx * pb_z

                     + 9.0 * pa_yy * fx * fx * fz * pb_yy

                     + 12.0 * pa_yzz * fx * fx * fz * pb_y

                     + 48.0 * pa_yz * fx * fx * fz * pb_yz

                     + 9.0 * fx * fx * pa_zz * fz * pb_zz

                     - 0.25 * fx * fx * pb_yy * fz * fgb

                     - 0.25 * fx * fx * fz * fgb * pb_zz

                     - 0.5 * pa_yy * fx * pb_yy * fz * fgb

                     - 0.5 * pa_yy * fx * fz * fgb * pb_zz

                     - 0.5 * pa_yy * fz * fga * pb_yy * fx

                     - 0.5 * pa_yy * fz * fga * fx * pb_zz

                     - 2.0 * pa_y * fx * fz * fga * pb_yzz

                     - 0.5 * fx * pa_zz * pb_yy * fz * fgb

                     - 0.5 * fx * pa_zz * fz * fgb * pb_zz

                     - 0.5 * fz * fga * pa_zz * pb_yy * fx

                     - 0.5 * fz * fga * pa_zz * fx * pb_zz

                     - 2.0 * fz * fga * pa_z * fx * pb_yyz

                     - pa_yyzz * pb_yy * fz * fgb

                     - pa_yyzz * fz * fgb * pb_zz

                     + 3.0 * pa_yy * fz * fx * fx * pb_zz

                     + 12.0 * pa_y * fx * fx * fz * pb_yzz

                     + 3.0 * fx * fx * pa_zz * fz * pb_yy

                     + 12.0 * fx * fx * pa_z * fz * pb_yyz

                     + 7.0 * pa_yyzz * fz * pb_yy * fx

                     + 7.0 * pa_yyzz * fz * fx * pb_zz

                     + 28.0 * pa_yyz * fz * fx * pb_yyz

                     + 28.0 * pa_yzz * fx * fz * pb_yzz

                     - fx * fz * fga * pb_yyzz

                     - pa_yy * fz * fga * pb_yyzz

                     + 3.0 * fx * fx * fz * pb_yyzz

                     - fz * fga * pa_zz * pb_yyzz

                     + 7.0 * pa_yy * fz * fx * pb_yyzz

                     + 7.0 * fx * pa_zz * fz * pb_yyzz

                     + 16.0 * pa_yyzz * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_yzzz_s_0(double fx,
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
                                     double pb_yzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.5 * pa_yz * fx * fx * fx

                     + 2.25 * pa_y * fx * fx * fx * pb_z

                     + 0.75 * fx * fx * fx * pa_z * pb_y

                     + 1.125 * fx * fx * fx * pb_yz

                     + 1.5 * pa_yyz * fx * fx * pb_y

                     + 2.25 * pa_yy * fx * fx * pb_yz

                     + 1.5 * pa_yzz * fx * fx * pb_z

                     + 3.0 * pa_yz * fx * fx * pb_zz

                     + 0.5 * pa_y * fx * fx * pb_zzz

                     + 0.75 * fx * fx * pa_zz * pb_yz

                     + 1.5 * fx * fx * pa_z * pb_yzz

                     + 1.5 * pa_yyzz * pb_yz * fx

                     + 3.0 * pa_yyz * fx * pb_yzz

                     + pa_yzz * fx * pb_zzz

                     + 0.25 * fx * fx * pb_yzzz

                     + 0.5 * pa_yy * fx * pb_yzzz

                     + 0.5 * fx * pa_zz * pb_yzzz

                     + pa_yyzz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_yzzz_r_0(double fga,
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
                                     double pb_yzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_yz * fx * fx * fz * fgb

                     + 15.0 * pa_yz * fx * fx * fx * fz

                     + 22.5 * pa_y * fx * fx * fx * fz * pb_z

                     - 1.5 * pa_y * fx * fx * pb_z * fz * fgb

                     - 1.5 * pa_y * fx * fx * fz * fga * pb_z

                     - 1.5 * fx * fx * pa_z * pb_y * fz * fgb

                     - 1.5 * fz * fga * pa_z * fx * fx * pb_y

                     - 3.0 * fz * fga * fx * fx * pb_yz

                     - 3.0 * pa_yyz * fx * pb_y * fz * fgb

                     - 3.0 * pa_yzz * fx * pb_z * fz * fgb

                     + 7.5 * fx * fx * fx * pa_z * fz * pb_y

                     + 11.25 * fx * fx * fx * fz * pb_yz

                     + 18.0 * pa_yyz * fz * fx * fx * pb_y

                     + 27.0 * pa_yy * fx * fx * fz * pb_yz

                     + 18.0 * pa_yzz * fx * fx * fz * pb_z

                     + 36.0 * pa_yz * fx * fx * fz * pb_zz

                     - 0.75 * fx * fx * pb_yz * fz * fgb

                     - 1.5 * pa_yy * fx * pb_yz * fz * fgb

                     - 1.5 * pa_yy * fz * fga * pb_yz * fx

                     - pa_y * fx * fz * fga * pb_zzz

                     - 1.5 * fx * pa_zz * pb_yz * fz * fgb

                     - 1.5 * fz * fga * pa_zz * pb_yz * fx

                     - 3.0 * fz * fga * pa_z * fx * pb_yzz

                     - 3.0 * pa_yyzz * pb_yz * fz * fgb

                     + 6.0 * pa_y * fx * fx * fz * pb_zzz

                     + 9.0 * fx * fx * pa_zz * fz * pb_yz

                     + 18.0 * fx * fx * pa_z * fz * pb_yzz

                     + 21.0 * pa_yyzz * fz * pb_yz * fx

                     + 42.0 * pa_yyz * fz * fx * pb_yzz

                     + 14.0 * pa_yzz * fx * fz * pb_zzz

                     - fx * fz * fga * pb_yzzz

                     - pa_yy * fz * fga * pb_yzzz

                     + 3.0 * fx * fx * fz * pb_yzzz

                     - fz * fga * pa_zz * pb_yzzz

                     + 7.0 * pa_yy * fz * fx * pb_yzzz

                     + 7.0 * fx * pa_zz * fz * pb_yzzz

                     + 16.0 * pa_yyzz * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_zzzz_s_0(double fx,
                                     double pa_yy,
                                     double pa_yyz,
                                     double pa_yyzz,
                                     double pa_z,
                                     double pa_zz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double pb_zzzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.9375 * fx * fx * fx * fx

                     + 1.875 * pa_yy * fx * fx * fx

                     + 0.375 * fx * fx * fx * pa_zz

                     + 3.0 * fx * fx * fx * pa_z * pb_z

                     + 2.25 * fx * fx * fx * pb_zz

                     + 0.75 * pa_yyzz * fx * fx

                     + 6.0 * pa_yyz * fx * fx * pb_z

                     + 4.5 * pa_yy * fx * fx * pb_zz

                     + 1.5 * fx * fx * pa_zz * pb_zz

                     + 2.0 * fx * fx * pa_z * pb_zzz

                     + 3.0 * pa_yyzz * pb_zz * fx

                     + 4.0 * pa_yyz * fx * pb_zzz

                     + 0.25 * fx * fx * pb_zzzz

                     + 0.5 * pa_yy * fx * pb_zzzz

                     + 0.5 * fx * pa_zz * pb_zzzz

                     + pa_yyzz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_zzzz_r_0(double fga,
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
                                     double pb_zzzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * fx * fz * fgb

                     - 2.25 * fz * fga * fx * fx * fx

                     - 4.5 * pa_yy * fx * fx * fz * fgb

                     + 7.5 * fx * fx * fx * fx * fz

                     + 18.75 * pa_yy * fx * fx * fx * fz

                     - 0.75 * pa_yy * fz * fga * fx * fx

                     - 1.5 * fx * fx * pa_zz * fz * fgb

                     - 6.0 * fx * fx * pa_z * pb_z * fz * fgb

                     - 0.75 * fz * fga * pa_zz * fx * fx

                     - 6.0 * fz * fga * pa_z * fx * fx * pb_z

                     - 6.0 * fz * fga * fx * fx * pb_zz

                     - 3.0 * pa_yyzz * fx * fz * fgb

                     - 12.0 * pa_yyz * fx * pb_z * fz * fgb

                     + 3.75 * fx * fx * fx * pa_zz * fz

                     + 30.0 * fx * fx * fx * pa_z * fz * pb_z

                     + 22.5 * fx * fx * fx * fz * pb_zz

                     + 9.0 * pa_yyzz * fz * fx * fx

                     + 72.0 * pa_yyz * fz * fx * fx * pb_z

                     + 54.0 * pa_yy * fx * fx * fz * pb_zz

                     - 1.5 * fx * fx * pb_zz * fz * fgb

                     - 3.0 * pa_yy * fx * pb_zz * fz * fgb

                     - 3.0 * pa_yy * fz * fga * pb_zz * fx

                     - 3.0 * fx * pa_zz * pb_zz * fz * fgb

                     - 3.0 * fz * fga * pa_zz * pb_zz * fx

                     - 4.0 * fz * fga * pa_z * fx * pb_zzz

                     - 6.0 * pa_yyzz * pb_zz * fz * fgb

                     + 18.0 * fx * fx * pa_zz * fz * pb_zz

                     + 24.0 * fx * fx * pa_z * fz * pb_zzz

                     + 42.0 * pa_yyzz * fz * pb_zz * fx

                     + 56.0 * pa_yyz * fz * fx * pb_zzz

                     - fx * fz * fga * pb_zzzz

                     - pa_yy * fz * fga * pb_zzzz

                     + 3.0 * fx * fx * fz * pb_zzzz

                     - fz * fga * pa_zz * pb_zzzz

                     + 7.0 * pa_yy * fz * fx * pb_zzzz

                     + 7.0 * fx * pa_zz * fz * pb_zzzz

                     + 16.0 * pa_yyzz * fz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xxxx_s_0(double fx,
                                     double pa_yz,
                                     double pa_yzzz,
                                     double pb_xx,
                                     double pb_xxxx,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_yz * fx * fx * fx

                     + 0.75 * pa_yzzz * fx * fx

                     + 4.5 * pa_yz * fx * fx * pb_xx

                     + 3.0 * pa_yzzz * pb_xx * fx

                     + 1.5 * pa_yz * fx * pb_xxxx

                     + pa_yzzz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xxxx_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_yz,
                                     double pa_yzzz,
                                     double pb_xx,
                                     double pb_xxxx,
                                     double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_yz * fx * fx * fz * fgb

                     - 2.25 * pa_yz * fz * fga * fx * fx

                     - 3.0 * pa_yzzz * fx * fz * fgb

                     + 11.25 * pa_yz * fz * fx * fx * fx

                     + 9.0 * pa_yzzz * fz * fx * fx

                     - 9.0 * pa_yz * fx * pb_xx * fz * fgb

                     - 9.0 * pa_yz * fz * fga * pb_xx * fx

                     - 6.0 * pa_yzzz * pb_xx * fz * fgb

                     + 54.0 * pa_yz * fz * fx * fx * pb_xx

                     + 42.0 * pa_yzzz * fz * pb_xx * fx

                     - 3.0 * pa_yz * fz * fga * pb_xxxx

                     + 21.0 * pa_yz * fz * fx * pb_xxxx

                     + 16.0 * pa_yzzz * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xxxy_s_0(double fx,
                                     double pa_yz,
                                     double pa_yzzz,
                                     double pa_z,
                                     double pa_zzz,
                                     double pb_x,
                                     double pb_xxx,
                                     double pb_xxxy,
                                     double pb_xy,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pa_z * pb_x

                     + 0.75 * fx * fx * pa_zzz * pb_x

                     + 2.25 * pa_yz * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_z * pb_xxx

                     + 1.5 * pa_yzzz * pb_xy * fx

                     + 0.5 * fx * pa_zzz * pb_xxx

                     + 1.5 * pa_yz * fx * pb_xxxy

                     + pa_yzzz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xxxy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_yz,
                                     double pa_yzzz,
                                     double pa_z,
                                     double pa_zzz,
                                     double pb_x,
                                     double pb_xxx,
                                     double pb_xxxy,
                                     double pb_xy,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pa_z * pb_x * fz * fgb

                     - 2.25 * fx * fx * pa_z * fz * fga * pb_x

                     - 1.5 * fx * pa_zzz * pb_x * fz * fgb

                     + 11.25 * fx * fx * fx * pa_z * fz * pb_x

                     + 9.0 * fx * fx * pa_zzz * fz * pb_x

                     - 4.5 * pa_yz * fx * pb_xy * fz * fgb

                     - 4.5 * pa_yz * fz * fga * pb_xy * fx

                     - 1.5 * fx * pa_z * fz * fga * pb_xxx

                     - 3.0 * pa_yzzz * pb_xy * fz * fgb

                     + 27.0 * pa_yz * fz * fx * fx * pb_xy

                     + 9.0 * fx * fx * pa_z * fz * pb_xxx

                     + 21.0 * pa_yzzz * fz * pb_xy * fx

                     + 7.0 * fx * pa_zzz * fz * pb_xxx

                     - 3.0 * pa_yz * fz * fga * pb_xxxy

                     + 21.0 * pa_yz * fz * fx * pb_xxxy

                     + 16.0 * pa_yzzz * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xxxz_s_0(double fx,
                                     double pa_y,
                                     double pa_yz,
                                     double pa_yzz,
                                     double pa_yzzz,
                                     double pb_x,
                                     double pb_xxx,
                                     double pb_xxxz,
                                     double pb_xz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_y * fx * fx * fx * pb_x

                     + 2.25 * pa_yzz * fx * fx * pb_x

                     + 2.25 * pa_yz * fx * fx * pb_xz

                     + 0.75 * pa_y * fx * fx * pb_xxx

                     + 1.5 * pa_yzzz * pb_xz * fx

                     + 1.5 * pa_yzz * fx * pb_xxx

                     + 1.5 * pa_yz * fx * pb_xxxz

                     + pa_yzzz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xxxz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_y,
                                     double pa_yz,
                                     double pa_yzz,
                                     double pa_yzzz,
                                     double pb_x,
                                     double pb_xxx,
                                     double pb_xxxz,
                                     double pb_xz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * pa_y * fx * fx * pb_x * fz * fgb

                     - 2.25 * pa_y * fx * fx * fz * fga * pb_x

                     - 4.5 * pa_yzz * fx * pb_x * fz * fgb

                     + 11.25 * pa_y * fx * fx * fx * fz * pb_x

                     + 27.0 * pa_yzz * fz * fx * fx * pb_x

                     - 4.5 * pa_yz * fx * pb_xz * fz * fgb

                     - 4.5 * pa_yz * fz * fga * pb_xz * fx

                     - 1.5 * pa_y * fx * fz * fga * pb_xxx

                     - 3.0 * pa_yzzz * pb_xz * fz * fgb

                     + 27.0 * pa_yz * fz * fx * fx * pb_xz

                     + 9.0 * pa_y * fx * fx * fz * pb_xxx

                     + 21.0 * pa_yzzz * fz * pb_xz * fx

                     + 21.0 * pa_yzz * fz * fx * pb_xxx

                     - 3.0 * pa_yz * fz * fga * pb_xxxz

                     + 21.0 * pa_yz * fz * fx * pb_xxxz

                     + 16.0 * pa_yzzz * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xxyy_s_0(double fx,
                                     double pa_yz,
                                     double pa_yzzz,
                                     double pa_z,
                                     double pa_zzz,
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyy,
                                     double pb_y,
                                     double pb_yy,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * pa_yz * fx * fx * fx

                     + 0.75 * fx * fx * fx * pa_z * pb_y

                     + 0.25 * pa_yzzz * fx * fx

                     + 0.5 * fx * fx * pa_zzz * pb_y

                     + 0.75 * pa_yz * fx * fx * pb_xx

                     + 0.75 * pa_yz * fx * fx * pb_yy

                     + 1.5 * fx * fx * pa_z * pb_xxy

                     + 0.5 * pa_yzzz * pb_xx * fx

                     + 0.5 * pa_yzzz * fx * pb_yy

                     + fx * pa_zzz * pb_xxy

                     + 1.5 * pa_yz * fx * pb_xxyy

                     + pa_yzzz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xxyy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_yz,
                                     double pa_yzzz,
                                     double pa_z,
                                     double pa_zzz,
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyy,
                                     double pb_y,
                                     double pb_yy,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_yz * fx * fx * fz * fgb

                     - 0.75 * pa_yz * fz * fga * fx * fx

                     - 1.5 * fx * fx * pa_z * fz * fgb * pb_y

                     - 1.5 * fx * fx * pa_z * fz * fga * pb_y

                     - pa_yzzz * fx * fz * fgb

                     + 3.75 * pa_yz * fz * fx * fx * fx

                     - fx * pa_zzz * fz * fgb * pb_y

                     + 7.5 * fx * fx * fx * pa_z * fz * pb_y

                     + 3.0 * pa_yzzz * fz * fx * fx

                     + 6.0 * fx * fx * pa_zzz * fz * pb_y

                     - 1.5 * pa_yz * fx * pb_xx * fz * fgb

                     - 1.5 * pa_yz * fx * fz * fgb * pb_yy

                     - 1.5 * pa_yz * fz * fga * pb_xx * fx

                     - 1.5 * pa_yz * fz * fga * fx * pb_yy

                     - 3.0 * fx * pa_z * fz * fga * pb_xxy

                     - pa_yzzz * pb_xx * fz * fgb

                     - pa_yzzz * fz * fgb * pb_yy

                     + 9.0 * pa_yz * fz * fx * fx * pb_xx

                     + 9.0 * pa_yz * fz * fx * fx * pb_yy

                     + 18.0 * fx * fx * pa_z * fz * pb_xxy

                     + 7.0 * pa_yzzz * fz * pb_xx * fx

                     + 7.0 * pa_yzzz * fz * fx * pb_yy

                     + 14.0 * fx * pa_zzz * fz * pb_xxy

                     - 3.0 * pa_yz * fz * fga * pb_xxyy

                     + 21.0 * pa_yz * fz * fx * pb_xxyy

                     + 16.0 * pa_yzzz * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xxyz_s_0(double fx,
                                     double pa_y,
                                     double pa_yz,
                                     double pa_yzz,
                                     double pa_yzzz,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyz,
                                     double pb_xxz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.375 * fx * fx * fx * pa_zz

                     + 0.375 * pa_y * fx * fx * fx * pb_y

                     + 0.375 * fx * fx * fx * pa_z * pb_z

                     + 0.375 * fx * fx * fx * pb_xx

                     + 0.75 * pa_yzz * fx * fx * pb_y

                     + 0.25 * fx * fx * pa_zzz * pb_z

                     + 0.75 * fx * fx * pa_zz * pb_xx

                     + 0.75 * pa_yz * fx * fx * pb_yz

                     + 0.75 * pa_y * fx * fx * pb_xxy

                     + 0.75 * fx * fx * pa_z * pb_xxz

                     + 0.5 * pa_yzzz * fx * pb_yz

                     + 1.5 * pa_yzz * fx * pb_xxy

                     + 0.5 * fx * pa_zzz * pb_xxz

                     + 1.5 * pa_yz * fx * pb_xxyz

                     + pa_yzzz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xxyz_r_0(double fga,
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
                                     double pb_xx,
                                     double pb_xxy,
                                     double pb_xxyz,
                                     double pb_xxz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.375 * fx * fx * fx * fz * fgb

                     - 0.375 * fx * fx * fx * fz * fga

                     - 0.75 * fx * fx * pa_zz * fz * fgb

                     + 1.5 * fx * fx * fx * fx * fz

                     + 3.75 * fx * fx * fx * pa_zz * fz

                     - 0.75 * pa_y * fx * fx * fz * fgb * pb_y

                     - 0.75 * pa_y * fx * fx * fz * fga * pb_y

                     - 0.75 * fx * fx * pa_z * fz * fgb * pb_z

                     - 0.75 * fx * fx * pa_z * fz * fga * pb_z

                     - 0.75 * fx * fx * fz * fga * pb_xx

                     - 1.5 * pa_yzz * fx * fz * fgb * pb_y

                     + 3.75 * pa_y * fx * fx * fx * fz * pb_y

                     - 0.5 * fx * pa_zzz * fz * fgb * pb_z

                     + 3.75 * fx * fx * fx * pa_z * fz * pb_z

                     + 3.75 * fx * fx * fx * fz * pb_xx

                     + 9.0 * pa_yzz * fz * fx * fx * pb_y

                     + 3.0 * fx * fx * pa_zzz * fz * pb_z

                     + 9.0 * fx * fx * pa_zz * fz * pb_xx

                     - 1.5 * pa_yz * fx * fz * fgb * pb_yz

                     - 1.5 * pa_yz * fz * fga * fx * pb_yz

                     - 1.5 * pa_y * fx * fz * fga * pb_xxy

                     - 1.5 * fx * pa_z * fz * fga * pb_xxz

                     - pa_yzzz * fz * fgb * pb_yz

                     + 9.0 * pa_yz * fz * fx * fx * pb_yz

                     + 9.0 * pa_y * fx * fx * fz * pb_xxy

                     + 9.0 * fx * fx * pa_z * fz * pb_xxz

                     + 7.0 * pa_yzzz * fz * fx * pb_yz

                     + 21.0 * pa_yzz * fz * fx * pb_xxy

                     + 7.0 * fx * pa_zzz * fz * pb_xxz

                     - 3.0 * pa_yz * fz * fga * pb_xxyz

                     + 21.0 * pa_yz * fz * fx * pb_xxyz

                     + 16.0 * pa_yzzz * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xxzz_s_0(double fx,
                                     double pa_y,
                                     double pa_yz,
                                     double pa_yzz,
                                     double pa_yzzz,
                                     double pb_xx,
                                     double pb_xxz,
                                     double pb_xxzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_yz * fx * fx * fx

                     + 0.75 * pa_y * fx * fx * fx * pb_z

                     + 0.25 * pa_yzzz * fx * fx

                     + 1.5 * pa_yzz * fx * fx * pb_z

                     + 2.25 * pa_yz * fx * fx * pb_xx

                     + 0.75 * pa_yz * fx * fx * pb_zz

                     + 1.5 * pa_y * fx * fx * pb_xxz

                     + 0.5 * pa_yzzz * pb_xx * fx

                     + 0.5 * pa_yzzz * fx * pb_zz

                     + 3.0 * pa_yzz * fx * pb_xxz

                     + 1.5 * pa_yz * fx * pb_xxzz

                     + pa_yzzz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xxzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_y,
                                     double pa_yz,
                                     double pa_yzz,
                                     double pa_yzzz,
                                     double pb_xx,
                                     double pb_xxz,
                                     double pb_xxzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_yz * fx * fx * fz * fgb

                     + 11.25 * pa_yz * fx * fx * fx * fz

                     - 0.75 * pa_yz * fz * fga * fx * fx

                     - 1.5 * pa_y * fx * fx * fz * fgb * pb_z

                     - 1.5 * pa_y * fx * fx * fz * fga * pb_z

                     - pa_yzzz * fx * fz * fgb

                     - 3.0 * pa_yzz * fx * fz * fgb * pb_z

                     + 7.5 * pa_y * fx * fx * fx * fz * pb_z

                     + 3.0 * pa_yzzz * fz * fx * fx

                     + 18.0 * pa_yzz * fz * fx * fx * pb_z

                     + 27.0 * pa_yz * fx * fx * fz * pb_xx

                     - 1.5 * pa_yz * fx * pb_xx * fz * fgb

                     - 1.5 * pa_yz * fx * fz * fgb * pb_zz

                     - 1.5 * pa_yz * fz * fga * pb_xx * fx

                     - 1.5 * pa_yz * fz * fga * fx * pb_zz

                     - 3.0 * pa_y * fx * fz * fga * pb_xxz

                     - pa_yzzz * pb_xx * fz * fgb

                     - pa_yzzz * fz * fgb * pb_zz

                     + 9.0 * pa_yz * fz * fx * fx * pb_zz

                     + 18.0 * pa_y * fx * fx * fz * pb_xxz

                     + 7.0 * pa_yzzz * fz * pb_xx * fx

                     + 7.0 * pa_yzzz * fz * fx * pb_zz

                     + 42.0 * pa_yzz * fz * fx * pb_xxz

                     - 3.0 * pa_yz * fz * fga * pb_xxzz

                     + 21.0 * pa_yz * fz * fx * pb_xxzz

                     + 16.0 * pa_yzzz * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xyyy_s_0(double fx,
                                     double pa_yz,
                                     double pa_yzzz,
                                     double pa_z,
                                     double pa_zzz,
                                     double pb_x,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_xyyy,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pa_z * pb_x

                     + 0.75 * fx * fx * pa_zzz * pb_x

                     + 2.25 * pa_yz * fx * fx * pb_xy

                     + 2.25 * fx * fx * pa_z * pb_xyy

                     + 1.5 * pa_yzzz * pb_xy * fx

                     + 1.5 * fx * pa_zzz * pb_xyy

                     + 1.5 * pa_yz * fx * pb_xyyy

                     + pa_yzzz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xyyy_r_0(double fga,
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
                                     double pb_xyyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pa_z * pb_x * fz * fgb

                     - 2.25 * fx * fx * pa_z * fz * fga * pb_x

                     - 1.5 * fx * pa_zzz * pb_x * fz * fgb

                     + 11.25 * fx * fx * fx * pa_z * fz * pb_x

                     + 9.0 * fx * fx * pa_zzz * fz * pb_x

                     - 4.5 * pa_yz * fx * pb_xy * fz * fgb

                     - 4.5 * pa_yz * fz * fga * pb_xy * fx

                     - 4.5 * fx * pa_z * fz * fga * pb_xyy

                     - 3.0 * pa_yzzz * pb_xy * fz * fgb

                     + 27.0 * pa_yz * fz * fx * fx * pb_xy

                     + 27.0 * fx * fx * pa_z * fz * pb_xyy

                     + 21.0 * pa_yzzz * fz * pb_xy * fx

                     + 21.0 * fx * pa_zzz * fz * pb_xyy

                     - 3.0 * pa_yz * fz * fga * pb_xyyy

                     + 21.0 * pa_yz * fz * fx * pb_xyyy

                     + 16.0 * pa_yzzz * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xyyz_s_0(double fx,
                                     double pa_y,
                                     double pa_yz,
                                     double pa_yzz,
                                     double pa_yzzz,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pb_x,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_xyyz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double s_0_0)
    {
        return s_0_0 * (0.375 * pa_y * fx * fx * fx * pb_x

                     + 0.75 * fx * fx * fx * pb_xy

                     + 0.75 * pa_yzz * fx * fx * pb_x

                     + 1.5 * fx * fx * pa_zz * pb_xy

                     + 0.75 * pa_yz * fx * fx * pb_xz

                     + 0.75 * pa_y * fx * fx * pb_xyy

                     + 1.5 * fx * fx * pa_z * pb_xyz

                     + 0.5 * pa_yzzz * pb_xz * fx

                     + 1.5 * pa_yzz * fx * pb_xyy

                     + fx * pa_zzz * pb_xyz

                     + 1.5 * pa_yz * fx * pb_xyyz

                     + pa_yzzz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xyyz_r_0(double fga,
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
                                     double pb_x,
                                     double pb_xy,
                                     double pb_xyy,
                                     double pb_xyyz,
                                     double pb_xyz,
                                     double pb_xz,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.75 * pa_y * fx * fx * pb_x * fz * fgb

                     - 0.75 * pa_y * fx * fx * fz * fga * pb_x

                     - 1.5 * fx * fx * fz * fga * pb_xy

                     - 1.5 * pa_yzz * fx * pb_x * fz * fgb

                     + 3.75 * pa_y * fx * fx * fx * fz * pb_x

                     + 7.5 * fx * fx * fx * fz * pb_xy

                     + 9.0 * pa_yzz * fz * fx * fx * pb_x

                     + 18.0 * fx * fx * pa_zz * fz * pb_xy

                     - 1.5 * pa_yz * fx * pb_xz * fz * fgb

                     - 1.5 * pa_yz * fz * fga * pb_xz * fx

                     - 1.5 * pa_y * fx * fz * fga * pb_xyy

                     - 3.0 * fx * pa_z * fz * fga * pb_xyz

                     - pa_yzzz * pb_xz * fz * fgb

                     + 9.0 * pa_yz * fz * fx * fx * pb_xz

                     + 9.0 * pa_y * fx * fx * fz * pb_xyy

                     + 18.0 * fx * fx * pa_z * fz * pb_xyz

                     + 7.0 * pa_yzzz * fz * pb_xz * fx

                     + 21.0 * pa_yzz * fz * fx * pb_xyy

                     + 14.0 * fx * pa_zzz * fz * pb_xyz

                     - 3.0 * pa_yz * fz * fga * pb_xyyz

                     + 21.0 * pa_yz * fz * fx * pb_xyyz

                     + 16.0 * pa_yzzz * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xyzz_s_0(double fx,
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
                                     double pb_xyzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pa_z * pb_x

                     + 0.75 * fx * fx * fx * pb_xz

                     + 2.25 * pa_yz * fx * fx * pb_xy

                     + 0.25 * fx * fx * pa_zzz * pb_x

                     + 1.5 * fx * fx * pa_zz * pb_xz

                     + 1.5 * pa_y * fx * fx * pb_xyz

                     + 0.75 * fx * fx * pa_z * pb_xzz

                     + 0.5 * pa_yzzz * pb_xy * fx

                     + 3.0 * pa_yzz * fx * pb_xyz

                     + 0.5 * fx * pa_zzz * pb_xzz

                     + 1.5 * pa_yz * fx * pb_xyzz

                     + pa_yzzz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xyzz_r_0(double fga,
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
                                     double pb_x,
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_xyzz,
                                     double pb_xz,
                                     double pb_xzz,
                                     double r_0_0)
    {
        return r_0_0 * (11.25 * fx * fx * fx * pa_z * fz * pb_x

                     - 0.75 * fx * fx * pa_z * pb_x * fz * fgb

                     - 0.75 * fx * fx * pa_z * fz * fga * pb_x

                     - 1.5 * fx * fx * fz * fga * pb_xz

                     - 0.5 * fx * pa_zzz * pb_x * fz * fgb

                     + 7.5 * fx * fx * fx * fz * pb_xz

                     + 27.0 * pa_yz * fx * fx * fz * pb_xy

                     + 3.0 * fx * fx * pa_zzz * fz * pb_x

                     + 18.0 * fx * fx * pa_zz * fz * pb_xz

                     - 1.5 * pa_yz * fx * pb_xy * fz * fgb

                     - 1.5 * pa_yz * fz * fga * pb_xy * fx

                     - 3.0 * pa_y * fx * fz * fga * pb_xyz

                     - 1.5 * fx * pa_z * fz * fga * pb_xzz

                     - pa_yzzz * pb_xy * fz * fgb

                     + 18.0 * pa_y * fx * fx * fz * pb_xyz

                     + 9.0 * fx * fx * pa_z * fz * pb_xzz

                     + 7.0 * pa_yzzz * fz * pb_xy * fx

                     + 42.0 * pa_yzz * fz * fx * pb_xyz

                     + 7.0 * fx * pa_zzz * fz * pb_xzz

                     - 3.0 * pa_yz * fz * fga * pb_xyzz

                     + 21.0 * pa_yz * fz * fx * pb_xyzz

                     + 16.0 * pa_yzzz * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xzzz_s_0(double fx,
                                     double pa_y,
                                     double pa_yz,
                                     double pa_yzz,
                                     double pa_yzzz,
                                     double pb_x,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_xzzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.875 * pa_y * fx * fx * fx * pb_x

                     + 2.25 * pa_yzz * fx * fx * pb_x

                     + 6.75 * pa_yz * fx * fx * pb_xz

                     + 2.25 * pa_y * fx * fx * pb_xzz

                     + 1.5 * pa_yzzz * pb_xz * fx

                     + 4.5 * pa_yzz * fx * pb_xzz

                     + 1.5 * pa_yz * fx * pb_xzzz

                     + pa_yzzz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xzzz_r_0(double fga,
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
                                     double pb_xzzz,
                                     double r_0_0)
    {
        return r_0_0 * (18.75 * pa_y * fx * fx * fx * fz * pb_x

                     - 2.25 * pa_y * fx * fx * pb_x * fz * fgb

                     - 2.25 * pa_y * fx * fx * fz * fga * pb_x

                     - 4.5 * pa_yzz * fx * pb_x * fz * fgb

                     + 27.0 * pa_yzz * fz * fx * fx * pb_x

                     + 81.0 * pa_yz * fx * fx * fz * pb_xz

                     - 4.5 * pa_yz * fx * pb_xz * fz * fgb

                     - 4.5 * pa_yz * fz * fga * pb_xz * fx

                     - 4.5 * pa_y * fx * fz * fga * pb_xzz

                     - 3.0 * pa_yzzz * pb_xz * fz * fgb

                     + 27.0 * pa_y * fx * fx * fz * pb_xzz

                     + 21.0 * pa_yzzz * fz * pb_xz * fx

                     + 63.0 * pa_yzz * fz * fx * pb_xzz

                     - 3.0 * pa_yz * fz * fga * pb_xzzz

                     + 21.0 * pa_yz * fz * fx * pb_xzzz

                     + 16.0 * pa_yzzz * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_yyyy_s_0(double fx,
                                     double pa_yz,
                                     double pa_yzzz,
                                     double pa_z,
                                     double pa_zzz,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double pb_yyyy,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_yz * fx * fx * fx

                     + 4.5 * fx * fx * fx * pa_z * pb_y

                     + 0.75 * pa_yzzz * fx * fx

                     + 3.0 * fx * fx * pa_zzz * pb_y

                     + 4.5 * pa_yz * fx * fx * pb_yy

                     + 3.0 * fx * fx * pa_z * pb_yyy

                     + 3.0 * pa_yzzz * pb_yy * fx

                     + 2.0 * fx * pa_zzz * pb_yyy

                     + 1.5 * pa_yz * fx * pb_yyyy

                     + pa_yzzz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_yyyy_r_0(double fga,
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
                                     double pb_yyyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_yz * fx * fx * fz * fgb

                     - 2.25 * pa_yz * fz * fga * fx * fx

                     - 9.0 * fx * fx * pa_z * pb_y * fz * fgb

                     - 9.0 * fx * fx * pa_z * fz * fga * pb_y

                     - 3.0 * pa_yzzz * fx * fz * fgb

                     + 11.25 * pa_yz * fz * fx * fx * fx

                     - 6.0 * fx * pa_zzz * pb_y * fz * fgb

                     + 45.0 * fx * fx * fx * pa_z * fz * pb_y

                     + 9.0 * pa_yzzz * fz * fx * fx

                     + 36.0 * fx * fx * pa_zzz * fz * pb_y

                     - 9.0 * pa_yz * fx * pb_yy * fz * fgb

                     - 9.0 * pa_yz * fz * fga * pb_yy * fx

                     - 6.0 * fx * pa_z * fz * fga * pb_yyy

                     - 6.0 * pa_yzzz * pb_yy * fz * fgb

                     + 54.0 * pa_yz * fz * fx * fx * pb_yy

                     + 36.0 * fx * fx * pa_z * fz * pb_yyy

                     + 42.0 * pa_yzzz * fz * pb_yy * fx

                     + 28.0 * fx * pa_zzz * fz * pb_yyy

                     - 3.0 * pa_yz * fz * fga * pb_yyyy

                     + 21.0 * pa_yz * fz * fx * pb_yyyy

                     + 16.0 * pa_yzzz * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_yyyz_s_0(double fx,
                                     double pa_y,
                                     double pa_yz,
                                     double pa_yzz,
                                     double pa_yzzz,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pb_y,
                                     double pb_yy,
                                     double pb_yyy,
                                     double pb_yyyz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double s_0_0)
    {
        return s_0_0 * (0.5625 * fx * fx * fx * fx

                     + 1.125 * fx * fx * fx * pa_zz

                     + 1.125 * pa_y * fx * fx * fx * pb_y

                     + 1.125 * fx * fx * fx * pa_z * pb_z

                     + 1.125 * fx * fx * fx * pb_yy

                     + 2.25 * pa_yzz * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_zzz * pb_z

                     + 2.25 * fx * fx * pa_zz * pb_yy

                     + 2.25 * pa_yz * fx * fx * pb_yz

                     + 0.75 * pa_y * fx * fx * pb_yyy

                     + 2.25 * fx * fx * pa_z * pb_yyz

                     + 1.5 * pa_yzzz * pb_yz * fx

                     + 1.5 * pa_yzz * fx * pb_yyy

                     + 1.5 * fx * pa_zzz * pb_yyz

                     + 1.5 * pa_yz * fx * pb_yyyz

                     + pa_yzzz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_yyyz_r_0(double fga,
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
                                     double pb_yyy,
                                     double pb_yyyz,
                                     double pb_yyz,
                                     double pb_yz,
                                     double pb_z,
                                     double r_0_0)
    {
        return r_0_0 * (- 1.125 * fx * fx * fx * fz * fgb

                     - 1.125 * fx * fx * fx * fz * fga

                     - 2.25 * fx * fx * pa_zz * fz * fgb

                     + 4.5 * fx * fx * fx * fx * fz

                     + 11.25 * fx * fx * fx * pa_zz * fz

                     - 2.25 * pa_y * fx * fx * pb_y * fz * fgb

                     - 2.25 * pa_y * fx * fx * fz * fga * pb_y

                     - 2.25 * fx * fx * pa_z * fz * fgb * pb_z

                     - 2.25 * fx * fx * pa_z * fz * fga * pb_z

                     - 2.25 * fx * fx * fz * fga * pb_yy

                     - 4.5 * pa_yzz * fx * pb_y * fz * fgb

                     + 11.25 * pa_y * fx * fx * fx * fz * pb_y

                     - 1.5 * fx * pa_zzz * fz * fgb * pb_z

                     + 11.25 * fx * fx * fx * pa_z * fz * pb_z

                     + 11.25 * fx * fx * fx * fz * pb_yy

                     + 27.0 * pa_yzz * fz * fx * fx * pb_y

                     + 9.0 * fx * fx * pa_zzz * fz * pb_z

                     + 27.0 * fx * fx * pa_zz * fz * pb_yy

                     - 4.5 * pa_yz * fx * pb_yz * fz * fgb

                     - 4.5 * pa_yz * fz * fga * pb_yz * fx

                     - 1.5 * pa_y * fx * fz * fga * pb_yyy

                     - 4.5 * fx * pa_z * fz * fga * pb_yyz

                     - 3.0 * pa_yzzz * pb_yz * fz * fgb

                     + 27.0 * pa_yz * fz * fx * fx * pb_yz

                     + 9.0 * pa_y * fx * fx * fz * pb_yyy

                     + 27.0 * fx * fx * pa_z * fz * pb_yyz

                     + 21.0 * pa_yzzz * fz * pb_yz * fx

                     + 21.0 * pa_yzz * fz * fx * pb_yyy

                     + 21.0 * fx * pa_zzz * fz * pb_yyz

                     - 3.0 * pa_yz * fz * fga * pb_yyyz

                     + 21.0 * pa_yz * fz * fx * pb_yyyz

                     + 16.0 * pa_yzzz * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_yyzz_s_0(double fx,
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
                                     double pb_yyzz,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * pa_yz * fx * fx * fx

                     + 2.25 * fx * fx * fx * pa_z * pb_y

                     + 0.75 * pa_y * fx * fx * fx * pb_z

                     + 1.5 * fx * fx * fx * pb_yz

                     + 0.25 * pa_yzzz * fx * fx

                     + 1.5 * pa_yzz * fx * fx * pb_z

                     + 2.25 * pa_yz * fx * fx * pb_yy

                     + 0.5 * fx * fx * pa_zzz * pb_y

                     + 3.0 * fx * fx * pa_zz * pb_yz

                     + 0.75 * pa_yz * fx * fx * pb_zz

                     + 1.5 * pa_y * fx * fx * pb_yyz

                     + 1.5 * fx * fx * pa_z * pb_yzz

                     + 0.5 * pa_yzzz * pb_yy * fx

                     + 0.5 * pa_yzzz * fx * pb_zz

                     + 3.0 * pa_yzz * fx * pb_yyz

                     + fx * pa_zzz * pb_yzz

                     + 1.5 * pa_yz * fx * pb_yyzz

                     + pa_yzzz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_yyzz_r_0(double fga,
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
                                     double pb_yyzz,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_yz * fx * fx * fz * fgb

                     + 11.25 * pa_yz * fx * fx * fx * fz

                     + 22.5 * fx * fx * fx * pa_z * fz * pb_y

                     - 0.75 * pa_yz * fz * fga * fx * fx

                     - 1.5 * pa_y * fx * fx * fz * fgb * pb_z

                     - 1.5 * pa_y * fx * fx * fz * fga * pb_z

                     - 1.5 * fx * fx * pa_z * pb_y * fz * fgb

                     - 1.5 * fx * fx * pa_z * fz * fga * pb_y

                     - 3.0 * fx * fx * fz * fga * pb_yz

                     - pa_yzzz * fx * fz * fgb

                     - 3.0 * pa_yzz * fx * fz * fgb * pb_z

                     + 7.5 * pa_y * fx * fx * fx * fz * pb_z

                     - fx * pa_zzz * pb_y * fz * fgb

                     + 15.0 * fx * fx * fx * fz * pb_yz

                     + 3.0 * pa_yzzz * fz * fx * fx

                     + 18.0 * pa_yzz * fz * fx * fx * pb_z

                     + 27.0 * pa_yz * fx * fx * fz * pb_yy

                     + 6.0 * fx * fx * pa_zzz * fz * pb_y

                     + 36.0 * fx * fx * pa_zz * fz * pb_yz

                     - 1.5 * pa_yz * fx * pb_yy * fz * fgb

                     - 1.5 * pa_yz * fx * fz * fgb * pb_zz

                     - 1.5 * pa_yz * fz * fga * pb_yy * fx

                     - 1.5 * pa_yz * fz * fga * fx * pb_zz

                     - 3.0 * pa_y * fx * fz * fga * pb_yyz

                     - 3.0 * fx * pa_z * fz * fga * pb_yzz

                     - pa_yzzz * pb_yy * fz * fgb

                     - pa_yzzz * fz * fgb * pb_zz

                     + 9.0 * pa_yz * fz * fx * fx * pb_zz

                     + 18.0 * pa_y * fx * fx * fz * pb_yyz

                     + 18.0 * fx * fx * pa_z * fz * pb_yzz

                     + 7.0 * pa_yzzz * fz * pb_yy * fx

                     + 7.0 * pa_yzzz * fz * fx * pb_zz

                     + 42.0 * pa_yzz * fz * fx * pb_yyz

                     + 14.0 * fx * pa_zzz * fz * pb_yzz

                     - 3.0 * pa_yz * fz * fga * pb_yyzz

                     + 21.0 * pa_yz * fz * fx * pb_yyzz

                     + 16.0 * pa_yzzz * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_yzzz_s_0(double fx,
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
                                     double pb_yzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double s_0_0)
    {
        return s_0_0 * (0.9375 * fx * fx * fx * fx

                     + 1.875 * pa_y * fx * fx * fx * pb_y

                     + 1.125 * fx * fx * fx * pa_zz

                     + 3.375 * fx * fx * fx * pa_z * pb_z

                     + 1.125 * fx * fx * fx * pb_zz

                     + 2.25 * pa_yzz * fx * fx * pb_y

                     + 6.75 * pa_yz * fx * fx * pb_yz

                     + 0.75 * fx * fx * pa_zzz * pb_z

                     + 2.25 * fx * fx * pa_zz * pb_zz

                     + 2.25 * pa_y * fx * fx * pb_yzz

                     + 0.75 * fx * fx * pa_z * pb_zzz

                     + 1.5 * pa_yzzz * pb_yz * fx

                     + 4.5 * pa_yzz * fx * pb_yzz

                     + 0.5 * fx * pa_zzz * pb_zzz

                     + 1.5 * pa_yz * fx * pb_yzzz

                     + pa_yzzz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_yzzz_r_0(double fga,
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
                                     double pb_yzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double r_0_0)
    {
        return r_0_0 * (7.5 * fx * fx * fx * fx * fz

                     - 1.125 * fx * fx * fx * fz * fgb

                     - 1.125 * fx * fx * fx * fz * fga

                     - 2.25 * fx * fx * pa_zz * fz * fgb

                     + 18.75 * pa_y * fx * fx * fx * fz * pb_y

                     + 11.25 * fx * fx * fx * pa_zz * fz

                     + 33.75 * fx * fx * fx * pa_z * fz * pb_z

                     - 2.25 * pa_y * fx * fx * pb_y * fz * fgb

                     - 2.25 * pa_y * fx * fx * fz * fga * pb_y

                     - 2.25 * fx * fx * pa_z * pb_z * fz * fgb

                     - 2.25 * fx * fx * pa_z * fz * fga * pb_z

                     - 2.25 * fx * fx * fz * fga * pb_zz

                     - 4.5 * pa_yzz * fx * pb_y * fz * fgb

                     - 1.5 * fx * pa_zzz * pb_z * fz * fgb

                     + 11.25 * fx * fx * fx * fz * pb_zz

                     + 27.0 * pa_yzz * fz * fx * fx * pb_y

                     + 81.0 * pa_yz * fx * fx * fz * pb_yz

                     + 9.0 * fx * fx * pa_zzz * fz * pb_z

                     + 27.0 * fx * fx * pa_zz * fz * pb_zz

                     - 4.5 * pa_yz * fx * pb_yz * fz * fgb

                     - 4.5 * pa_yz * fz * fga * pb_yz * fx

                     - 4.5 * pa_y * fx * fz * fga * pb_yzz

                     - 1.5 * fx * pa_z * fz * fga * pb_zzz

                     - 3.0 * pa_yzzz * pb_yz * fz * fgb

                     + 27.0 * pa_y * fx * fx * fz * pb_yzz

                     + 9.0 * fx * fx * pa_z * fz * pb_zzz

                     + 21.0 * pa_yzzz * fz * pb_yz * fx

                     + 63.0 * pa_yzz * fz * fx * pb_yzz

                     + 7.0 * fx * pa_zzz * fz * pb_zzz

                     - 3.0 * pa_yz * fz * fga * pb_yzzz

                     + 21.0 * pa_yz * fz * fx * pb_yzzz

                     + 16.0 * pa_yzzz * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_zzzz_s_0(double fx,
                                     double pa_y,
                                     double pa_yz,
                                     double pa_yzz,
                                     double pa_yzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double pb_zzzz,
                                     double s_0_0)
    {
        return s_0_0 * (5.625 * pa_yz * fx * fx * fx

                     + 7.5 * pa_y * fx * fx * fx * pb_z

                     + 0.75 * pa_yzzz * fx * fx

                     + 9.0 * pa_yzz * fx * fx * pb_z

                     + 13.5 * pa_yz * fx * fx * pb_zz

                     + 3.0 * pa_y * fx * fx * pb_zzz

                     + 3.0 * pa_yzzz * pb_zz * fx

                     + 6.0 * pa_yzz * fx * pb_zzz

                     + 1.5 * pa_yz * fx * pb_zzzz

                     + pa_yzzz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_zzzz_r_0(double fga,
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
                                     double pb_zzzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 13.5 * pa_yz * fx * fx * fz * fgb

                     + 56.25 * pa_yz * fx * fx * fx * fz

                     + 75.0 * pa_y * fx * fx * fx * fz * pb_z

                     - 2.25 * pa_yz * fz * fga * fx * fx

                     - 9.0 * pa_y * fx * fx * pb_z * fz * fgb

                     - 9.0 * pa_y * fx * fx * fz * fga * pb_z

                     - 3.0 * pa_yzzz * fx * fz * fgb

                     - 18.0 * pa_yzz * fx * pb_z * fz * fgb

                     + 9.0 * pa_yzzz * fz * fx * fx

                     + 108.0 * pa_yzz * fz * fx * fx * pb_z

                     + 162.0 * pa_yz * fx * fx * fz * pb_zz

                     - 9.0 * pa_yz * fx * pb_zz * fz * fgb

                     - 9.0 * pa_yz * fz * fga * pb_zz * fx

                     - 6.0 * pa_y * fx * fz * fga * pb_zzz

                     - 6.0 * pa_yzzz * pb_zz * fz * fgb

                     + 36.0 * pa_y * fx * fx * fz * pb_zzz

                     + 42.0 * pa_yzzz * fz * pb_zz * fx

                     + 84.0 * pa_yzz * fz * fx * pb_zzz

                     - 3.0 * pa_yz * fz * fga * pb_zzzz

                     + 21.0 * pa_yz * fz * fx * pb_zzzz

                     + 16.0 * pa_yzzz * fz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xxxx_s_0(double fx,
                                     double pa_zz,
                                     double pa_zzzz,
                                     double pb_xx,
                                     double pb_xxxx,
                                     double s_0_0)
    {
        return s_0_0 * (0.5625 * fx * fx * fx * fx

                     + 2.25 * pa_zz * fx * fx * fx

                     + 0.75 * pa_zzzz * fx * fx

                     + 2.25 * fx * fx * fx * pb_xx

                     + 9.0 * pa_zz * fx * fx * pb_xx

                     + 3.0 * pa_zzzz * pb_xx * fx

                     + 0.75 * fx * fx * pb_xxxx

                     + 3.0 * pa_zz * fx * pb_xxxx

                     + pa_zzzz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xxxx_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_zz,
                                     double pa_zzzz,
                                     double pb_xx,
                                     double pb_xxxx,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * fx * fz * fgb

                     - 2.25 * fx * fx * fx * fz * fga

                     - 9.0 * pa_zz * fx * fx * fz * fgb

                     - 4.5 * pa_zz * fz * fga * fx * fx

                     + 4.5 * fx * fx * fx * fx * fz

                     - 3.0 * pa_zzzz * fx * fz * fgb

                     + 22.5 * pa_zz * fz * fx * fx * fx

                     + 9.0 * pa_zzzz * fz * fx * fx

                     - 4.5 * fx * fx * pb_xx * fz * fgb

                     - 9.0 * fx * fx * fz * fga * pb_xx

                     - 18.0 * pa_zz * fx * pb_xx * fz * fgb

                     - 18.0 * pa_zz * fz * fga * pb_xx * fx

                     + 22.5 * fx * fx * fx * fz * pb_xx

                     - 6.0 * pa_zzzz * pb_xx * fz * fgb

                     + 108.0 * pa_zz * fz * fx * fx * pb_xx

                     + 42.0 * pa_zzzz * fz * pb_xx * fx

                     - 3.0 * fx * fz * fga * pb_xxxx

                     - 6.0 * pa_zz * fz * fga * pb_xxxx

                     + 9.0 * fx * fx * fz * pb_xxxx

                     + 42.0 * pa_zz * fz * fx * pb_xxxx

                     + 16.0 * pa_zzzz * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xxxy_s_0(double fx,
                                     double pa_zz,
                                     double pa_zzzz,
                                     double pb_xxxy,
                                     double pb_xy,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pb_xy

                     + 4.5 * pa_zz * fx * fx * pb_xy

                     + 1.5 * pa_zzzz * pb_xy * fx

                     + 0.75 * fx * fx * pb_xxxy

                     + 3.0 * pa_zz * fx * pb_xxxy

                     + pa_zzzz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xxxy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_zz,
                                     double pa_zzzz,
                                     double pb_xxxy,
                                     double pb_xy,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pb_xy * fz * fgb

                     - 4.5 * fx * fx * fz * fga * pb_xy

                     - 9.0 * pa_zz * fx * pb_xy * fz * fgb

                     - 9.0 * pa_zz * fz * fga * pb_xy * fx

                     + 11.25 * fx * fx * fx * fz * pb_xy

                     - 3.0 * pa_zzzz * pb_xy * fz * fgb

                     + 54.0 * pa_zz * fz * fx * fx * pb_xy

                     + 21.0 * pa_zzzz * fz * pb_xy * fx

                     - 3.0 * fx * fz * fga * pb_xxxy

                     - 6.0 * pa_zz * fz * fga * pb_xxxy

                     + 9.0 * fx * fx * fz * pb_xxxy

                     + 42.0 * pa_zz * fz * fx * pb_xxxy

                     + 16.0 * pa_zzzz * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xxxz_s_0(double fx,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pa_zzzz,
                                     double pb_x,
                                     double pb_xxx,
                                     double pb_xxxz,
                                     double pb_xz,
                                     double s_0_0)
    {
        return s_0_0 * (4.5 * pa_z * fx * fx * fx * pb_x

                     + 3.0 * pa_zzz * fx * fx * pb_x

                     + 1.125 * fx * fx * fx * pb_xz

                     + 4.5 * pa_zz * fx * fx * pb_xz

                     + 3.0 * pa_z * fx * fx * pb_xxx

                     + 1.5 * pa_zzzz * pb_xz * fx

                     + 2.0 * pa_zzz * fx * pb_xxx

                     + 0.75 * fx * fx * pb_xxxz

                     + 3.0 * pa_zz * fx * pb_xxxz

                     + pa_zzzz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xxxz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pa_zzzz,
                                     double pb_x,
                                     double pb_xxx,
                                     double pb_xxxz,
                                     double pb_xz,
                                     double r_0_0)
    {
        return r_0_0 * (- 9.0 * pa_z * fx * fx * pb_x * fz * fgb

                     - 9.0 * pa_z * fx * fx * fz * fga * pb_x

                     - 6.0 * pa_zzz * fx * pb_x * fz * fgb

                     + 45.0 * pa_z * fx * fx * fx * fz * pb_x

                     + 36.0 * pa_zzz * fz * fx * fx * pb_x

                     - 2.25 * fx * fx * pb_xz * fz * fgb

                     - 4.5 * fx * fx * fz * fga * pb_xz

                     - 9.0 * pa_zz * fx * pb_xz * fz * fgb

                     - 9.0 * pa_zz * fz * fga * pb_xz * fx

                     - 6.0 * pa_z * fx * fz * fga * pb_xxx

                     + 11.25 * fx * fx * fx * fz * pb_xz

                     - 3.0 * pa_zzzz * pb_xz * fz * fgb

                     + 54.0 * pa_zz * fz * fx * fx * pb_xz

                     + 36.0 * pa_z * fx * fx * fz * pb_xxx

                     + 21.0 * pa_zzzz * fz * pb_xz * fx

                     + 28.0 * pa_zzz * fz * fx * pb_xxx

                     - 3.0 * fx * fz * fga * pb_xxxz

                     - 6.0 * pa_zz * fz * fga * pb_xxxz

                     + 9.0 * fx * fx * fz * pb_xxxz

                     + 42.0 * pa_zz * fz * fx * pb_xxxz

                     + 16.0 * pa_zzzz * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xxyy_s_0(double fx,
                                     double pa_zz,
                                     double pa_zzzz,
                                     double pb_xx,
                                     double pb_xxyy,
                                     double pb_yy,
                                     double s_0_0)
    {
        return s_0_0 * (0.1875 * fx * fx * fx * fx

                     + 0.75 * pa_zz * fx * fx * fx

                     + 0.25 * pa_zzzz * fx * fx

                     + 0.375 * fx * fx * fx * pb_xx

                     + 0.375 * fx * fx * fx * pb_yy

                     + 1.5 * pa_zz * fx * fx * pb_xx

                     + 1.5 * pa_zz * fx * fx * pb_yy

                     + 0.5 * pa_zzzz * pb_xx * fx

                     + 0.5 * pa_zzzz * fx * pb_yy

                     + 0.75 * fx * fx * pb_xxyy

                     + 3.0 * pa_zz * fx * pb_xxyy

                     + pa_zzzz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xxyy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_zz,
                                     double pa_zzzz,
                                     double pb_xx,
                                     double pb_xxyy,
                                     double pb_yy,
                                     double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fx * fz * fgb

                     - 0.75 * fx * fx * fx * fz * fga

                     - 3.0 * pa_zz * fx * fx * fz * fgb

                     - 1.5 * pa_zz * fz * fga * fx * fx

                     + 1.5 * fx * fx * fx * fx * fz

                     - pa_zzzz * fx * fz * fgb

                     + 7.5 * pa_zz * fz * fx * fx * fx

                     + 3.0 * pa_zzzz * fz * fx * fx

                     - 0.75 * fx * fx * pb_xx * fz * fgb

                     - 0.75 * fx * fx * fz * fgb * pb_yy

                     - 1.5 * fx * fx * fz * fga * pb_xx

                     - 1.5 * fx * fx * fz * fga * pb_yy

                     - 3.0 * pa_zz * fx * pb_xx * fz * fgb

                     - 3.0 * pa_zz * fx * fz * fgb * pb_yy

                     - 3.0 * pa_zz * fz * fga * pb_xx * fx

                     - 3.0 * pa_zz * fz * fga * fx * pb_yy

                     + 3.75 * fx * fx * fx * fz * pb_xx

                     + 3.75 * fx * fx * fx * fz * pb_yy

                     - pa_zzzz * pb_xx * fz * fgb

                     - pa_zzzz * fz * fgb * pb_yy

                     + 18.0 * pa_zz * fz * fx * fx * pb_xx

                     + 18.0 * pa_zz * fz * fx * fx * pb_yy

                     + 7.0 * pa_zzzz * fz * pb_xx * fx

                     + 7.0 * pa_zzzz * fz * fx * pb_yy

                     - 3.0 * fx * fz * fga * pb_xxyy

                     - 6.0 * pa_zz * fz * fga * pb_xxyy

                     + 9.0 * fx * fx * fz * pb_xxyy

                     + 42.0 * pa_zz * fz * fx * pb_xxyy

                     + 16.0 * pa_zzzz * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xxyz_s_0(double fx,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pa_zzzz,
                                     double pb_xxy,
                                     double pb_xxyz,
                                     double pb_y,
                                     double pb_yz,
                                     double s_0_0)
    {
        return s_0_0 * (1.5 * pa_z * fx * fx * fx * pb_y

                     + pa_zzz * fx * fx * pb_y

                     + 0.375 * fx * fx * fx * pb_yz

                     + 1.5 * pa_zz * fx * fx * pb_yz

                     + 3.0 * pa_z * fx * fx * pb_xxy

                     + 0.5 * pa_zzzz * fx * pb_yz

                     + 2.0 * pa_zzz * fx * pb_xxy

                     + 0.75 * fx * fx * pb_xxyz

                     + 3.0 * pa_zz * fx * pb_xxyz

                     + pa_zzzz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xxyz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pa_zzzz,
                                     double pb_xxy,
                                     double pb_xxyz,
                                     double pb_y,
                                     double pb_yz,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_z * fx * fx * fz * fgb * pb_y

                     - 3.0 * pa_z * fx * fx * fz * fga * pb_y

                     - 2.0 * pa_zzz * fx * fz * fgb * pb_y

                     + 15.0 * pa_z * fx * fx * fx * fz * pb_y

                     + 12.0 * pa_zzz * fz * fx * fx * pb_y

                     - 0.75 * fx * fx * fz * fgb * pb_yz

                     - 1.5 * fx * fx * fz * fga * pb_yz

                     - 3.0 * pa_zz * fx * fz * fgb * pb_yz

                     - 3.0 * pa_zz * fz * fga * fx * pb_yz

                     - 6.0 * pa_z * fx * fz * fga * pb_xxy

                     + 3.75 * fx * fx * fx * fz * pb_yz

                     - pa_zzzz * fz * fgb * pb_yz

                     + 18.0 * pa_zz * fz * fx * fx * pb_yz

                     + 36.0 * pa_z * fx * fx * fz * pb_xxy

                     + 7.0 * pa_zzzz * fz * fx * pb_yz

                     + 28.0 * pa_zzz * fz * fx * pb_xxy

                     - 3.0 * fx * fz * fga * pb_xxyz

                     - 6.0 * pa_zz * fz * fga * pb_xxyz

                     + 9.0 * fx * fx * fz * pb_xxyz

                     + 42.0 * pa_zz * fz * fx * pb_xxyz

                     + 16.0 * pa_zzzz * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xxzz_s_0(double fx,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pa_zzzz,
                                     double pb_xx,
                                     double pb_xxz,
                                     double pb_xxzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.9375 * fx * fx * fx * fx

                     + 2.25 * pa_zz * fx * fx * fx

                     + 3.0 * pa_z * fx * fx * fx * pb_z

                     + 1.875 * fx * fx * fx * pb_xx

                     + 0.25 * pa_zzzz * fx * fx

                     + 2.0 * pa_zzz * fx * fx * pb_z

                     + 4.5 * pa_zz * fx * fx * pb_xx

                     + 0.375 * fx * fx * fx * pb_zz

                     + 1.5 * pa_zz * fx * fx * pb_zz

                     + 6.0 * pa_z * fx * fx * pb_xxz

                     + 0.5 * pa_zzzz * pb_xx * fx

                     + 0.5 * pa_zzzz * fx * pb_zz

                     + 4.0 * pa_zzz * fx * pb_xxz

                     + 0.75 * fx * fx * pb_xxzz

                     + 3.0 * pa_zz * fx * pb_xxzz

                     + pa_zzzz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xxzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pa_zzzz,
                                     double pb_xx,
                                     double pb_xxz,
                                     double pb_xxzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * fx * fz * fgb

                     - 2.25 * fx * fx * fx * fz * fga

                     - 6.0 * pa_zz * fx * fx * fz * fgb

                     + 7.5 * fx * fx * fx * fx * fz

                     + 22.5 * pa_zz * fx * fx * fx * fz

                     - 1.5 * pa_zz * fz * fga * fx * fx

                     - 6.0 * pa_z * fx * fx * fz * fgb * pb_z

                     - 6.0 * pa_z * fx * fx * fz * fga * pb_z

                     - 4.5 * fx * fx * fz * fga * pb_xx

                     - pa_zzzz * fx * fz * fgb

                     - 4.0 * pa_zzz * fx * fz * fgb * pb_z

                     + 30.0 * pa_z * fx * fx * fx * fz * pb_z

                     + 18.75 * fx * fx * fx * fz * pb_xx

                     + 3.0 * pa_zzzz * fz * fx * fx

                     + 24.0 * pa_zzz * fz * fx * fx * pb_z

                     + 54.0 * pa_zz * fx * fx * fz * pb_xx

                     - 0.75 * fx * fx * pb_xx * fz * fgb

                     - 0.75 * fx * fx * fz * fgb * pb_zz

                     - 1.5 * fx * fx * fz * fga * pb_zz

                     - 3.0 * pa_zz * fx * pb_xx * fz * fgb

                     - 3.0 * pa_zz * fx * fz * fgb * pb_zz

                     - 3.0 * pa_zz * fz * fga * pb_xx * fx

                     - 3.0 * pa_zz * fz * fga * fx * pb_zz

                     - 12.0 * pa_z * fx * fz * fga * pb_xxz

                     + 3.75 * fx * fx * fx * fz * pb_zz

                     - pa_zzzz * pb_xx * fz * fgb

                     - pa_zzzz * fz * fgb * pb_zz

                     + 18.0 * pa_zz * fz * fx * fx * pb_zz

                     + 72.0 * pa_z * fx * fx * fz * pb_xxz

                     + 7.0 * pa_zzzz * fz * pb_xx * fx

                     + 7.0 * pa_zzzz * fz * fx * pb_zz

                     + 56.0 * pa_zzz * fz * fx * pb_xxz

                     - 3.0 * fx * fz * fga * pb_xxzz

                     - 6.0 * pa_zz * fz * fga * pb_xxzz

                     + 9.0 * fx * fx * fz * pb_xxzz

                     + 42.0 * pa_zz * fz * fx * pb_xxzz

                     + 16.0 * pa_zzzz * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xyyy_s_0(double fx,
                                     double pa_zz,
                                     double pa_zzzz,
                                     double pb_xy,
                                     double pb_xyyy,
                                     double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pb_xy

                     + 4.5 * pa_zz * fx * fx * pb_xy

                     + 1.5 * pa_zzzz * pb_xy * fx

                     + 0.75 * fx * fx * pb_xyyy

                     + 3.0 * pa_zz * fx * pb_xyyy

                     + pa_zzzz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xyyy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_zz,
                                     double pa_zzzz,
                                     double pb_xy,
                                     double pb_xyyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pb_xy * fz * fgb

                     - 4.5 * fx * fx * fz * fga * pb_xy

                     - 9.0 * pa_zz * fx * pb_xy * fz * fgb

                     - 9.0 * pa_zz * fz * fga * pb_xy * fx

                     + 11.25 * fx * fx * fx * fz * pb_xy

                     - 3.0 * pa_zzzz * pb_xy * fz * fgb

                     + 54.0 * pa_zz * fz * fx * fx * pb_xy

                     + 21.0 * pa_zzzz * fz * pb_xy * fx

                     - 3.0 * fx * fz * fga * pb_xyyy

                     - 6.0 * pa_zz * fz * fga * pb_xyyy

                     + 9.0 * fx * fx * fz * pb_xyyy

                     + 42.0 * pa_zz * fz * fx * pb_xyyy

                     + 16.0 * pa_zzzz * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xyyz_s_0(double fx,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pa_zzzz,
                                     double pb_x,
                                     double pb_xyy,
                                     double pb_xyyz,
                                     double pb_xz,
                                     double s_0_0)
    {
        return s_0_0 * (1.5 * pa_z * fx * fx * fx * pb_x

                     + pa_zzz * fx * fx * pb_x

                     + 0.375 * fx * fx * fx * pb_xz

                     + 1.5 * pa_zz * fx * fx * pb_xz

                     + 3.0 * pa_z * fx * fx * pb_xyy

                     + 0.5 * pa_zzzz * pb_xz * fx

                     + 2.0 * pa_zzz * fx * pb_xyy

                     + 0.75 * fx * fx * pb_xyyz

                     + 3.0 * pa_zz * fx * pb_xyyz

                     + pa_zzzz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xyyz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pa_zzzz,
                                     double pb_x,
                                     double pb_xyy,
                                     double pb_xyyz,
                                     double pb_xz,
                                     double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_z * fx * fx * pb_x * fz * fgb

                     - 3.0 * pa_z * fx * fx * fz * fga * pb_x

                     - 2.0 * pa_zzz * fx * pb_x * fz * fgb

                     + 15.0 * pa_z * fx * fx * fx * fz * pb_x

                     + 12.0 * pa_zzz * fz * fx * fx * pb_x

                     - 0.75 * fx * fx * pb_xz * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_xz

                     - 3.0 * pa_zz * fx * pb_xz * fz * fgb

                     - 3.0 * pa_zz * fz * fga * pb_xz * fx

                     - 6.0 * pa_z * fx * fz * fga * pb_xyy

                     + 3.75 * fx * fx * fx * fz * pb_xz

                     - pa_zzzz * pb_xz * fz * fgb

                     + 18.0 * pa_zz * fz * fx * fx * pb_xz

                     + 36.0 * pa_z * fx * fx * fz * pb_xyy

                     + 7.0 * pa_zzzz * fz * pb_xz * fx

                     + 28.0 * pa_zzz * fz * fx * pb_xyy

                     - 3.0 * fx * fz * fga * pb_xyyz

                     - 6.0 * pa_zz * fz * fga * pb_xyyz

                     + 9.0 * fx * fx * fz * pb_xyyz

                     + 42.0 * pa_zz * fz * fx * pb_xyyz

                     + 16.0 * pa_zzzz * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xyzz_s_0(double fx,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pa_zzzz,
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_xyzz,
                                     double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pb_xy

                     + 4.5 * pa_zz * fx * fx * pb_xy

                     + 6.0 * pa_z * fx * fx * pb_xyz

                     + 0.5 * pa_zzzz * pb_xy * fx

                     + 4.0 * pa_zzz * fx * pb_xyz

                     + 0.75 * fx * fx * pb_xyzz

                     + 3.0 * pa_zz * fx * pb_xyzz

                     + pa_zzzz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xyzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pa_zzzz,
                                     double pb_xy,
                                     double pb_xyz,
                                     double pb_xyzz,
                                     double r_0_0)
    {
        return r_0_0 * (- 4.5 * fx * fx * fz * fga * pb_xy

                     + 18.75 * fx * fx * fx * fz * pb_xy

                     + 54.0 * pa_zz * fx * fx * fz * pb_xy

                     - 0.75 * fx * fx * pb_xy * fz * fgb

                     - 3.0 * pa_zz * fx * pb_xy * fz * fgb

                     - 3.0 * pa_zz * fz * fga * pb_xy * fx

                     - 12.0 * pa_z * fx * fz * fga * pb_xyz

                     - pa_zzzz * pb_xy * fz * fgb

                     + 72.0 * pa_z * fx * fx * fz * pb_xyz

                     + 7.0 * pa_zzzz * fz * pb_xy * fx

                     + 56.0 * pa_zzz * fz * fx * pb_xyz

                     - 3.0 * fx * fz * fga * pb_xyzz

                     - 6.0 * pa_zz * fz * fga * pb_xyzz

                     + 9.0 * fx * fx * fz * pb_xyzz

                     + 42.0 * pa_zz * fz * fx * pb_xyzz

                     + 16.0 * pa_zzzz * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xzzz_s_0(double fx,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pa_zzzz,
                                     double pb_x,
                                     double pb_xz,
                                     double pb_xzz,
                                     double pb_xzzz,
                                     double s_0_0)
    {
        return s_0_0 * (7.5 * pa_z * fx * fx * fx * pb_x

                     + 5.625 * fx * fx * fx * pb_xz

                     + 3.0 * pa_zzz * fx * fx * pb_x

                     + 13.5 * pa_zz * fx * fx * pb_xz

                     + 9.0 * pa_z * fx * fx * pb_xzz

                     + 1.5 * pa_zzzz * pb_xz * fx

                     + 6.0 * pa_zzz * fx * pb_xzz

                     + 0.75 * fx * fx * pb_xzzz

                     + 3.0 * pa_zz * fx * pb_xzzz

                     + pa_zzzz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xzzz_r_0(double fga,
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
                                     double pb_xzzz,
                                     double r_0_0)
    {
        return r_0_0 * (75.0 * pa_z * fx * fx * fx * fz * pb_x

                     - 9.0 * pa_z * fx * fx * pb_x * fz * fgb

                     - 9.0 * pa_z * fx * fx * fz * fga * pb_x

                     - 13.5 * fx * fx * fz * fga * pb_xz

                     - 6.0 * pa_zzz * fx * pb_x * fz * fgb

                     + 56.25 * fx * fx * fx * fz * pb_xz

                     + 36.0 * pa_zzz * fz * fx * fx * pb_x

                     + 162.0 * pa_zz * fx * fx * fz * pb_xz

                     - 2.25 * fx * fx * pb_xz * fz * fgb

                     - 9.0 * pa_zz * fx * pb_xz * fz * fgb

                     - 9.0 * pa_zz * fz * fga * pb_xz * fx

                     - 18.0 * pa_z * fx * fz * fga * pb_xzz

                     - 3.0 * pa_zzzz * pb_xz * fz * fgb

                     + 108.0 * pa_z * fx * fx * fz * pb_xzz

                     + 21.0 * pa_zzzz * fz * pb_xz * fx

                     + 84.0 * pa_zzz * fz * fx * pb_xzz

                     - 3.0 * fx * fz * fga * pb_xzzz

                     - 6.0 * pa_zz * fz * fga * pb_xzzz

                     + 9.0 * fx * fx * fz * pb_xzzz

                     + 42.0 * pa_zz * fz * fx * pb_xzzz

                     + 16.0 * pa_zzzz * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_yyyy_s_0(double fx,
                                     double pa_zz,
                                     double pa_zzzz,
                                     double pb_yy,
                                     double pb_yyyy,
                                     double s_0_0)
    {
        return s_0_0 * (0.5625 * fx * fx * fx * fx

                     + 2.25 * pa_zz * fx * fx * fx

                     + 0.75 * pa_zzzz * fx * fx

                     + 2.25 * fx * fx * fx * pb_yy

                     + 9.0 * pa_zz * fx * fx * pb_yy

                     + 3.0 * pa_zzzz * pb_yy * fx

                     + 0.75 * fx * fx * pb_yyyy

                     + 3.0 * pa_zz * fx * pb_yyyy

                     + pa_zzzz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_yyyy_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_zz,
                                     double pa_zzzz,
                                     double pb_yy,
                                     double pb_yyyy,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * fx * fz * fgb

                     - 2.25 * fx * fx * fx * fz * fga

                     - 9.0 * pa_zz * fx * fx * fz * fgb

                     - 4.5 * pa_zz * fz * fga * fx * fx

                     + 4.5 * fx * fx * fx * fx * fz

                     - 3.0 * pa_zzzz * fx * fz * fgb

                     + 22.5 * pa_zz * fz * fx * fx * fx

                     + 9.0 * pa_zzzz * fz * fx * fx

                     - 4.5 * fx * fx * pb_yy * fz * fgb

                     - 9.0 * fx * fx * fz * fga * pb_yy

                     - 18.0 * pa_zz * fx * pb_yy * fz * fgb

                     - 18.0 * pa_zz * fz * fga * pb_yy * fx

                     + 22.5 * fx * fx * fx * fz * pb_yy

                     - 6.0 * pa_zzzz * pb_yy * fz * fgb

                     + 108.0 * pa_zz * fz * fx * fx * pb_yy

                     + 42.0 * pa_zzzz * fz * pb_yy * fx

                     - 3.0 * fx * fz * fga * pb_yyyy

                     - 6.0 * pa_zz * fz * fga * pb_yyyy

                     + 9.0 * fx * fx * fz * pb_yyyy

                     + 42.0 * pa_zz * fz * fx * pb_yyyy

                     + 16.0 * pa_zzzz * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_yyyz_s_0(double fx,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pa_zzzz,
                                     double pb_y,
                                     double pb_yyy,
                                     double pb_yyyz,
                                     double pb_yz,
                                     double s_0_0)
    {
        return s_0_0 * (4.5 * pa_z * fx * fx * fx * pb_y

                     + 3.0 * pa_zzz * fx * fx * pb_y

                     + 1.125 * fx * fx * fx * pb_yz

                     + 4.5 * pa_zz * fx * fx * pb_yz

                     + 3.0 * pa_z * fx * fx * pb_yyy

                     + 1.5 * pa_zzzz * pb_yz * fx

                     + 2.0 * pa_zzz * fx * pb_yyy

                     + 0.75 * fx * fx * pb_yyyz

                     + 3.0 * pa_zz * fx * pb_yyyz

                     + pa_zzzz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_yyyz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pa_zzzz,
                                     double pb_y,
                                     double pb_yyy,
                                     double pb_yyyz,
                                     double pb_yz,
                                     double r_0_0)
    {
        return r_0_0 * (- 9.0 * pa_z * fx * fx * pb_y * fz * fgb

                     - 9.0 * pa_z * fx * fx * fz * fga * pb_y

                     - 6.0 * pa_zzz * fx * pb_y * fz * fgb

                     + 45.0 * pa_z * fx * fx * fx * fz * pb_y

                     + 36.0 * pa_zzz * fz * fx * fx * pb_y

                     - 2.25 * fx * fx * pb_yz * fz * fgb

                     - 4.5 * fx * fx * fz * fga * pb_yz

                     - 9.0 * pa_zz * fx * pb_yz * fz * fgb

                     - 9.0 * pa_zz * fz * fga * pb_yz * fx

                     - 6.0 * pa_z * fx * fz * fga * pb_yyy

                     + 11.25 * fx * fx * fx * fz * pb_yz

                     - 3.0 * pa_zzzz * pb_yz * fz * fgb

                     + 54.0 * pa_zz * fz * fx * fx * pb_yz

                     + 36.0 * pa_z * fx * fx * fz * pb_yyy

                     + 21.0 * pa_zzzz * fz * pb_yz * fx

                     + 28.0 * pa_zzz * fz * fx * pb_yyy

                     - 3.0 * fx * fz * fga * pb_yyyz

                     - 6.0 * pa_zz * fz * fga * pb_yyyz

                     + 9.0 * fx * fx * fz * pb_yyyz

                     + 42.0 * pa_zz * fz * fx * pb_yyyz

                     + 16.0 * pa_zzzz * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_yyzz_s_0(double fx,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pa_zzzz,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_yyzz,
                                     double pb_z,
                                     double pb_zz,
                                     double s_0_0)
    {
        return s_0_0 * (0.9375 * fx * fx * fx * fx

                     + 2.25 * pa_zz * fx * fx * fx

                     + 3.0 * pa_z * fx * fx * fx * pb_z

                     + 1.875 * fx * fx * fx * pb_yy

                     + 0.25 * pa_zzzz * fx * fx

                     + 2.0 * pa_zzz * fx * fx * pb_z

                     + 4.5 * pa_zz * fx * fx * pb_yy

                     + 0.375 * fx * fx * fx * pb_zz

                     + 1.5 * pa_zz * fx * fx * pb_zz

                     + 6.0 * pa_z * fx * fx * pb_yyz

                     + 0.5 * pa_zzzz * pb_yy * fx

                     + 0.5 * pa_zzzz * fx * pb_zz

                     + 4.0 * pa_zzz * fx * pb_yyz

                     + 0.75 * fx * fx * pb_yyzz

                     + 3.0 * pa_zz * fx * pb_yyzz

                     + pa_zzzz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_yyzz_r_0(double fga,
                                     double fgb,
                                     double fx,
                                     double fz,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pa_zzzz,
                                     double pb_yy,
                                     double pb_yyz,
                                     double pb_yyzz,
                                     double pb_z,
                                     double pb_zz,
                                     double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * fx * fz * fgb

                     - 2.25 * fx * fx * fx * fz * fga

                     - 6.0 * pa_zz * fx * fx * fz * fgb

                     + 7.5 * fx * fx * fx * fx * fz

                     + 22.5 * pa_zz * fx * fx * fx * fz

                     - 1.5 * pa_zz * fz * fga * fx * fx

                     - 6.0 * pa_z * fx * fx * fz * fgb * pb_z

                     - 6.0 * pa_z * fx * fx * fz * fga * pb_z

                     - 4.5 * fx * fx * fz * fga * pb_yy

                     - pa_zzzz * fx * fz * fgb

                     - 4.0 * pa_zzz * fx * fz * fgb * pb_z

                     + 30.0 * pa_z * fx * fx * fx * fz * pb_z

                     + 18.75 * fx * fx * fx * fz * pb_yy

                     + 3.0 * pa_zzzz * fz * fx * fx

                     + 24.0 * pa_zzz * fz * fx * fx * pb_z

                     + 54.0 * pa_zz * fx * fx * fz * pb_yy

                     - 0.75 * fx * fx * pb_yy * fz * fgb

                     - 0.75 * fx * fx * fz * fgb * pb_zz

                     - 1.5 * fx * fx * fz * fga * pb_zz

                     - 3.0 * pa_zz * fx * pb_yy * fz * fgb

                     - 3.0 * pa_zz * fx * fz * fgb * pb_zz

                     - 3.0 * pa_zz * fz * fga * pb_yy * fx

                     - 3.0 * pa_zz * fz * fga * fx * pb_zz

                     - 12.0 * pa_z * fx * fz * fga * pb_yyz

                     + 3.75 * fx * fx * fx * fz * pb_zz

                     - pa_zzzz * pb_yy * fz * fgb

                     - pa_zzzz * fz * fgb * pb_zz

                     + 18.0 * pa_zz * fz * fx * fx * pb_zz

                     + 72.0 * pa_z * fx * fx * fz * pb_yyz

                     + 7.0 * pa_zzzz * fz * pb_yy * fx

                     + 7.0 * pa_zzzz * fz * fx * pb_zz

                     + 56.0 * pa_zzz * fz * fx * pb_yyz

                     - 3.0 * fx * fz * fga * pb_yyzz

                     - 6.0 * pa_zz * fz * fga * pb_yyzz

                     + 9.0 * fx * fx * fz * pb_yyzz

                     + 42.0 * pa_zz * fz * fx * pb_yyzz

                     + 16.0 * pa_zzzz * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_yzzz_s_0(double fx,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pa_zzzz,
                                     double pb_y,
                                     double pb_yz,
                                     double pb_yzz,
                                     double pb_yzzz,
                                     double s_0_0)
    {
        return s_0_0 * (7.5 * pa_z * fx * fx * fx * pb_y

                     + 5.625 * fx * fx * fx * pb_yz

                     + 3.0 * pa_zzz * fx * fx * pb_y

                     + 13.5 * pa_zz * fx * fx * pb_yz

                     + 9.0 * pa_z * fx * fx * pb_yzz

                     + 1.5 * pa_zzzz * pb_yz * fx

                     + 6.0 * pa_zzz * fx * pb_yzz

                     + 0.75 * fx * fx * pb_yzzz

                     + 3.0 * pa_zz * fx * pb_yzzz

                     + pa_zzzz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_yzzz_r_0(double fga,
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
                                     double pb_yzzz,
                                     double r_0_0)
    {
        return r_0_0 * (75.0 * pa_z * fx * fx * fx * fz * pb_y

                     - 9.0 * pa_z * fx * fx * pb_y * fz * fgb

                     - 9.0 * pa_z * fx * fx * fz * fga * pb_y

                     - 13.5 * fx * fx * fz * fga * pb_yz

                     - 6.0 * pa_zzz * fx * pb_y * fz * fgb

                     + 56.25 * fx * fx * fx * fz * pb_yz

                     + 36.0 * pa_zzz * fz * fx * fx * pb_y

                     + 162.0 * pa_zz * fx * fx * fz * pb_yz

                     - 2.25 * fx * fx * pb_yz * fz * fgb

                     - 9.0 * pa_zz * fx * pb_yz * fz * fgb

                     - 9.0 * pa_zz * fz * fga * pb_yz * fx

                     - 18.0 * pa_z * fx * fz * fga * pb_yzz

                     - 3.0 * pa_zzzz * pb_yz * fz * fgb

                     + 108.0 * pa_z * fx * fx * fz * pb_yzz

                     + 21.0 * pa_zzzz * fz * pb_yz * fx

                     + 84.0 * pa_zzz * fz * fx * pb_yzz

                     - 3.0 * fx * fz * fga * pb_yzzz

                     - 6.0 * pa_zz * fz * fga * pb_yzzz

                     + 9.0 * fx * fx * fz * pb_yzzz

                     + 42.0 * pa_zz * fz * fx * pb_yzzz

                     + 16.0 * pa_zzzz * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_zzzz_s_0(double fx,
                                     double pa_z,
                                     double pa_zz,
                                     double pa_zzz,
                                     double pa_zzzz,
                                     double pb_z,
                                     double pb_zz,
                                     double pb_zzz,
                                     double pb_zzzz,
                                     double s_0_0)
    {
        return s_0_0 * (6.5625 * fx * fx * fx * fx

                     + 11.25 * pa_zz * fx * fx * fx

                     + 30.0 * pa_z * fx * fx * fx * pb_z

                     + 11.25 * fx * fx * fx * pb_zz

                     + 0.75 * pa_zzzz * fx * fx

                     + 12.0 * pa_zzz * fx * fx * pb_z

                     + 27.0 * pa_zz * fx * fx * pb_zz

                     + 12.0 * pa_z * fx * fx * pb_zzz

                     + 3.0 * pa_zzzz * pb_zz * fx

                     + 8.0 * pa_zzz * fx * pb_zzz

                     + 0.75 * fx * fx * pb_zzzz

                     + 3.0 * pa_zz * fx * pb_zzzz

                     + pa_zzzz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_zzzz_r_0(double fga,
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
                                     double pb_zzzz,
                                     double r_0_0)
    {
        return r_0_0 * (52.5 * fx * fx * fx * fx * fz

                     - 11.25 * fx * fx * fx * fz * fgb

                     - 11.25 * fx * fx * fx * fz * fga

                     - 27.0 * pa_zz * fx * fx * fz * fgb

                     + 112.5 * pa_zz * fx * fx * fx * fz

                     + 300.0 * pa_z * fx * fx * fx * fz * pb_z

                     - 4.5 * pa_zz * fz * fga * fx * fx

                     - 36.0 * pa_z * fx * fx * pb_z * fz * fgb

                     - 36.0 * pa_z * fx * fx * fz * fga * pb_z

                     - 27.0 * fx * fx * fz * fga * pb_zz

                     - 3.0 * pa_zzzz * fx * fz * fgb

                     - 24.0 * pa_zzz * fx * pb_z * fz * fgb

                     + 112.5 * fx * fx * fx * fz * pb_zz

                     + 9.0 * pa_zzzz * fz * fx * fx

                     + 144.0 * pa_zzz * fz * fx * fx * pb_z

                     + 324.0 * pa_zz * fx * fx * fz * pb_zz

                     - 4.5 * fx * fx * pb_zz * fz * fgb

                     - 18.0 * pa_zz * fx * pb_zz * fz * fgb

                     - 18.0 * pa_zz * fz * fga * pb_zz * fx

                     - 24.0 * pa_z * fx * fz * fga * pb_zzz

                     - 6.0 * pa_zzzz * pb_zz * fz * fgb

                     + 144.0 * pa_z * fx * fx * fz * pb_zzz

                     + 42.0 * pa_zzzz * fz * pb_zz * fx

                     + 112.0 * pa_zzz * fz * fx * pb_zzz

                     - 3.0 * fx * fz * fga * pb_zzzz

                     - 6.0 * pa_zz * fz * fga * pb_zzzz

                     + 9.0 * fx * fx * fz * pb_zzzz

                     + 42.0 * pa_zz * fz * fx * pb_zzzz

                     + 16.0 * pa_zzzz * fz * pb_zzzz);

    }


} // kinrecfunc namespace

