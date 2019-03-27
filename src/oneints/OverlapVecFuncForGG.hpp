//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

namespace ovlvecfunc { // ovlvecfunc namespace

    // SIMD elementary functions for (G||G) integrals

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


} // ovlrecfunc namespace

