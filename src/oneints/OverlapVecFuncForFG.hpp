//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

namespace ovlvecfunc { // ovlvecfunc namespace

    // SIMD elementary functions for (F||G) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xxxx_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double pb_xxxx,
                                    double s_0_0)
    {
        return s_0_0 * (5.625 * pa_x * fx * fx * fx

                     + 7.5 * fx * fx * fx * pb_x

                     + 0.75 * pa_xxx * fx * fx

                     + 9.0 * pa_xx * fx * fx * pb_x

                     + 13.5 * pa_x * fx * fx * pb_xx

                     + 3.0 * fx * fx * pb_xxx

                     + 3.0 * pa_xxx * pb_xx * fx

                     + 6.0 * pa_xx * fx * pb_xxx

                     + 1.5 * pa_x * fx * pb_xxxx

                     + pa_xxx * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xxxy_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pb_xxxy,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pb_y

                     + 2.25 * pa_xx * fx * fx * pb_y

                     + 6.75 * pa_x * fx * fx * pb_xy

                     + 2.25 * fx * fx * pb_xxy

                     + 1.5 * pa_xxx * pb_xy * fx

                     + 4.5 * pa_xx * fx * pb_xxy

                     + 1.5 * pa_x * fx * pb_xxxy

                     + pa_xxx * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xxxz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pb_xxxz,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pb_z

                     + 2.25 * pa_xx * fx * fx * pb_z

                     + 6.75 * pa_x * fx * fx * pb_xz

                     + 2.25 * fx * fx * pb_xxz

                     + 1.5 * pa_xxx * pb_xz * fx

                     + 4.5 * pa_xx * fx * pb_xxz

                     + 1.5 * pa_x * fx * pb_xxxz

                     + pa_xxx * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xxyy_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxyy,
                                    double pb_xyy,
                                    double pb_yy,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * pa_x * fx * fx * fx

                     + 0.75 * fx * fx * fx * pb_x

                     + 0.25 * pa_xxx * fx * fx

                     + 1.5 * pa_xx * fx * fx * pb_x

                     + 2.25 * pa_x * fx * fx * pb_yy

                     + 0.75 * pa_x * fx * fx * pb_xx

                     + 1.5 * fx * fx * pb_xyy

                     + 0.5 * pa_xxx * pb_xx * fx

                     + 0.5 * pa_xxx * fx * pb_yy

                     + 3.0 * pa_xx * fx * pb_xyy

                     + 1.5 * pa_x * fx * pb_xxyy

                     + pa_xxx * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xxyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pb_xxyz,
                                    double pb_xyz,
                                    double pb_yz,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_x * fx * fx * pb_yz

                     + 1.5 * fx * fx * pb_xyz

                     + 0.5 * pa_xxx * fx * pb_yz

                     + 3.0 * pa_xx * fx * pb_xyz

                     + 1.5 * pa_x * fx * pb_xxyz

                     + pa_xxx * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xxzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxzz,
                                    double pb_xzz,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * pa_x * fx * fx * fx

                     + 0.75 * fx * fx * fx * pb_x

                     + 0.25 * pa_xxx * fx * fx

                     + 1.5 * pa_xx * fx * fx * pb_x

                     + 2.25 * pa_x * fx * fx * pb_zz

                     + 0.75 * pa_x * fx * fx * pb_xx

                     + 1.5 * fx * fx * pb_xzz

                     + 0.5 * pa_xxx * pb_xx * fx

                     + 0.5 * pa_xxx * fx * pb_zz

                     + 3.0 * pa_xx * fx * pb_xzz

                     + 1.5 * pa_x * fx * pb_xxzz

                     + pa_xxx * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xyyy_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pb_xy,
                                    double pb_xyyy,
                                    double pb_y,
                                    double pb_yyy,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pb_y

                     + 2.25 * pa_xx * fx * fx * pb_y

                     + 2.25 * pa_x * fx * fx * pb_xy

                     + 0.75 * fx * fx * pb_yyy

                     + 1.5 * pa_xxx * pb_xy * fx

                     + 1.5 * pa_xx * fx * pb_yyy

                     + 1.5 * pa_x * fx * pb_xyyy

                     + pa_xxx * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xyyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pb_xyyz,
                                    double pb_xz,
                                    double pb_yyz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_z

                     + 0.75 * pa_xx * fx * fx * pb_z

                     + 0.75 * pa_x * fx * fx * pb_xz

                     + 0.75 * fx * fx * pb_yyz

                     + 0.5 * pa_xxx * pb_xz * fx

                     + 1.5 * pa_xx * fx * pb_yyz

                     + 1.5 * pa_x * fx * pb_xyyz

                     + pa_xxx * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xyzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pb_xy,
                                    double pb_xyzz,
                                    double pb_y,
                                    double pb_yzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_y

                     + 0.75 * pa_xx * fx * fx * pb_y

                     + 0.75 * pa_x * fx * fx * pb_xy

                     + 0.75 * fx * fx * pb_yzz

                     + 0.5 * pa_xxx * pb_xy * fx

                     + 1.5 * pa_xx * fx * pb_yzz

                     + 1.5 * pa_x * fx * pb_xyzz

                     + pa_xxx * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xzzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pb_xz,
                                    double pb_xzzz,
                                    double pb_z,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pb_z

                     + 2.25 * pa_xx * fx * fx * pb_z

                     + 2.25 * pa_x * fx * fx * pb_xz

                     + 0.75 * fx * fx * pb_zzz

                     + 1.5 * pa_xxx * pb_xz * fx

                     + 1.5 * pa_xx * fx * pb_zzz

                     + 1.5 * pa_x * fx * pb_xzzz

                     + pa_xxx * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_yyyy_s_0(double fx,
                                    double pa_x,
                                    double pa_xxx,
                                    double pb_yy,
                                    double pb_yyyy,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * pa_x * fx * fx * fx

                     + 0.75 * pa_xxx * fx * fx

                     + 4.5 * pa_x * fx * fx * pb_yy

                     + 3.0 * pa_xxx * pb_yy * fx

                     + 1.5 * pa_x * fx * pb_yyyy

                     + pa_xxx * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_yyyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xxx,
                                    double pb_yyyz,
                                    double pb_yz,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_x * fx * fx * pb_yz

                     + 1.5 * pa_xxx * pb_yz * fx

                     + 1.5 * pa_x * fx * pb_yyyz

                     + pa_xxx * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_yyzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xxx,
                                    double pb_yy,
                                    double pb_yyzz,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 0.25 * pa_xxx * fx * fx

                     + 0.75 * pa_x * fx * fx * pb_yy

                     + 0.75 * pa_x * fx * fx * pb_zz

                     + 0.5 * pa_xxx * pb_yy * fx

                     + 0.5 * pa_xxx * fx * pb_zz

                     + 1.5 * pa_x * fx * pb_yyzz

                     + pa_xxx * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_yzzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xxx,
                                    double pb_yz,
                                    double pb_yzzz,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_x * fx * fx * pb_yz

                     + 1.5 * pa_xxx * pb_yz * fx

                     + 1.5 * pa_x * fx * pb_yzzz

                     + pa_xxx * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_zzzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xxx,
                                    double pb_zz,
                                    double pb_zzzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * pa_x * fx * fx * fx

                     + 0.75 * pa_xxx * fx * fx

                     + 4.5 * pa_x * fx * fx * pb_zz

                     + 3.0 * pa_xxx * pb_zz * fx

                     + 1.5 * pa_x * fx * pb_zzzz

                     + pa_xxx * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xxxx_s_0(double fx,
                                    double pa_xxy,
                                    double pa_xy,
                                    double pa_y,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double pb_xxxx,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pa_y

                     + 0.75 * pa_xxy * fx * fx

                     + 6.0 * pa_xy * fx * fx * pb_x

                     + 4.5 * fx * fx * pa_y * pb_xx

                     + 3.0 * pa_xxy * pb_xx * fx

                     + 4.0 * pa_xy * fx * pb_xxx

                     + 0.5 * fx * pa_y * pb_xxxx

                     + pa_xxy * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xxxy_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
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
        return s_0_0 * (0.75 * pa_x * fx * fx * fx

                     + 1.125 * fx * fx * fx * pb_x

                     + 0.75 * pa_xx * fx * fx * pb_x

                     + 1.5 * pa_xy * fx * fx * pb_y

                     + 1.5 * pa_x * fx * fx * pb_xx

                     + 2.25 * fx * fx * pa_y * pb_xy

                     + 0.25 * fx * fx * pb_xxx

                     + 1.5 * pa_xxy * pb_xy * fx

                     + 0.5 * pa_xx * fx * pb_xxx

                     + 3.0 * pa_xy * fx * pb_xxy

                     + 0.5 * fx * pa_y * pb_xxxy

                     + pa_xxy * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xxxz_s_0(double fx,
                                    double pa_xxy,
                                    double pa_xy,
                                    double pa_y,
                                    double pb_xxxz,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xy * fx * fx * pb_z

                     + 2.25 * fx * fx * pa_y * pb_xz

                     + 1.5 * pa_xxy * pb_xz * fx

                     + 3.0 * pa_xy * fx * pb_xxz

                     + 0.5 * fx * pa_y * pb_xxxz

                     + pa_xxy * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xxyy_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
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
        return s_0_0 * (0.375 * fx * fx * fx * pa_y

                     + 0.75 * fx * fx * fx * pb_y

                     + 0.25 * pa_xxy * fx * fx

                     + 0.5 * pa_xx * fx * fx * pb_y

                     + pa_xy * fx * fx * pb_x

                     + 2.0 * pa_x * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_y * pb_yy

                     + 0.25 * fx * fx * pa_y * pb_xx

                     + 0.5 * fx * fx * pb_xxy

                     + 0.5 * pa_xxy * pb_xx * fx

                     + 0.5 * pa_xxy * fx * pb_yy

                     + pa_xx * fx * pb_xxy

                     + 2.0 * pa_xy * fx * pb_xyy

                     + 0.5 * fx * pa_y * pb_xxyy

                     + pa_xxy * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xxyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
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
        return s_0_0 * (0.375 * fx * fx * fx * pb_z

                     + 0.25 * pa_xx * fx * fx * pb_z

                     + pa_x * fx * fx * pb_xz

                     + 0.75 * fx * fx * pa_y * pb_yz

                     + 0.25 * fx * fx * pb_xxz

                     + 0.5 * pa_xxy * fx * pb_yz

                     + 0.5 * pa_xx * fx * pb_xxz

                     + 2.0 * pa_xy * fx * pb_xyz

                     + 0.5 * fx * pa_y * pb_xxyz

                     + pa_xxy * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xxzz_s_0(double fx,
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
        return s_0_0 * (0.375 * fx * fx * fx * pa_y

                     + 0.25 * pa_xxy * fx * fx

                     + pa_xy * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_y * pb_zz

                     + 0.25 * fx * fx * pa_y * pb_xx

                     + 0.5 * pa_xxy * pb_xx * fx

                     + 0.5 * pa_xxy * fx * pb_zz

                     + 2.0 * pa_xy * fx * pb_xzz

                     + 0.5 * fx * pa_y * pb_xxzz

                     + pa_xxy * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xyyy_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
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
        return s_0_0 * (0.75 * pa_x * fx * fx * fx

                     + 0.375 * fx * fx * fx * pb_x

                     + 0.75 * pa_xx * fx * fx * pb_x

                     + 1.5 * pa_xy * fx * fx * pb_y

                     + 1.5 * pa_x * fx * fx * pb_yy

                     + 0.75 * fx * fx * pa_y * pb_xy

                     + 0.75 * fx * fx * pb_xyy

                     + 1.5 * pa_xxy * pb_xy * fx

                     + 1.5 * pa_xx * fx * pb_xyy

                     + pa_xy * fx * pb_yyy

                     + 0.5 * fx * pa_y * pb_xyyy

                     + pa_xxy * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xyyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
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
        return s_0_0 * (0.5 * pa_xy * fx * fx * pb_z

                     + pa_x * fx * fx * pb_yz

                     + 0.25 * fx * fx * pa_y * pb_xz

                     + 0.5 * fx * fx * pb_xyz

                     + 0.5 * pa_xxy * pb_xz * fx

                     + pa_xx * fx * pb_xyz

                     + pa_xy * fx * pb_yyz

                     + 0.5 * fx * pa_y * pb_xyyz

                     + pa_xxy * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xyzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
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
        return s_0_0 * (0.25 * pa_x * fx * fx * fx

                     + 0.125 * fx * fx * fx * pb_x

                     + 0.25 * pa_xx * fx * fx * pb_x

                     + 0.5 * pa_xy * fx * fx * pb_y

                     + 0.5 * pa_x * fx * fx * pb_zz

                     + 0.25 * fx * fx * pa_y * pb_xy

                     + 0.25 * fx * fx * pb_xzz

                     + 0.5 * pa_xxy * pb_xy * fx

                     + 0.5 * pa_xx * fx * pb_xzz

                     + pa_xy * fx * pb_yzz

                     + 0.5 * fx * pa_y * pb_xyzz

                     + pa_xxy * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xzzz_s_0(double fx,
                                    double pa_xxy,
                                    double pa_xy,
                                    double pa_y,
                                    double pb_xz,
                                    double pb_xzzz,
                                    double pb_z,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xy * fx * fx * pb_z

                     + 0.75 * fx * fx * pa_y * pb_xz

                     + 1.5 * pa_xxy * pb_xz * fx

                     + pa_xy * fx * pb_zzz

                     + 0.5 * fx * pa_y * pb_xzzz

                     + pa_xxy * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_yyyy_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_y,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double pb_yyyy,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_y

                     + 1.5 * fx * fx * fx * pb_y

                     + 0.75 * pa_xxy * fx * fx

                     + 3.0 * pa_xx * fx * fx * pb_y

                     + 1.5 * fx * fx * pa_y * pb_yy

                     + fx * fx * pb_yyy

                     + 3.0 * pa_xxy * pb_yy * fx

                     + 2.0 * pa_xx * fx * pb_yyy

                     + 0.5 * fx * pa_y * pb_yyyy

                     + pa_xxy * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_yyyz_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_y,
                                    double pb_yyyz,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_z

                     + 0.75 * pa_xx * fx * fx * pb_z

                     + 0.75 * fx * fx * pa_y * pb_yz

                     + 0.75 * fx * fx * pb_yyz

                     + 1.5 * pa_xxy * pb_yz * fx

                     + 1.5 * pa_xx * fx * pb_yyz

                     + 0.5 * fx * pa_y * pb_yyyz

                     + pa_xxy * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_yyzz_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_y,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyzz,
                                    double pb_yzz,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx * pa_y

                     + 0.25 * fx * fx * fx * pb_y

                     + 0.25 * pa_xxy * fx * fx

                     + 0.5 * pa_xx * fx * fx * pb_y

                     + 0.25 * fx * fx * pa_y * pb_yy

                     + 0.25 * fx * fx * pa_y * pb_zz

                     + 0.5 * fx * fx * pb_yzz

                     + 0.5 * pa_xxy * pb_yy * fx

                     + 0.5 * pa_xxy * fx * pb_zz

                     + pa_xx * fx * pb_yzz

                     + 0.5 * fx * pa_y * pb_yyzz

                     + pa_xxy * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_yzzz_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_y,
                                    double pb_yz,
                                    double pb_yzzz,
                                    double pb_z,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_z

                     + 0.75 * pa_xx * fx * fx * pb_z

                     + 0.75 * fx * fx * pa_y * pb_yz

                     + 0.25 * fx * fx * pb_zzz

                     + 1.5 * pa_xxy * pb_yz * fx

                     + 0.5 * pa_xx * fx * pb_zzz

                     + 0.5 * fx * pa_y * pb_yzzz

                     + pa_xxy * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_zzzz_s_0(double fx,
                                    double pa_xxy,
                                    double pa_y,
                                    double pb_zz,
                                    double pb_zzzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_y

                     + 0.75 * pa_xxy * fx * fx

                     + 1.5 * fx * fx * pa_y * pb_zz

                     + 3.0 * pa_xxy * pb_zz * fx

                     + 0.5 * fx * pa_y * pb_zzzz

                     + pa_xxy * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xxxx_s_0(double fx,
                                    double pa_xxz,
                                    double pa_xz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double pb_xxxx,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pa_z

                     + 0.75 * pa_xxz * fx * fx

                     + 6.0 * pa_xz * fx * fx * pb_x

                     + 4.5 * fx * fx * pa_z * pb_xx

                     + 3.0 * pa_xxz * pb_xx * fx

                     + 4.0 * pa_xz * fx * pb_xxx

                     + 0.5 * fx * pa_z * pb_xxxx

                     + pa_xxz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xxxy_s_0(double fx,
                                    double pa_xxz,
                                    double pa_xz,
                                    double pa_z,
                                    double pb_xxxy,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xz * fx * fx * pb_y

                     + 2.25 * fx * fx * pa_z * pb_xy

                     + 1.5 * pa_xxz * pb_xy * fx

                     + 3.0 * pa_xz * fx * pb_xxy

                     + 0.5 * fx * pa_z * pb_xxxy

                     + pa_xxz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xxxz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
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
        return s_0_0 * (0.75 * pa_x * fx * fx * fx

                     + 1.125 * fx * fx * fx * pb_x

                     + 0.75 * pa_xx * fx * fx * pb_x

                     + 1.5 * pa_xz * fx * fx * pb_z

                     + 1.5 * pa_x * fx * fx * pb_xx

                     + 2.25 * fx * fx * pa_z * pb_xz

                     + 0.25 * fx * fx * pb_xxx

                     + 1.5 * pa_xxz * pb_xz * fx

                     + 0.5 * pa_xx * fx * pb_xxx

                     + 3.0 * pa_xz * fx * pb_xxz

                     + 0.5 * fx * pa_z * pb_xxxz

                     + pa_xxz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xxyy_s_0(double fx,
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
        return s_0_0 * (0.375 * fx * fx * fx * pa_z

                     + 0.25 * pa_xxz * fx * fx

                     + pa_xz * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_z * pb_yy

                     + 0.25 * fx * fx * pa_z * pb_xx

                     + 0.5 * pa_xxz * pb_xx * fx

                     + 0.5 * pa_xxz * fx * pb_yy

                     + 2.0 * pa_xz * fx * pb_xyy

                     + 0.5 * fx * pa_z * pb_xxyy

                     + pa_xxz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xxyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
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
        return s_0_0 * (0.375 * fx * fx * fx * pb_y

                     + 0.25 * pa_xx * fx * fx * pb_y

                     + pa_x * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_z * pb_yz

                     + 0.25 * fx * fx * pb_xxy

                     + 0.5 * pa_xxz * fx * pb_yz

                     + 0.5 * pa_xx * fx * pb_xxy

                     + 2.0 * pa_xz * fx * pb_xyz

                     + 0.5 * fx * pa_z * pb_xxyz

                     + pa_xxz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xxzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
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
        return s_0_0 * (0.375 * fx * fx * fx * pa_z

                     + 0.75 * fx * fx * fx * pb_z

                     + 0.25 * pa_xxz * fx * fx

                     + 0.5 * pa_xx * fx * fx * pb_z

                     + pa_xz * fx * fx * pb_x

                     + 2.0 * pa_x * fx * fx * pb_xz

                     + 0.75 * fx * fx * pa_z * pb_zz

                     + 0.25 * fx * fx * pa_z * pb_xx

                     + 0.5 * fx * fx * pb_xxz

                     + 0.5 * pa_xxz * pb_xx * fx

                     + 0.5 * pa_xxz * fx * pb_zz

                     + pa_xx * fx * pb_xxz

                     + 2.0 * pa_xz * fx * pb_xzz

                     + 0.5 * fx * pa_z * pb_xxzz

                     + pa_xxz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xyyy_s_0(double fx,
                                    double pa_xxz,
                                    double pa_xz,
                                    double pa_z,
                                    double pb_xy,
                                    double pb_xyyy,
                                    double pb_y,
                                    double pb_yyy,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xz * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_z * pb_xy

                     + 1.5 * pa_xxz * pb_xy * fx

                     + pa_xz * fx * pb_yyy

                     + 0.5 * fx * pa_z * pb_xyyy

                     + pa_xxz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xyyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
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
        return s_0_0 * (0.25 * pa_x * fx * fx * fx

                     + 0.125 * fx * fx * fx * pb_x

                     + 0.25 * pa_xx * fx * fx * pb_x

                     + 0.5 * pa_xz * fx * fx * pb_z

                     + 0.5 * pa_x * fx * fx * pb_yy

                     + 0.25 * fx * fx * pa_z * pb_xz

                     + 0.25 * fx * fx * pb_xyy

                     + 0.5 * pa_xxz * pb_xz * fx

                     + 0.5 * pa_xx * fx * pb_xyy

                     + pa_xz * fx * pb_yyz

                     + 0.5 * fx * pa_z * pb_xyyz

                     + pa_xxz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xyzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
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
        return s_0_0 * (0.5 * pa_xz * fx * fx * pb_y

                     + pa_x * fx * fx * pb_yz

                     + 0.25 * fx * fx * pa_z * pb_xy

                     + 0.5 * fx * fx * pb_xyz

                     + 0.5 * pa_xxz * pb_xy * fx

                     + pa_xx * fx * pb_xyz

                     + pa_xz * fx * pb_yzz

                     + 0.5 * fx * pa_z * pb_xyzz

                     + pa_xxz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xzzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xx,
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
        return s_0_0 * (0.75 * pa_x * fx * fx * fx

                     + 0.375 * fx * fx * fx * pb_x

                     + 0.75 * pa_xx * fx * fx * pb_x

                     + 1.5 * pa_xz * fx * fx * pb_z

                     + 1.5 * pa_x * fx * fx * pb_zz

                     + 0.75 * fx * fx * pa_z * pb_xz

                     + 0.75 * fx * fx * pb_xzz

                     + 1.5 * pa_xxz * pb_xz * fx

                     + 1.5 * pa_xx * fx * pb_xzz

                     + pa_xz * fx * pb_zzz

                     + 0.5 * fx * pa_z * pb_xzzz

                     + pa_xxz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_yyyy_s_0(double fx,
                                    double pa_xxz,
                                    double pa_z,
                                    double pb_yy,
                                    double pb_yyyy,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_z

                     + 0.75 * pa_xxz * fx * fx

                     + 1.5 * fx * fx * pa_z * pb_yy

                     + 3.0 * pa_xxz * pb_yy * fx

                     + 0.5 * fx * pa_z * pb_yyyy

                     + pa_xxz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_yyyz_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxz,
                                    double pa_z,
                                    double pb_y,
                                    double pb_yyy,
                                    double pb_yyyz,
                                    double pb_yz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_y

                     + 0.75 * pa_xx * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_z * pb_yz

                     + 0.25 * fx * fx * pb_yyy

                     + 1.5 * pa_xxz * pb_yz * fx

                     + 0.5 * pa_xx * fx * pb_yyy

                     + 0.5 * fx * pa_z * pb_yyyz

                     + pa_xxz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_yyzz_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxz,
                                    double pa_z,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_yyzz,
                                    double pb_z,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx * pa_z

                     + 0.25 * fx * fx * fx * pb_z

                     + 0.25 * pa_xxz * fx * fx

                     + 0.5 * pa_xx * fx * fx * pb_z

                     + 0.25 * fx * fx * pa_z * pb_yy

                     + 0.25 * fx * fx * pa_z * pb_zz

                     + 0.5 * fx * fx * pb_yyz

                     + 0.5 * pa_xxz * pb_yy * fx

                     + 0.5 * pa_xxz * fx * pb_zz

                     + pa_xx * fx * pb_yyz

                     + 0.5 * fx * pa_z * pb_yyzz

                     + pa_xxz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_yzzz_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxz,
                                    double pa_z,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double pb_yzzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_y

                     + 0.75 * pa_xx * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_z * pb_yz

                     + 0.75 * fx * fx * pb_yzz

                     + 1.5 * pa_xxz * pb_yz * fx

                     + 1.5 * pa_xx * fx * pb_yzz

                     + 0.5 * fx * pa_z * pb_yzzz

                     + pa_xxz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_zzzz_s_0(double fx,
                                    double pa_xx,
                                    double pa_xxz,
                                    double pa_z,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double pb_zzzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_z

                     + 1.5 * fx * fx * fx * pb_z

                     + 0.75 * pa_xxz * fx * fx

                     + 3.0 * pa_xx * fx * fx * pb_z

                     + 1.5 * fx * fx * pa_z * pb_zz

                     + fx * fx * pb_zzz

                     + 3.0 * pa_xxz * pb_zz * fx

                     + 2.0 * pa_xx * fx * pb_zzz

                     + 0.5 * fx * pa_z * pb_zzzz

                     + pa_xxz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xxxx_s_0(double fx,
                                    double pa_x,
                                    double pa_xyy,
                                    double pa_yy,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double pb_xxxx,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 1.5 * fx * fx * fx * pb_x

                     + 0.75 * pa_xyy * fx * fx

                     + 3.0 * fx * fx * pa_yy * pb_x

                     + 1.5 * pa_x * fx * fx * pb_xx

                     + fx * fx * pb_xxx

                     + 3.0 * pa_xyy * pb_xx * fx

                     + 2.0 * fx * pa_yy * pb_xxx

                     + 0.5 * pa_x * fx * pb_xxxx

                     + pa_xyy * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xxxy_s_0(double fx,
                                    double pa_x,
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
        return s_0_0 * (0.75 * fx * fx * fx * pa_y

                     + 0.375 * fx * fx * fx * pb_y

                     + 1.5 * pa_xy * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_yy * pb_y

                     + 1.5 * fx * fx * pa_y * pb_xx

                     + 0.75 * pa_x * fx * fx * pb_xy

                     + 0.75 * fx * fx * pb_xxy

                     + 1.5 * pa_xyy * pb_xy * fx

                     + pa_xy * fx * pb_xxx

                     + 1.5 * fx * pa_yy * pb_xxy

                     + 0.5 * pa_x * fx * pb_xxxy

                     + pa_xyy * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xxxz_s_0(double fx,
                                    double pa_x,
                                    double pa_xyy,
                                    double pa_yy,
                                    double pb_xxxz,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_z

                     + 0.75 * fx * fx * pa_yy * pb_z

                     + 0.75 * pa_x * fx * fx * pb_xz

                     + 0.75 * fx * fx * pb_xxz

                     + 1.5 * pa_xyy * pb_xz * fx

                     + 1.5 * fx * pa_yy * pb_xxz

                     + 0.5 * pa_x * fx * pb_xxxz

                     + pa_xyy * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xxyy_s_0(double fx,
                                    double pa_x,
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
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 0.75 * fx * fx * fx * pb_x

                     + 0.25 * pa_xyy * fx * fx

                     + pa_xy * fx * fx * pb_y

                     + 0.75 * pa_x * fx * fx * pb_xx

                     + 0.5 * fx * fx * pa_yy * pb_x

                     + 2.0 * fx * fx * pa_y * pb_xy

                     + 0.25 * pa_x * fx * fx * pb_yy

                     + 0.5 * fx * fx * pb_xyy

                     + 0.5 * pa_xyy * pb_xx * fx

                     + 0.5 * pa_xyy * fx * pb_yy

                     + 2.0 * pa_xy * fx * pb_xxy

                     + fx * pa_yy * pb_xyy

                     + 0.5 * pa_x * fx * pb_xxyy

                     + pa_xyy * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xxyz_s_0(double fx,
                                    double pa_x,
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
        return s_0_0 * (0.5 * pa_xy * fx * fx * pb_z

                     + fx * fx * pa_y * pb_xz

                     + 0.25 * pa_x * fx * fx * pb_yz

                     + 0.5 * fx * fx * pb_xyz

                     + 0.5 * pa_xyy * fx * pb_yz

                     + pa_xy * fx * pb_xxz

                     + fx * pa_yy * pb_xyz

                     + 0.5 * pa_x * fx * pb_xxyz

                     + pa_xyy * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xxzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xyy,
                                    double pa_yy,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxzz,
                                    double pb_xzz,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (0.125 * pa_x * fx * fx * fx

                     + 0.25 * fx * fx * fx * pb_x

                     + 0.25 * pa_xyy * fx * fx

                     + 0.5 * fx * fx * pa_yy * pb_x

                     + 0.25 * pa_x * fx * fx * pb_xx

                     + 0.25 * pa_x * fx * fx * pb_zz

                     + 0.5 * fx * fx * pb_xzz

                     + 0.5 * pa_xyy * pb_xx * fx

                     + 0.5 * pa_xyy * fx * pb_zz

                     + fx * pa_yy * pb_xzz

                     + 0.5 * pa_x * fx * pb_xxzz

                     + pa_xyy * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xyyy_s_0(double fx,
                                    double pa_x,
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
        return s_0_0 * (0.75 * fx * fx * fx * pa_y

                     + 1.125 * fx * fx * fx * pb_y

                     + 1.5 * pa_xy * fx * fx * pb_x

                     + 2.25 * pa_x * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_yy * pb_y

                     + 1.5 * fx * fx * pa_y * pb_yy

                     + 0.25 * fx * fx * pb_yyy

                     + 1.5 * pa_xyy * pb_xy * fx

                     + 3.0 * pa_xy * fx * pb_xyy

                     + 0.5 * fx * pa_yy * pb_yyy

                     + 0.5 * pa_x * fx * pb_xyyy

                     + pa_xyy * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xyyz_s_0(double fx,
                                    double pa_x,
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
        return s_0_0 * (0.375 * fx * fx * fx * pb_z

                     + 0.75 * pa_x * fx * fx * pb_xz

                     + 0.25 * fx * fx * pa_yy * pb_z

                     + fx * fx * pa_y * pb_yz

                     + 0.25 * fx * fx * pb_yyz

                     + 0.5 * pa_xyy * pb_xz * fx

                     + 2.0 * pa_xy * fx * pb_xyz

                     + 0.5 * fx * pa_yy * pb_yyz

                     + 0.5 * pa_x * fx * pb_xyyz

                     + pa_xyy * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xyzz_s_0(double fx,
                                    double pa_x,
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
        return s_0_0 * (0.25 * fx * fx * fx * pa_y

                     + 0.125 * fx * fx * fx * pb_y

                     + 0.5 * pa_xy * fx * fx * pb_x

                     + 0.25 * fx * fx * pa_yy * pb_y

                     + 0.5 * fx * fx * pa_y * pb_zz

                     + 0.25 * pa_x * fx * fx * pb_xy

                     + 0.25 * fx * fx * pb_yzz

                     + 0.5 * pa_xyy * pb_xy * fx

                     + pa_xy * fx * pb_xzz

                     + 0.5 * fx * pa_yy * pb_yzz

                     + 0.5 * pa_x * fx * pb_xyzz

                     + pa_xyy * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xzzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xyy,
                                    double pa_yy,
                                    double pb_xz,
                                    double pb_xzzz,
                                    double pb_z,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_z

                     + 0.75 * fx * fx * pa_yy * pb_z

                     + 0.75 * pa_x * fx * fx * pb_xz

                     + 0.25 * fx * fx * pb_zzz

                     + 1.5 * pa_xyy * pb_xz * fx

                     + 0.5 * fx * pa_yy * pb_zzz

                     + 0.5 * pa_x * fx * pb_xzzz

                     + pa_xyy * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_yyyy_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double pb_yyyy,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * pa_x * fx * fx * fx

                     + 0.75 * pa_xyy * fx * fx

                     + 6.0 * pa_xy * fx * fx * pb_y

                     + 4.5 * pa_x * fx * fx * pb_yy

                     + 3.0 * pa_xyy * pb_yy * fx

                     + 4.0 * pa_xy * fx * pb_yyy

                     + 0.5 * pa_x * fx * pb_yyyy

                     + pa_xyy * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_yyyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pb_yyyz,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xy * fx * fx * pb_z

                     + 2.25 * pa_x * fx * fx * pb_yz

                     + 1.5 * pa_xyy * pb_yz * fx

                     + 3.0 * pa_xy * fx * pb_yyz

                     + 0.5 * pa_x * fx * pb_yyyz

                     + pa_xyy * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_yyzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyzz,
                                    double pb_yzz,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 0.25 * pa_xyy * fx * fx

                     + pa_xy * fx * fx * pb_y

                     + 0.75 * pa_x * fx * fx * pb_zz

                     + 0.25 * pa_x * fx * fx * pb_yy

                     + 0.5 * pa_xyy * pb_yy * fx

                     + 0.5 * pa_xyy * fx * pb_zz

                     + 2.0 * pa_xy * fx * pb_yzz

                     + 0.5 * pa_x * fx * pb_yyzz

                     + pa_xyy * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_yzzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pb_yz,
                                    double pb_yzzz,
                                    double pb_z,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xy * fx * fx * pb_z

                     + 0.75 * pa_x * fx * fx * pb_yz

                     + 1.5 * pa_xyy * pb_yz * fx

                     + pa_xy * fx * pb_zzz

                     + 0.5 * pa_x * fx * pb_yzzz

                     + pa_xyy * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_zzzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xyy,
                                    double pb_zz,
                                    double pb_zzzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 0.75 * pa_xyy * fx * fx

                     + 1.5 * pa_x * fx * fx * pb_zz

                     + 3.0 * pa_xyy * pb_zz * fx

                     + 0.5 * pa_x * fx * pb_zzzz

                     + pa_xyy * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xxxx_s_0(double fx,
                                    double pa_xyz,
                                    double pa_yz,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double pb_xxxx,
                                    double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xyz * fx * fx

                     + 3.0 * fx * fx * pa_yz * pb_x

                     + 3.0 * pa_xyz * pb_xx * fx

                     + 2.0 * fx * pa_yz * pb_xxx

                     + pa_xyz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xxxy_s_0(double fx,
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
        return s_0_0 * (0.375 * fx * fx * fx * pa_z

                     + 0.75 * pa_xz * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_yz * pb_y

                     + 0.75 * fx * fx * pa_z * pb_xx

                     + 1.5 * pa_xyz * pb_xy * fx

                     + 0.5 * pa_xz * fx * pb_xxx

                     + 1.5 * fx * pa_yz * pb_xxy

                     + pa_xyz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xxxz_s_0(double fx,
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
        return s_0_0 * (0.375 * fx * fx * fx * pa_y

                     + 0.75 * pa_xy * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_yz * pb_z

                     + 0.75 * fx * fx * pa_y * pb_xx

                     + 1.5 * pa_xyz * pb_xz * fx

                     + 0.5 * pa_xy * fx * pb_xxx

                     + 1.5 * fx * pa_yz * pb_xxz

                     + pa_xyz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xxyy_s_0(double fx,
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
        return s_0_0 * (0.25 * pa_xyz * fx * fx

                     + 0.5 * pa_xz * fx * fx * pb_y

                     + 0.5 * fx * fx * pa_yz * pb_x

                     + fx * fx * pa_z * pb_xy

                     + 0.5 * pa_xyz * pb_xx * fx

                     + 0.5 * pa_xyz * fx * pb_yy

                     + pa_xz * fx * pb_xxy

                     + fx * pa_yz * pb_xyy

                     + pa_xyz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xxyz_s_0(double fx,
                                    double pa_x,
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
        return s_0_0 * (0.125 * pa_x * fx * fx * fx

                     + 0.25 * fx * fx * fx * pb_x

                     + 0.25 * pa_xy * fx * fx * pb_y

                     + 0.25 * pa_xz * fx * fx * pb_z

                     + 0.25 * pa_x * fx * fx * pb_xx

                     + 0.5 * fx * fx * pa_y * pb_xy

                     + 0.5 * fx * fx * pa_z * pb_xz

                     + 0.5 * pa_xyz * fx * pb_yz

                     + 0.5 * pa_xy * fx * pb_xxy

                     + 0.5 * pa_xz * fx * pb_xxz

                     + fx * pa_yz * pb_xyz

                     + pa_xyz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xxzz_s_0(double fx,
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
        return s_0_0 * (0.25 * pa_xyz * fx * fx

                     + 0.5 * pa_xy * fx * fx * pb_z

                     + 0.5 * fx * fx * pa_yz * pb_x

                     + fx * fx * pa_y * pb_xz

                     + 0.5 * pa_xyz * pb_xx * fx

                     + 0.5 * pa_xyz * fx * pb_zz

                     + pa_xy * fx * pb_xxz

                     + fx * pa_yz * pb_xzz

                     + pa_xyz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xyyy_s_0(double fx,
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
        return s_0_0 * (0.375 * fx * fx * fx * pa_z

                     + 0.75 * pa_xz * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_yz * pb_y

                     + 0.75 * fx * fx * pa_z * pb_yy

                     + 1.5 * pa_xyz * pb_xy * fx

                     + 1.5 * pa_xz * fx * pb_xyy

                     + 0.5 * fx * pa_yz * pb_yyy

                     + pa_xyz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xyyz_s_0(double fx,
                                    double pa_x,
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
        return s_0_0 * (0.125 * fx * fx * fx * pa_y

                     + 0.25 * fx * fx * fx * pb_y

                     + 0.25 * pa_xy * fx * fx * pb_x

                     + 0.5 * pa_x * fx * fx * pb_xy

                     + 0.25 * fx * fx * pa_yz * pb_z

                     + 0.25 * fx * fx * pa_y * pb_yy

                     + 0.5 * fx * fx * pa_z * pb_yz

                     + 0.5 * pa_xyz * pb_xz * fx

                     + 0.5 * pa_xy * fx * pb_xyy

                     + pa_xz * fx * pb_xyz

                     + 0.5 * fx * pa_yz * pb_yyz

                     + pa_xyz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xyzz_s_0(double fx,
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
        return s_0_0 * (0.125 * fx * fx * fx * pa_z

                     + 0.25 * fx * fx * fx * pb_z

                     + 0.25 * pa_xz * fx * fx * pb_x

                     + 0.5 * pa_x * fx * fx * pb_xz

                     + 0.25 * fx * fx * pa_yz * pb_y

                     + 0.5 * fx * fx * pa_y * pb_yz

                     + 0.25 * fx * fx * pa_z * pb_zz

                     + 0.5 * pa_xyz * pb_xy * fx

                     + pa_xy * fx * pb_xyz

                     + 0.5 * pa_xz * fx * pb_xzz

                     + 0.5 * fx * pa_yz * pb_yzz

                     + pa_xyz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xzzz_s_0(double fx,
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
        return s_0_0 * (0.375 * fx * fx * fx * pa_y

                     + 0.75 * pa_xy * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_yz * pb_z

                     + 0.75 * fx * fx * pa_y * pb_zz

                     + 1.5 * pa_xyz * pb_xz * fx

                     + 1.5 * pa_xy * fx * pb_xzz

                     + 0.5 * fx * pa_yz * pb_zzz

                     + pa_xyz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_yyyy_s_0(double fx,
                                    double pa_xyz,
                                    double pa_xz,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double pb_yyyy,
                                    double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xyz * fx * fx

                     + 3.0 * pa_xz * fx * fx * pb_y

                     + 3.0 * pa_xyz * pb_yy * fx

                     + 2.0 * pa_xz * fx * pb_yyy

                     + pa_xyz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_yyyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
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
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 0.75 * pa_xy * fx * fx * pb_y

                     + 0.75 * pa_xz * fx * fx * pb_z

                     + 0.75 * pa_x * fx * fx * pb_yy

                     + 1.5 * pa_xyz * pb_yz * fx

                     + 0.5 * pa_xy * fx * pb_yyy

                     + 1.5 * pa_xz * fx * pb_yyz

                     + pa_xyz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_yyzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
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
        return s_0_0 * (0.25 * pa_xyz * fx * fx

                     + 0.5 * pa_xy * fx * fx * pb_z

                     + 0.5 * pa_xz * fx * fx * pb_y

                     + pa_x * fx * fx * pb_yz

                     + 0.5 * pa_xyz * pb_yy * fx

                     + 0.5 * pa_xyz * fx * pb_zz

                     + pa_xy * fx * pb_yyz

                     + pa_xz * fx * pb_yzz

                     + pa_xyz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_yzzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xy,
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
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 0.75 * pa_xy * fx * fx * pb_y

                     + 0.75 * pa_xz * fx * fx * pb_z

                     + 0.75 * pa_x * fx * fx * pb_zz

                     + 1.5 * pa_xyz * pb_yz * fx

                     + 1.5 * pa_xy * fx * pb_yzz

                     + 0.5 * pa_xz * fx * pb_zzz

                     + pa_xyz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_zzzz_s_0(double fx,
                                    double pa_xy,
                                    double pa_xyz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double pb_zzzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xyz * fx * fx

                     + 3.0 * pa_xy * fx * fx * pb_z

                     + 3.0 * pa_xyz * pb_zz * fx

                     + 2.0 * pa_xy * fx * pb_zzz

                     + pa_xyz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xxxx_s_0(double fx,
                                    double pa_x,
                                    double pa_xzz,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double pb_xxxx,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 1.5 * fx * fx * fx * pb_x

                     + 0.75 * pa_xzz * fx * fx

                     + 3.0 * fx * fx * pa_zz * pb_x

                     + 1.5 * pa_x * fx * fx * pb_xx

                     + fx * fx * pb_xxx

                     + 3.0 * pa_xzz * pb_xx * fx

                     + 2.0 * fx * pa_zz * pb_xxx

                     + 0.5 * pa_x * fx * pb_xxxx

                     + pa_xzz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xxxy_s_0(double fx,
                                    double pa_x,
                                    double pa_xzz,
                                    double pa_zz,
                                    double pb_xxxy,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_zz * pb_y

                     + 0.75 * pa_x * fx * fx * pb_xy

                     + 0.75 * fx * fx * pb_xxy

                     + 1.5 * pa_xzz * pb_xy * fx

                     + 1.5 * fx * pa_zz * pb_xxy

                     + 0.5 * pa_x * fx * pb_xxxy

                     + pa_xzz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xxxz_s_0(double fx,
                                    double pa_x,
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
        return s_0_0 * (0.75 * fx * fx * fx * pa_z

                     + 0.375 * fx * fx * fx * pb_z

                     + 1.5 * pa_xz * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_zz * pb_z

                     + 1.5 * fx * fx * pa_z * pb_xx

                     + 0.75 * pa_x * fx * fx * pb_xz

                     + 0.75 * fx * fx * pb_xxz

                     + 1.5 * pa_xzz * pb_xz * fx

                     + pa_xz * fx * pb_xxx

                     + 1.5 * fx * pa_zz * pb_xxz

                     + 0.5 * pa_x * fx * pb_xxxz

                     + pa_xzz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xxyy_s_0(double fx,
                                    double pa_x,
                                    double pa_xzz,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxyy,
                                    double pb_xyy,
                                    double pb_yy,
                                    double s_0_0)
    {
        return s_0_0 * (0.125 * pa_x * fx * fx * fx

                     + 0.25 * fx * fx * fx * pb_x

                     + 0.25 * pa_xzz * fx * fx

                     + 0.5 * fx * fx * pa_zz * pb_x

                     + 0.25 * pa_x * fx * fx * pb_xx

                     + 0.25 * pa_x * fx * fx * pb_yy

                     + 0.5 * fx * fx * pb_xyy

                     + 0.5 * pa_xzz * pb_xx * fx

                     + 0.5 * pa_xzz * fx * pb_yy

                     + fx * pa_zz * pb_xyy

                     + 0.5 * pa_x * fx * pb_xxyy

                     + pa_xzz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xxyz_s_0(double fx,
                                    double pa_x,
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
        return s_0_0 * (0.5 * pa_xz * fx * fx * pb_y

                     + fx * fx * pa_z * pb_xy

                     + 0.25 * pa_x * fx * fx * pb_yz

                     + 0.5 * fx * fx * pb_xyz

                     + 0.5 * pa_xzz * fx * pb_yz

                     + pa_xz * fx * pb_xxy

                     + fx * pa_zz * pb_xyz

                     + 0.5 * pa_x * fx * pb_xxyz

                     + pa_xzz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xxzz_s_0(double fx,
                                    double pa_x,
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
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 0.75 * fx * fx * fx * pb_x

                     + 0.25 * pa_xzz * fx * fx

                     + pa_xz * fx * fx * pb_z

                     + 0.75 * pa_x * fx * fx * pb_xx

                     + 0.5 * fx * fx * pa_zz * pb_x

                     + 2.0 * fx * fx * pa_z * pb_xz

                     + 0.25 * pa_x * fx * fx * pb_zz

                     + 0.5 * fx * fx * pb_xzz

                     + 0.5 * pa_xzz * pb_xx * fx

                     + 0.5 * pa_xzz * fx * pb_zz

                     + 2.0 * pa_xz * fx * pb_xxz

                     + fx * pa_zz * pb_xzz

                     + 0.5 * pa_x * fx * pb_xxzz

                     + pa_xzz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xyyy_s_0(double fx,
                                    double pa_x,
                                    double pa_xzz,
                                    double pa_zz,
                                    double pb_xy,
                                    double pb_xyyy,
                                    double pb_y,
                                    double pb_yyy,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_zz * pb_y

                     + 0.75 * pa_x * fx * fx * pb_xy

                     + 0.25 * fx * fx * pb_yyy

                     + 1.5 * pa_xzz * pb_xy * fx

                     + 0.5 * fx * pa_zz * pb_yyy

                     + 0.5 * pa_x * fx * pb_xyyy

                     + pa_xzz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xyyz_s_0(double fx,
                                    double pa_x,
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
        return s_0_0 * (0.25 * fx * fx * fx * pa_z

                     + 0.125 * fx * fx * fx * pb_z

                     + 0.5 * pa_xz * fx * fx * pb_x

                     + 0.25 * fx * fx * pa_zz * pb_z

                     + 0.5 * fx * fx * pa_z * pb_yy

                     + 0.25 * pa_x * fx * fx * pb_xz

                     + 0.25 * fx * fx * pb_yyz

                     + 0.5 * pa_xzz * pb_xz * fx

                     + pa_xz * fx * pb_xyy

                     + 0.5 * fx * pa_zz * pb_yyz

                     + 0.5 * pa_x * fx * pb_xyyz

                     + pa_xzz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xyzz_s_0(double fx,
                                    double pa_x,
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
        return s_0_0 * (0.375 * fx * fx * fx * pb_y

                     + 0.75 * pa_x * fx * fx * pb_xy

                     + 0.25 * fx * fx * pa_zz * pb_y

                     + fx * fx * pa_z * pb_yz

                     + 0.25 * fx * fx * pb_yzz

                     + 0.5 * pa_xzz * pb_xy * fx

                     + 2.0 * pa_xz * fx * pb_xyz

                     + 0.5 * fx * pa_zz * pb_yzz

                     + 0.5 * pa_x * fx * pb_xyzz

                     + pa_xzz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xzzz_s_0(double fx,
                                    double pa_x,
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
        return s_0_0 * (0.75 * fx * fx * fx * pa_z

                     + 1.125 * fx * fx * fx * pb_z

                     + 1.5 * pa_xz * fx * fx * pb_x

                     + 2.25 * pa_x * fx * fx * pb_xz

                     + 0.75 * fx * fx * pa_zz * pb_z

                     + 1.5 * fx * fx * pa_z * pb_zz

                     + 0.25 * fx * fx * pb_zzz

                     + 1.5 * pa_xzz * pb_xz * fx

                     + 3.0 * pa_xz * fx * pb_xzz

                     + 0.5 * fx * pa_zz * pb_zzz

                     + 0.5 * pa_x * fx * pb_xzzz

                     + pa_xzz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_yyyy_s_0(double fx,
                                    double pa_x,
                                    double pa_xzz,
                                    double pb_yy,
                                    double pb_yyyy,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 0.75 * pa_xzz * fx * fx

                     + 1.5 * pa_x * fx * fx * pb_yy

                     + 3.0 * pa_xzz * pb_yy * fx

                     + 0.5 * pa_x * fx * pb_yyyy

                     + pa_xzz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_yyyz_s_0(double fx,
                                    double pa_x,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pb_y,
                                    double pb_yyy,
                                    double pb_yyyz,
                                    double pb_yz,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xz * fx * fx * pb_y

                     + 0.75 * pa_x * fx * fx * pb_yz

                     + 1.5 * pa_xzz * pb_yz * fx

                     + pa_xz * fx * pb_yyy

                     + 0.5 * pa_x * fx * pb_yyyz

                     + pa_xzz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_yyzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_yyzz,
                                    double pb_z,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_x * fx * fx * fx

                     + 0.25 * pa_xzz * fx * fx

                     + pa_xz * fx * fx * pb_z

                     + 0.75 * pa_x * fx * fx * pb_yy

                     + 0.25 * pa_x * fx * fx * pb_zz

                     + 0.5 * pa_xzz * pb_yy * fx

                     + 0.5 * pa_xzz * fx * pb_zz

                     + 2.0 * pa_xz * fx * pb_yyz

                     + 0.5 * pa_x * fx * pb_yyzz

                     + pa_xzz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_yzzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double pb_yzzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xz * fx * fx * pb_y

                     + 2.25 * pa_x * fx * fx * pb_yz

                     + 1.5 * pa_xzz * pb_yz * fx

                     + 3.0 * pa_xz * fx * pb_yzz

                     + 0.5 * pa_x * fx * pb_yzzz

                     + pa_xzz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_zzzz_s_0(double fx,
                                    double pa_x,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double pb_zzzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * pa_x * fx * fx * fx

                     + 0.75 * pa_xzz * fx * fx

                     + 6.0 * pa_xz * fx * fx * pb_z

                     + 4.5 * pa_x * fx * fx * pb_zz

                     + 3.0 * pa_xzz * pb_zz * fx

                     + 4.0 * pa_xz * fx * pb_zzz

                     + 0.5 * pa_x * fx * pb_zzzz

                     + pa_xzz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xxxx_s_0(double fx,
                                    double pa_y,
                                    double pa_yyy,
                                    double pb_xx,
                                    double pb_xxxx,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * pa_y * fx * fx * fx

                     + 0.75 * pa_yyy * fx * fx

                     + 4.5 * pa_y * fx * fx * pb_xx

                     + 3.0 * pa_yyy * pb_xx * fx

                     + 1.5 * pa_y * fx * pb_xxxx

                     + pa_yyy * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xxxy_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_x,
                                    double pb_xxx,
                                    double pb_xxxy,
                                    double pb_xy,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pb_x

                     + 2.25 * pa_yy * fx * fx * pb_x

                     + 2.25 * pa_y * fx * fx * pb_xy

                     + 0.75 * fx * fx * pb_xxx

                     + 1.5 * pa_yyy * pb_xy * fx

                     + 1.5 * pa_yy * fx * pb_xxx

                     + 1.5 * pa_y * fx * pb_xxxy

                     + pa_yyy * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xxxz_s_0(double fx,
                                    double pa_y,
                                    double pa_yyy,
                                    double pb_xxxz,
                                    double pb_xz,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_y * fx * fx * pb_xz

                     + 1.5 * pa_yyy * pb_xz * fx

                     + 1.5 * pa_y * fx * pb_xxxz

                     + pa_yyy * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xxyy_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_xx,
                                    double pb_xxy,
                                    double pb_xxyy,
                                    double pb_y,
                                    double pb_yy,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * pa_y * fx * fx * fx

                     + 0.75 * fx * fx * fx * pb_y

                     + 0.25 * pa_yyy * fx * fx

                     + 1.5 * pa_yy * fx * fx * pb_y

                     + 2.25 * pa_y * fx * fx * pb_xx

                     + 0.75 * pa_y * fx * fx * pb_yy

                     + 1.5 * fx * fx * pb_xxy

                     + 0.5 * pa_yyy * pb_xx * fx

                     + 0.5 * pa_yyy * fx * pb_yy

                     + 3.0 * pa_yy * fx * pb_xxy

                     + 1.5 * pa_y * fx * pb_xxyy

                     + pa_yyy * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xxyz_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_xxyz,
                                    double pb_xxz,
                                    double pb_yz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_z

                     + 0.75 * pa_yy * fx * fx * pb_z

                     + 0.75 * pa_y * fx * fx * pb_yz

                     + 0.75 * fx * fx * pb_xxz

                     + 0.5 * pa_yyy * fx * pb_yz

                     + 1.5 * pa_yy * fx * pb_xxz

                     + 1.5 * pa_y * fx * pb_xxyz

                     + pa_yyy * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xxzz_s_0(double fx,
                                    double pa_y,
                                    double pa_yyy,
                                    double pb_xx,
                                    double pb_xxzz,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_y * fx * fx * fx

                     + 0.25 * pa_yyy * fx * fx

                     + 0.75 * pa_y * fx * fx * pb_xx

                     + 0.75 * pa_y * fx * fx * pb_zz

                     + 0.5 * pa_yyy * pb_xx * fx

                     + 0.5 * pa_yyy * fx * pb_zz

                     + 1.5 * pa_y * fx * pb_xxzz

                     + pa_yyy * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xyyy_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyy,
                                    double pb_xyyy,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pb_x

                     + 2.25 * pa_yy * fx * fx * pb_x

                     + 6.75 * pa_y * fx * fx * pb_xy

                     + 2.25 * fx * fx * pb_xyy

                     + 1.5 * pa_yyy * pb_xy * fx

                     + 4.5 * pa_yy * fx * pb_xyy

                     + 1.5 * pa_y * fx * pb_xyyy

                     + pa_yyy * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xyyz_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_xyyz,
                                    double pb_xyz,
                                    double pb_xz,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_y * fx * fx * pb_xz

                     + 1.5 * fx * fx * pb_xyz

                     + 0.5 * pa_yyy * pb_xz * fx

                     + 3.0 * pa_yy * fx * pb_xyz

                     + 1.5 * pa_y * fx * pb_xyyz

                     + pa_yyy * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xyzz_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyzz,
                                    double pb_xzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_x

                     + 0.75 * pa_yy * fx * fx * pb_x

                     + 0.75 * pa_y * fx * fx * pb_xy

                     + 0.75 * fx * fx * pb_xzz

                     + 0.5 * pa_yyy * pb_xy * fx

                     + 1.5 * pa_yy * fx * pb_xzz

                     + 1.5 * pa_y * fx * pb_xyzz

                     + pa_yyy * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xzzz_s_0(double fx,
                                    double pa_y,
                                    double pa_yyy,
                                    double pb_xz,
                                    double pb_xzzz,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_y * fx * fx * pb_xz

                     + 1.5 * pa_yyy * pb_xz * fx

                     + 1.5 * pa_y * fx * pb_xzzz

                     + pa_yyy * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_yyyy_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double pb_yyyy,
                                    double s_0_0)
    {
        return s_0_0 * (5.625 * pa_y * fx * fx * fx

                     + 7.5 * fx * fx * fx * pb_y

                     + 0.75 * pa_yyy * fx * fx

                     + 9.0 * pa_yy * fx * fx * pb_y

                     + 13.5 * pa_y * fx * fx * pb_yy

                     + 3.0 * fx * fx * pb_yyy

                     + 3.0 * pa_yyy * pb_yy * fx

                     + 6.0 * pa_yy * fx * pb_yyy

                     + 1.5 * pa_y * fx * pb_yyyy

                     + pa_yyy * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_yyyz_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_yyyz,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pb_z

                     + 2.25 * pa_yy * fx * fx * pb_z

                     + 6.75 * pa_y * fx * fx * pb_yz

                     + 2.25 * fx * fx * pb_yyz

                     + 1.5 * pa_yyy * pb_yz * fx

                     + 4.5 * pa_yy * fx * pb_yyz

                     + 1.5 * pa_y * fx * pb_yyyz

                     + pa_yyy * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_yyzz_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyzz,
                                    double pb_yzz,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * pa_y * fx * fx * fx

                     + 0.75 * fx * fx * fx * pb_y

                     + 0.25 * pa_yyy * fx * fx

                     + 1.5 * pa_yy * fx * fx * pb_y

                     + 2.25 * pa_y * fx * fx * pb_zz

                     + 0.75 * pa_y * fx * fx * pb_yy

                     + 1.5 * fx * fx * pb_yzz

                     + 0.5 * pa_yyy * pb_yy * fx

                     + 0.5 * pa_yyy * fx * pb_zz

                     + 3.0 * pa_yy * fx * pb_yzz

                     + 1.5 * pa_y * fx * pb_yyzz

                     + pa_yyy * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_yzzz_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_yz,
                                    double pb_yzzz,
                                    double pb_z,
                                    double pb_zzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pb_z

                     + 2.25 * pa_yy * fx * fx * pb_z

                     + 2.25 * pa_y * fx * fx * pb_yz

                     + 0.75 * fx * fx * pb_zzz

                     + 1.5 * pa_yyy * pb_yz * fx

                     + 1.5 * pa_yy * fx * pb_zzz

                     + 1.5 * pa_y * fx * pb_yzzz

                     + pa_yyy * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_zzzz_s_0(double fx,
                                    double pa_y,
                                    double pa_yyy,
                                    double pb_zz,
                                    double pb_zzzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * pa_y * fx * fx * fx

                     + 0.75 * pa_yyy * fx * fx

                     + 4.5 * pa_y * fx * fx * pb_zz

                     + 3.0 * pa_yyy * pb_zz * fx

                     + 1.5 * pa_y * fx * pb_zzzz

                     + pa_yyy * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xxxx_s_0(double fx,
                                    double pa_yyz,
                                    double pa_z,
                                    double pb_xx,
                                    double pb_xxxx,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_z

                     + 0.75 * pa_yyz * fx * fx

                     + 1.5 * fx * fx * pa_z * pb_xx

                     + 3.0 * pa_yyz * pb_xx * fx

                     + 0.5 * fx * pa_z * pb_xxxx

                     + pa_yyz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xxxy_s_0(double fx,
                                    double pa_yyz,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xxx,
                                    double pb_xxxy,
                                    double pb_xy,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_yz * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_z * pb_xy

                     + 1.5 * pa_yyz * pb_xy * fx

                     + pa_yz * fx * pb_xxx

                     + 0.5 * fx * pa_z * pb_xxxy

                     + pa_yyz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xxxz_s_0(double fx,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xxx,
                                    double pb_xxxz,
                                    double pb_xz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_x

                     + 0.75 * pa_yy * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_z * pb_xz

                     + 0.25 * fx * fx * pb_xxx

                     + 1.5 * pa_yyz * pb_xz * fx

                     + 0.5 * pa_yy * fx * pb_xxx

                     + 0.5 * fx * pa_z * pb_xxxz

                     + pa_yyz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xxyy_s_0(double fx,
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
        return s_0_0 * (0.375 * fx * fx * fx * pa_z

                     + 0.25 * pa_yyz * fx * fx

                     + pa_yz * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_z * pb_xx

                     + 0.25 * fx * fx * pa_z * pb_yy

                     + 0.5 * pa_yyz * pb_xx * fx

                     + 0.5 * pa_yyz * fx * pb_yy

                     + 2.0 * pa_yz * fx * pb_xxy

                     + 0.5 * fx * pa_z * pb_xxyy

                     + pa_yyz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xxyz_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
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
        return s_0_0 * (0.25 * pa_y * fx * fx * fx

                     + 0.125 * fx * fx * fx * pb_y

                     + 0.25 * pa_yy * fx * fx * pb_y

                     + 0.5 * pa_yz * fx * fx * pb_z

                     + 0.5 * pa_y * fx * fx * pb_xx

                     + 0.25 * fx * fx * pa_z * pb_yz

                     + 0.25 * fx * fx * pb_xxy

                     + 0.5 * pa_yyz * fx * pb_yz

                     + 0.5 * pa_yy * fx * pb_xxy

                     + pa_yz * fx * pb_xxz

                     + 0.5 * fx * pa_z * pb_xxyz

                     + pa_yyz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xxzz_s_0(double fx,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_z,
                                    double pb_xx,
                                    double pb_xxz,
                                    double pb_xxzz,
                                    double pb_z,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx * pa_z

                     + 0.25 * fx * fx * fx * pb_z

                     + 0.25 * pa_yyz * fx * fx

                     + 0.5 * pa_yy * fx * fx * pb_z

                     + 0.25 * fx * fx * pa_z * pb_xx

                     + 0.25 * fx * fx * pa_z * pb_zz

                     + 0.5 * fx * fx * pb_xxz

                     + 0.5 * pa_yyz * pb_xx * fx

                     + 0.5 * pa_yyz * fx * pb_zz

                     + pa_yy * fx * pb_xxz

                     + 0.5 * fx * pa_z * pb_xxzz

                     + pa_yyz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xyyy_s_0(double fx,
                                    double pa_yyz,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyy,
                                    double pb_xyyy,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_yz * fx * fx * pb_x

                     + 2.25 * fx * fx * pa_z * pb_xy

                     + 1.5 * pa_yyz * pb_xy * fx

                     + 3.0 * pa_yz * fx * pb_xyy

                     + 0.5 * fx * pa_z * pb_xyyy

                     + pa_yyz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xyyz_s_0(double fx,
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
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_x

                     + 0.25 * pa_yy * fx * fx * pb_x

                     + pa_y * fx * fx * pb_xy

                     + 0.75 * fx * fx * pa_z * pb_xz

                     + 0.25 * fx * fx * pb_xyy

                     + 0.5 * pa_yyz * pb_xz * fx

                     + 0.5 * pa_yy * fx * pb_xyy

                     + 2.0 * pa_yz * fx * pb_xyz

                     + 0.5 * fx * pa_z * pb_xyyz

                     + pa_yyz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xyzz_s_0(double fx,
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
                                    double s_0_0)
    {
        return s_0_0 * (0.5 * pa_yz * fx * fx * pb_x

                     + pa_y * fx * fx * pb_xz

                     + 0.25 * fx * fx * pa_z * pb_xy

                     + 0.5 * fx * fx * pb_xyz

                     + 0.5 * pa_yyz * pb_xy * fx

                     + pa_yy * fx * pb_xyz

                     + pa_yz * fx * pb_xzz

                     + 0.5 * fx * pa_z * pb_xyzz

                     + pa_yyz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xzzz_s_0(double fx,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xz,
                                    double pb_xzz,
                                    double pb_xzzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_x

                     + 0.75 * pa_yy * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_z * pb_xz

                     + 0.75 * fx * fx * pb_xzz

                     + 1.5 * pa_yyz * pb_xz * fx

                     + 1.5 * pa_yy * fx * pb_xzz

                     + 0.5 * fx * pa_z * pb_xzzz

                     + pa_yyz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_yyyy_s_0(double fx,
                                    double pa_yyz,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double pb_yyyy,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pa_z

                     + 0.75 * pa_yyz * fx * fx

                     + 6.0 * pa_yz * fx * fx * pb_y

                     + 4.5 * fx * fx * pa_z * pb_yy

                     + 3.0 * pa_yyz * pb_yy * fx

                     + 4.0 * pa_yz * fx * pb_yyy

                     + 0.5 * fx * pa_z * pb_yyyy

                     + pa_yyz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_yyyz_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
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
        return s_0_0 * (0.75 * pa_y * fx * fx * fx

                     + 1.125 * fx * fx * fx * pb_y

                     + 0.75 * pa_yy * fx * fx * pb_y

                     + 1.5 * pa_yz * fx * fx * pb_z

                     + 1.5 * pa_y * fx * fx * pb_yy

                     + 2.25 * fx * fx * pa_z * pb_yz

                     + 0.25 * fx * fx * pb_yyy

                     + 1.5 * pa_yyz * pb_yz * fx

                     + 0.5 * pa_yy * fx * pb_yyy

                     + 3.0 * pa_yz * fx * pb_yyz

                     + 0.5 * fx * pa_z * pb_yyyz

                     + pa_yyz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_yyzz_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
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
        return s_0_0 * (0.375 * fx * fx * fx * pa_z

                     + 0.75 * fx * fx * fx * pb_z

                     + 0.25 * pa_yyz * fx * fx

                     + 0.5 * pa_yy * fx * fx * pb_z

                     + pa_yz * fx * fx * pb_y

                     + 2.0 * pa_y * fx * fx * pb_yz

                     + 0.75 * fx * fx * pa_z * pb_zz

                     + 0.25 * fx * fx * pa_z * pb_yy

                     + 0.5 * fx * fx * pb_yyz

                     + 0.5 * pa_yyz * pb_yy * fx

                     + 0.5 * pa_yyz * fx * pb_zz

                     + pa_yy * fx * pb_yyz

                     + 2.0 * pa_yz * fx * pb_yzz

                     + 0.5 * fx * pa_z * pb_yyzz

                     + pa_yyz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_yzzz_s_0(double fx,
                                    double pa_y,
                                    double pa_yy,
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
        return s_0_0 * (0.75 * pa_y * fx * fx * fx

                     + 0.375 * fx * fx * fx * pb_y

                     + 0.75 * pa_yy * fx * fx * pb_y

                     + 1.5 * pa_yz * fx * fx * pb_z

                     + 1.5 * pa_y * fx * fx * pb_zz

                     + 0.75 * fx * fx * pa_z * pb_yz

                     + 0.75 * fx * fx * pb_yzz

                     + 1.5 * pa_yyz * pb_yz * fx

                     + 1.5 * pa_yy * fx * pb_yzz

                     + pa_yz * fx * pb_zzz

                     + 0.5 * fx * pa_z * pb_yzzz

                     + pa_yyz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_zzzz_s_0(double fx,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_z,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double pb_zzzz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pa_z

                     + 1.5 * fx * fx * fx * pb_z

                     + 0.75 * pa_yyz * fx * fx

                     + 3.0 * pa_yy * fx * fx * pb_z

                     + 1.5 * fx * fx * pa_z * pb_zz

                     + fx * fx * pb_zzz

                     + 3.0 * pa_yyz * pb_zz * fx

                     + 2.0 * pa_yy * fx * pb_zzz

                     + 0.5 * fx * pa_z * pb_zzzz

                     + pa_yyz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xxxx_s_0(double fx,
                                    double pa_y,
                                    double pa_yzz,
                                    double pb_xx,
                                    double pb_xxxx,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_y * fx * fx * fx

                     + 0.75 * pa_yzz * fx * fx

                     + 1.5 * pa_y * fx * fx * pb_xx

                     + 3.0 * pa_yzz * pb_xx * fx

                     + 0.5 * pa_y * fx * pb_xxxx

                     + pa_yzz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xxxy_s_0(double fx,
                                    double pa_y,
                                    double pa_yzz,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xxx,
                                    double pb_xxxy,
                                    double pb_xy,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_zz * pb_x

                     + 0.75 * pa_y * fx * fx * pb_xy

                     + 0.25 * fx * fx * pb_xxx

                     + 1.5 * pa_yzz * pb_xy * fx

                     + 0.5 * fx * pa_zz * pb_xxx

                     + 0.5 * pa_y * fx * pb_xxxy

                     + pa_yzz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xxxz_s_0(double fx,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pb_x,
                                    double pb_xxx,
                                    double pb_xxxz,
                                    double pb_xz,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_yz * fx * fx * pb_x

                     + 0.75 * pa_y * fx * fx * pb_xz

                     + 1.5 * pa_yzz * pb_xz * fx

                     + pa_yz * fx * pb_xxx

                     + 0.5 * pa_y * fx * pb_xxxz

                     + pa_yzz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xxyy_s_0(double fx,
                                    double pa_y,
                                    double pa_yzz,
                                    double pa_zz,
                                    double pb_xx,
                                    double pb_xxy,
                                    double pb_xxyy,
                                    double pb_y,
                                    double pb_yy,
                                    double s_0_0)
    {
        return s_0_0 * (0.125 * pa_y * fx * fx * fx

                     + 0.25 * fx * fx * fx * pb_y

                     + 0.25 * pa_yzz * fx * fx

                     + 0.5 * fx * fx * pa_zz * pb_y

                     + 0.25 * pa_y * fx * fx * pb_xx

                     + 0.25 * pa_y * fx * fx * pb_yy

                     + 0.5 * fx * fx * pb_xxy

                     + 0.5 * pa_yzz * pb_xx * fx

                     + 0.5 * pa_yzz * fx * pb_yy

                     + fx * pa_zz * pb_xxy

                     + 0.5 * pa_y * fx * pb_xxyy

                     + pa_yzz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xxyz_s_0(double fx,
                                    double pa_y,
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
        return s_0_0 * (0.25 * fx * fx * fx * pa_z

                     + 0.125 * fx * fx * fx * pb_z

                     + 0.5 * pa_yz * fx * fx * pb_y

                     + 0.25 * fx * fx * pa_zz * pb_z

                     + 0.5 * fx * fx * pa_z * pb_xx

                     + 0.25 * pa_y * fx * fx * pb_yz

                     + 0.25 * fx * fx * pb_xxz

                     + 0.5 * pa_yzz * fx * pb_yz

                     + pa_yz * fx * pb_xxy

                     + 0.5 * fx * pa_zz * pb_xxz

                     + 0.5 * pa_y * fx * pb_xxyz

                     + pa_yzz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xxzz_s_0(double fx,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pb_xx,
                                    double pb_xxz,
                                    double pb_xxzz,
                                    double pb_z,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_y * fx * fx * fx

                     + 0.25 * pa_yzz * fx * fx

                     + pa_yz * fx * fx * pb_z

                     + 0.75 * pa_y * fx * fx * pb_xx

                     + 0.25 * pa_y * fx * fx * pb_zz

                     + 0.5 * pa_yzz * pb_xx * fx

                     + 0.5 * pa_yzz * fx * pb_zz

                     + 2.0 * pa_yz * fx * pb_xxz

                     + 0.5 * pa_y * fx * pb_xxzz

                     + pa_yzz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xyyy_s_0(double fx,
                                    double pa_y,
                                    double pa_yzz,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyy,
                                    double pb_xyyy,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_zz * pb_x

                     + 0.75 * pa_y * fx * fx * pb_xy

                     + 0.75 * fx * fx * pb_xyy

                     + 1.5 * pa_yzz * pb_xy * fx

                     + 1.5 * fx * pa_zz * pb_xyy

                     + 0.5 * pa_y * fx * pb_xyyy

                     + pa_yzz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xyyz_s_0(double fx,
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
                                    double s_0_0)
    {
        return s_0_0 * (0.5 * pa_yz * fx * fx * pb_x

                     + fx * fx * pa_z * pb_xy

                     + 0.25 * pa_y * fx * fx * pb_xz

                     + 0.5 * fx * fx * pb_xyz

                     + 0.5 * pa_yzz * pb_xz * fx

                     + pa_yz * fx * pb_xyy

                     + fx * pa_zz * pb_xyz

                     + 0.5 * pa_y * fx * pb_xyyz

                     + pa_yzz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xyzz_s_0(double fx,
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
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_x

                     + 0.75 * pa_y * fx * fx * pb_xy

                     + 0.25 * fx * fx * pa_zz * pb_x

                     + fx * fx * pa_z * pb_xz

                     + 0.25 * fx * fx * pb_xzz

                     + 0.5 * pa_yzz * pb_xy * fx

                     + 2.0 * pa_yz * fx * pb_xyz

                     + 0.5 * fx * pa_zz * pb_xzz

                     + 0.5 * pa_y * fx * pb_xyzz

                     + pa_yzz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xzzz_s_0(double fx,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pb_x,
                                    double pb_xz,
                                    double pb_xzz,
                                    double pb_xzzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.5 * pa_yz * fx * fx * pb_x

                     + 2.25 * pa_y * fx * fx * pb_xz

                     + 1.5 * pa_yzz * pb_xz * fx

                     + 3.0 * pa_yz * fx * pb_xzz

                     + 0.5 * pa_y * fx * pb_xzzz

                     + pa_yzz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_yyyy_s_0(double fx,
                                    double pa_y,
                                    double pa_yzz,
                                    double pa_zz,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double pb_yyyy,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_y * fx * fx * fx

                     + 1.5 * fx * fx * fx * pb_y

                     + 0.75 * pa_yzz * fx * fx

                     + 3.0 * fx * fx * pa_zz * pb_y

                     + 1.5 * pa_y * fx * fx * pb_yy

                     + fx * fx * pb_yyy

                     + 3.0 * pa_yzz * pb_yy * fx

                     + 2.0 * fx * pa_zz * pb_yyy

                     + 0.5 * pa_y * fx * pb_yyyy

                     + pa_yzz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_yyyz_s_0(double fx,
                                    double pa_y,
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
        return s_0_0 * (0.75 * fx * fx * fx * pa_z

                     + 0.375 * fx * fx * fx * pb_z

                     + 1.5 * pa_yz * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_zz * pb_z

                     + 1.5 * fx * fx * pa_z * pb_yy

                     + 0.75 * pa_y * fx * fx * pb_yz

                     + 0.75 * fx * fx * pb_yyz

                     + 1.5 * pa_yzz * pb_yz * fx

                     + pa_yz * fx * pb_yyy

                     + 1.5 * fx * pa_zz * pb_yyz

                     + 0.5 * pa_y * fx * pb_yyyz

                     + pa_yzz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_yyzz_s_0(double fx,
                                    double pa_y,
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
        return s_0_0 * (0.375 * pa_y * fx * fx * fx

                     + 0.75 * fx * fx * fx * pb_y

                     + 0.25 * pa_yzz * fx * fx

                     + pa_yz * fx * fx * pb_z

                     + 0.75 * pa_y * fx * fx * pb_yy

                     + 0.5 * fx * fx * pa_zz * pb_y

                     + 2.0 * fx * fx * pa_z * pb_yz

                     + 0.25 * pa_y * fx * fx * pb_zz

                     + 0.5 * fx * fx * pb_yzz

                     + 0.5 * pa_yzz * pb_yy * fx

                     + 0.5 * pa_yzz * fx * pb_zz

                     + 2.0 * pa_yz * fx * pb_yyz

                     + fx * pa_zz * pb_yzz

                     + 0.5 * pa_y * fx * pb_yyzz

                     + pa_yzz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_yzzz_s_0(double fx,
                                    double pa_y,
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
        return s_0_0 * (0.75 * fx * fx * fx * pa_z

                     + 1.125 * fx * fx * fx * pb_z

                     + 1.5 * pa_yz * fx * fx * pb_y

                     + 2.25 * pa_y * fx * fx * pb_yz

                     + 0.75 * fx * fx * pa_zz * pb_z

                     + 1.5 * fx * fx * pa_z * pb_zz

                     + 0.25 * fx * fx * pb_zzz

                     + 1.5 * pa_yzz * pb_yz * fx

                     + 3.0 * pa_yz * fx * pb_yzz

                     + 0.5 * fx * pa_zz * pb_zzz

                     + 0.5 * pa_y * fx * pb_yzzz

                     + pa_yzz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_zzzz_s_0(double fx,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double pb_zzzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * pa_y * fx * fx * fx

                     + 0.75 * pa_yzz * fx * fx

                     + 6.0 * pa_yz * fx * fx * pb_z

                     + 4.5 * pa_y * fx * fx * pb_zz

                     + 3.0 * pa_yzz * pb_zz * fx

                     + 4.0 * pa_yz * fx * pb_zzz

                     + 0.5 * pa_y * fx * pb_zzzz

                     + pa_yzz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xxxx_s_0(double fx,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_xx,
                                    double pb_xxxx,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * pa_z * fx * fx * fx

                     + 0.75 * pa_zzz * fx * fx

                     + 4.5 * pa_z * fx * fx * pb_xx

                     + 3.0 * pa_zzz * pb_xx * fx

                     + 1.5 * pa_z * fx * pb_xxxx

                     + pa_zzz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xxxy_s_0(double fx,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_xxxy,
                                    double pb_xy,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_z * fx * fx * pb_xy

                     + 1.5 * pa_zzz * pb_xy * fx

                     + 1.5 * pa_z * fx * pb_xxxy

                     + pa_zzz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xxxz_s_0(double fx,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_x,
                                    double pb_xxx,
                                    double pb_xxxz,
                                    double pb_xz,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pb_x

                     + 2.25 * pa_zz * fx * fx * pb_x

                     + 2.25 * pa_z * fx * fx * pb_xz

                     + 0.75 * fx * fx * pb_xxx

                     + 1.5 * pa_zzz * pb_xz * fx

                     + 1.5 * pa_zz * fx * pb_xxx

                     + 1.5 * pa_z * fx * pb_xxxz

                     + pa_zzz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xxyy_s_0(double fx,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_xx,
                                    double pb_xxyy,
                                    double pb_yy,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * pa_z * fx * fx * fx

                     + 0.25 * pa_zzz * fx * fx

                     + 0.75 * pa_z * fx * fx * pb_xx

                     + 0.75 * pa_z * fx * fx * pb_yy

                     + 0.5 * pa_zzz * pb_xx * fx

                     + 0.5 * pa_zzz * fx * pb_yy

                     + 1.5 * pa_z * fx * pb_xxyy

                     + pa_zzz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xxyz_s_0(double fx,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_xxy,
                                    double pb_xxyz,
                                    double pb_y,
                                    double pb_yz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_y

                     + 0.75 * pa_zz * fx * fx * pb_y

                     + 0.75 * pa_z * fx * fx * pb_yz

                     + 0.75 * fx * fx * pb_xxy

                     + 0.5 * pa_zzz * fx * pb_yz

                     + 1.5 * pa_zz * fx * pb_xxy

                     + 1.5 * pa_z * fx * pb_xxyz

                     + pa_zzz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xxzz_s_0(double fx,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_xx,
                                    double pb_xxz,
                                    double pb_xxzz,
                                    double pb_z,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * pa_z * fx * fx * fx

                     + 0.75 * fx * fx * fx * pb_z

                     + 0.25 * pa_zzz * fx * fx

                     + 1.5 * pa_zz * fx * fx * pb_z

                     + 2.25 * pa_z * fx * fx * pb_xx

                     + 0.75 * pa_z * fx * fx * pb_zz

                     + 1.5 * fx * fx * pb_xxz

                     + 0.5 * pa_zzz * pb_xx * fx

                     + 0.5 * pa_zzz * fx * pb_zz

                     + 3.0 * pa_zz * fx * pb_xxz

                     + 1.5 * pa_z * fx * pb_xxzz

                     + pa_zzz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xyyy_s_0(double fx,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_xy,
                                    double pb_xyyy,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_z * fx * fx * pb_xy

                     + 1.5 * pa_zzz * pb_xy * fx

                     + 1.5 * pa_z * fx * pb_xyyy

                     + pa_zzz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xyyz_s_0(double fx,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_x,
                                    double pb_xyy,
                                    double pb_xyyz,
                                    double pb_xz,
                                    double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx * pb_x

                     + 0.75 * pa_zz * fx * fx * pb_x

                     + 0.75 * pa_z * fx * fx * pb_xz

                     + 0.75 * fx * fx * pb_xyy

                     + 0.5 * pa_zzz * pb_xz * fx

                     + 1.5 * pa_zz * fx * pb_xyy

                     + 1.5 * pa_z * fx * pb_xyyz

                     + pa_zzz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xyzz_s_0(double fx,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_xy,
                                    double pb_xyz,
                                    double pb_xyzz,
                                    double s_0_0)
    {
        return s_0_0 * (2.25 * pa_z * fx * fx * pb_xy

                     + 1.5 * fx * fx * pb_xyz

                     + 0.5 * pa_zzz * pb_xy * fx

                     + 3.0 * pa_zz * fx * pb_xyz

                     + 1.5 * pa_z * fx * pb_xyzz

                     + pa_zzz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xzzz_s_0(double fx,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_x,
                                    double pb_xz,
                                    double pb_xzz,
                                    double pb_xzzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pb_x

                     + 2.25 * pa_zz * fx * fx * pb_x

                     + 6.75 * pa_z * fx * fx * pb_xz

                     + 2.25 * fx * fx * pb_xzz

                     + 1.5 * pa_zzz * pb_xz * fx

                     + 4.5 * pa_zz * fx * pb_xzz

                     + 1.5 * pa_z * fx * pb_xzzz

                     + pa_zzz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_yyyy_s_0(double fx,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_yy,
                                    double pb_yyyy,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * pa_z * fx * fx * fx

                     + 0.75 * pa_zzz * fx * fx

                     + 4.5 * pa_z * fx * fx * pb_yy

                     + 3.0 * pa_zzz * pb_yy * fx

                     + 1.5 * pa_z * fx * pb_yyyy

                     + pa_zzz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_yyyz_s_0(double fx,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_y,
                                    double pb_yyy,
                                    double pb_yyyz,
                                    double pb_yz,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * fx * fx * fx * pb_y

                     + 2.25 * pa_zz * fx * fx * pb_y

                     + 2.25 * pa_z * fx * fx * pb_yz

                     + 0.75 * fx * fx * pb_yyy

                     + 1.5 * pa_zzz * pb_yz * fx

                     + 1.5 * pa_zz * fx * pb_yyy

                     + 1.5 * pa_z * fx * pb_yyyz

                     + pa_zzz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_yyzz_s_0(double fx,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_yyzz,
                                    double pb_z,
                                    double pb_zz,
                                    double s_0_0)
    {
        return s_0_0 * (1.125 * pa_z * fx * fx * fx

                     + 0.75 * fx * fx * fx * pb_z

                     + 0.25 * pa_zzz * fx * fx

                     + 1.5 * pa_zz * fx * fx * pb_z

                     + 2.25 * pa_z * fx * fx * pb_yy

                     + 0.75 * pa_z * fx * fx * pb_zz

                     + 1.5 * fx * fx * pb_yyz

                     + 0.5 * pa_zzz * pb_yy * fx

                     + 0.5 * pa_zzz * fx * pb_zz

                     + 3.0 * pa_zz * fx * pb_yyz

                     + 1.5 * pa_z * fx * pb_yyzz

                     + pa_zzz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_yzzz_s_0(double fx,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double pb_yzzz,
                                    double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx * pb_y

                     + 2.25 * pa_zz * fx * fx * pb_y

                     + 6.75 * pa_z * fx * fx * pb_yz

                     + 2.25 * fx * fx * pb_yzz

                     + 1.5 * pa_zzz * pb_yz * fx

                     + 4.5 * pa_zz * fx * pb_yzz

                     + 1.5 * pa_z * fx * pb_yzzz

                     + pa_zzz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_zzzz_s_0(double fx,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double pb_zzzz,
                                    double s_0_0)
    {
        return s_0_0 * (5.625 * pa_z * fx * fx * fx

                     + 7.5 * fx * fx * fx * pb_z

                     + 0.75 * pa_zzz * fx * fx

                     + 9.0 * pa_zz * fx * fx * pb_z

                     + 13.5 * pa_z * fx * fx * pb_zz

                     + 3.0 * fx * fx * pb_zzz

                     + 3.0 * pa_zzz * pb_zz * fx

                     + 6.0 * pa_zz * fx * pb_zzz

                     + 1.5 * pa_z * fx * pb_zzzz

                     + pa_zzz * pb_zzzz);

    }


} // ovlrecfunc namespace

