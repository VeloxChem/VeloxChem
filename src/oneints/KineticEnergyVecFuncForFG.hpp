//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

namespace kinvecfunc { // kinvecfunc namespace

    // SIMD elementary functions for (F|T|G) integrals

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
    inline double fvec_xxx_xxxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double pb_xxxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 13.5 * pa_x * fx * fx * fz * fgb

                     + 45.0 * pa_x * fx * fx * fx * fz

                     + 60.0 * fx * fx * fx * fz * pb_x

                     - 2.25 * pa_x * fz * fga * fx * fx

                     - 9.0 * fx * fx * pb_x * fz * fgb

                     - 9.0 * fx * fx * fz * fga * pb_x

                     - 3.0 * pa_xxx * fx * fz * fgb

                     - 18.0 * pa_xx * fx * pb_x * fz * fgb

                     + 7.5 * pa_xxx * fz * fx * fx

                     + 90.0 * pa_xx * fz * fx * fx * pb_x

                     + 135.0 * pa_x * fx * fx * fz * pb_xx

                     - 9.0 * pa_x * fx * pb_xx * fz * fgb

                     - 9.0 * pa_x * fz * fga * pb_xx * fx

                     - 6.0 * fx * fz * fga * pb_xxx

                     - 6.0 * pa_xxx * pb_xx * fz * fgb

                     + 30.0 * fx * fx * fz * pb_xxx

                     + 36.0 * pa_xxx * fz * pb_xx * fx

                     + 72.0 * pa_xx * fz * fx * pb_xxx

                     - 3.0 * pa_x * fz * fga * pb_xxxx

                     + 18.0 * pa_x * fz * fx * pb_xxxx

                     + 14.0 * pa_xxx * fz * pb_xxxx);

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
    inline double fvec_xxx_xxxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pb_xxxy,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double r_0_0)
    {
        return r_0_0 * (15.0 * fx * fx * fx * fz * pb_y

                     - 2.25 * fx * fx * fz * fgb * pb_y

                     - 2.25 * fx * fx * fz * fga * pb_y

                     - 4.5 * pa_xx * fx * fz * fgb * pb_y

                     + 22.5 * pa_xx * fz * fx * fx * pb_y

                     + 67.5 * pa_x * fx * fx * fz * pb_xy

                     - 4.5 * pa_x * fx * pb_xy * fz * fgb

                     - 4.5 * pa_x * fz * fga * pb_xy * fx

                     - 4.5 * fx * fz * fga * pb_xxy

                     - 3.0 * pa_xxx * pb_xy * fz * fgb

                     + 22.5 * fx * fx * fz * pb_xxy

                     + 18.0 * pa_xxx * fz * pb_xy * fx

                     + 54.0 * pa_xx * fz * fx * pb_xxy

                     - 3.0 * pa_x * fz * fga * pb_xxxy

                     + 18.0 * pa_x * fz * fx * pb_xxxy

                     + 14.0 * pa_xxx * fz * pb_xxxy);

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
    inline double fvec_xxx_xxxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pb_xxxz,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (15.0 * fx * fx * fx * fz * pb_z

                     - 2.25 * fx * fx * fz * fgb * pb_z

                     - 2.25 * fx * fx * fz * fga * pb_z

                     - 4.5 * pa_xx * fx * fz * fgb * pb_z

                     + 22.5 * pa_xx * fz * fx * fx * pb_z

                     + 67.5 * pa_x * fx * fx * fz * pb_xz

                     - 4.5 * pa_x * fx * pb_xz * fz * fgb

                     - 4.5 * pa_x * fz * fga * pb_xz * fx

                     - 4.5 * fx * fz * fga * pb_xxz

                     - 3.0 * pa_xxx * pb_xz * fz * fgb

                     + 22.5 * fx * fx * fz * pb_xxz

                     + 18.0 * pa_xxx * fz * pb_xz * fx

                     + 54.0 * pa_xx * fz * fx * pb_xxz

                     - 3.0 * pa_x * fz * fga * pb_xxxz

                     + 18.0 * pa_x * fz * fx * pb_xxxz

                     + 14.0 * pa_xxx * fz * pb_xxxz);

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
    inline double fvec_xxx_xxyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxyy,
                                    double pb_xyy,
                                    double pb_yy,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_x * fx * fx * fz * fgb

                     + 9.0 * pa_x * fx * fx * fx * fz

                     - 0.75 * pa_x * fz * fga * fx * fx

                     - 1.5 * fx * fx * pb_x * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_x

                     - pa_xxx * fx * fz * fgb

                     - 3.0 * pa_xx * fx * pb_x * fz * fgb

                     + 6.0 * fx * fx * fx * fz * pb_x

                     + 2.5 * pa_xxx * fz * fx * fx

                     + 15.0 * pa_xx * fz * fx * fx * pb_x

                     + 22.5 * pa_x * fx * fx * fz * pb_yy

                     - 1.5 * pa_x * fx * pb_xx * fz * fgb

                     - 1.5 * pa_x * fx * fz * fgb * pb_yy

                     - 1.5 * pa_x * fz * fga * pb_xx * fx

                     - 1.5 * pa_x * fz * fga * fx * pb_yy

                     - 3.0 * fx * fz * fga * pb_xyy

                     - pa_xxx * pb_xx * fz * fgb

                     - pa_xxx * fz * fgb * pb_yy

                     + 7.5 * pa_x * fz * fx * fx * pb_xx

                     + 15.0 * fx * fx * fz * pb_xyy

                     + 6.0 * pa_xxx * fz * pb_xx * fx

                     + 6.0 * pa_xxx * fz * fx * pb_yy

                     + 36.0 * pa_xx * fz * fx * pb_xyy

                     - 3.0 * pa_x * fz * fga * pb_xxyy

                     + 18.0 * pa_x * fz * fx * pb_xxyy

                     + 14.0 * pa_xxx * fz * pb_xxyy);

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
    inline double fvec_xxx_xxyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pb_xxyz,
                                    double pb_xyz,
                                    double pb_yz,
                                    double r_0_0)
    {
        return r_0_0 * (22.5 * pa_x * fx * fx * fz * pb_yz

                     - 1.5 * pa_x * fx * fz * fgb * pb_yz

                     - 1.5 * pa_x * fz * fga * fx * pb_yz

                     - 3.0 * fx * fz * fga * pb_xyz

                     - pa_xxx * fz * fgb * pb_yz

                     + 15.0 * fx * fx * fz * pb_xyz

                     + 6.0 * pa_xxx * fz * fx * pb_yz

                     + 36.0 * pa_xx * fz * fx * pb_xyz

                     - 3.0 * pa_x * fz * fga * pb_xxyz

                     + 18.0 * pa_x * fz * fx * pb_xxyz

                     + 14.0 * pa_xxx * fz * pb_xxyz);

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
    inline double fvec_xxx_xxzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxzz,
                                    double pb_xzz,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_x * fx * fx * fz * fgb

                     + 9.0 * pa_x * fx * fx * fx * fz

                     - 0.75 * pa_x * fz * fga * fx * fx

                     - 1.5 * fx * fx * pb_x * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_x

                     - pa_xxx * fx * fz * fgb

                     - 3.0 * pa_xx * fx * pb_x * fz * fgb

                     + 6.0 * fx * fx * fx * fz * pb_x

                     + 2.5 * pa_xxx * fz * fx * fx

                     + 15.0 * pa_xx * fz * fx * fx * pb_x

                     + 22.5 * pa_x * fx * fx * fz * pb_zz

                     - 1.5 * pa_x * fx * pb_xx * fz * fgb

                     - 1.5 * pa_x * fx * fz * fgb * pb_zz

                     - 1.5 * pa_x * fz * fga * pb_xx * fx

                     - 1.5 * pa_x * fz * fga * fx * pb_zz

                     - 3.0 * fx * fz * fga * pb_xzz

                     - pa_xxx * pb_xx * fz * fgb

                     - pa_xxx * fz * fgb * pb_zz

                     + 7.5 * pa_x * fz * fx * fx * pb_xx

                     + 15.0 * fx * fx * fz * pb_xzz

                     + 6.0 * pa_xxx * fz * pb_xx * fx

                     + 6.0 * pa_xxx * fz * fx * pb_zz

                     + 36.0 * pa_xx * fz * fx * pb_xzz

                     - 3.0 * pa_x * fz * fga * pb_xxzz

                     + 18.0 * pa_x * fz * fx * pb_xxzz

                     + 14.0 * pa_xxx * fz * pb_xxzz);

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
    inline double fvec_xxx_xyyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pb_xy,
                                    double pb_xyyy,
                                    double pb_y,
                                    double pb_yyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pb_y * fz * fgb

                     - 2.25 * fx * fx * fz * fga * pb_y

                     - 4.5 * pa_xx * fx * pb_y * fz * fgb

                     + 9.0 * fx * fx * fx * fz * pb_y

                     + 22.5 * pa_xx * fz * fx * fx * pb_y

                     - 4.5 * pa_x * fx * pb_xy * fz * fgb

                     - 4.5 * pa_x * fz * fga * pb_xy * fx

                     - 1.5 * fx * fz * fga * pb_yyy

                     - 3.0 * pa_xxx * pb_xy * fz * fgb

                     + 22.5 * pa_x * fz * fx * fx * pb_xy

                     + 7.5 * fx * fx * fz * pb_yyy

                     + 18.0 * pa_xxx * fz * pb_xy * fx

                     + 18.0 * pa_xx * fz * fx * pb_yyy

                     - 3.0 * pa_x * fz * fga * pb_xyyy

                     + 18.0 * pa_x * fz * fx * pb_xyyy

                     + 14.0 * pa_xxx * fz * pb_xyyy);

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
    inline double fvec_xxx_xyyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pb_xyyz,
                                    double pb_xz,
                                    double pb_yyz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb * pb_z

                     - 0.75 * fx * fx * fz * fga * pb_z

                     - 1.5 * pa_xx * fx * fz * fgb * pb_z

                     + 3.0 * fx * fx * fx * fz * pb_z

                     + 7.5 * pa_xx * fz * fx * fx * pb_z

                     - 1.5 * pa_x * fx * pb_xz * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_xz * fx

                     - 1.5 * fx * fz * fga * pb_yyz

                     - pa_xxx * pb_xz * fz * fgb

                     + 7.5 * pa_x * fz * fx * fx * pb_xz

                     + 7.5 * fx * fx * fz * pb_yyz

                     + 6.0 * pa_xxx * fz * pb_xz * fx

                     + 18.0 * pa_xx * fz * fx * pb_yyz

                     - 3.0 * pa_x * fz * fga * pb_xyyz

                     + 18.0 * pa_x * fz * fx * pb_xyyz

                     + 14.0 * pa_xxx * fz * pb_xyyz);

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
    inline double fvec_xxx_xyzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pb_xy,
                                    double pb_xyzz,
                                    double pb_y,
                                    double pb_yzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pb_y * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pb_y

                     - 1.5 * pa_xx * fx * pb_y * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_y

                     + 7.5 * pa_xx * fz * fx * fx * pb_y

                     - 1.5 * pa_x * fx * pb_xy * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_xy * fx

                     - 1.5 * fx * fz * fga * pb_yzz

                     - pa_xxx * pb_xy * fz * fgb

                     + 7.5 * pa_x * fz * fx * fx * pb_xy

                     + 7.5 * fx * fx * fz * pb_yzz

                     + 6.0 * pa_xxx * fz * pb_xy * fx

                     + 18.0 * pa_xx * fz * fx * pb_yzz

                     - 3.0 * pa_x * fz * fga * pb_xyzz

                     + 18.0 * pa_x * fz * fx * pb_xyzz

                     + 14.0 * pa_xxx * fz * pb_xyzz);

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
    inline double fvec_xxx_xzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xx,
                                    double pa_xxx,
                                    double pb_xz,
                                    double pb_xzzz,
                                    double pb_z,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pb_z * fz * fgb

                     - 2.25 * fx * fx * fz * fga * pb_z

                     - 4.5 * pa_xx * fx * pb_z * fz * fgb

                     + 9.0 * fx * fx * fx * fz * pb_z

                     + 22.5 * pa_xx * fz * fx * fx * pb_z

                     - 4.5 * pa_x * fx * pb_xz * fz * fgb

                     - 4.5 * pa_x * fz * fga * pb_xz * fx

                     - 1.5 * fx * fz * fga * pb_zzz

                     - 3.0 * pa_xxx * pb_xz * fz * fgb

                     + 22.5 * pa_x * fz * fx * fx * pb_xz

                     + 7.5 * fx * fx * fz * pb_zzz

                     + 18.0 * pa_xxx * fz * pb_xz * fx

                     + 18.0 * pa_xx * fz * fx * pb_zzz

                     - 3.0 * pa_x * fz * fga * pb_xzzz

                     + 18.0 * pa_x * fz * fx * pb_xzzz

                     + 14.0 * pa_xxx * fz * pb_xzzz);

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
    inline double fvec_xxx_yyyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xxx,
                                    double pb_yy,
                                    double pb_yyyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_x * fx * fx * fz * fgb

                     - 2.25 * pa_x * fz * fga * fx * fx

                     - 3.0 * pa_xxx * fx * fz * fgb

                     + 9.0 * pa_x * fz * fx * fx * fx

                     + 7.5 * pa_xxx * fz * fx * fx

                     - 9.0 * pa_x * fx * pb_yy * fz * fgb

                     - 9.0 * pa_x * fz * fga * pb_yy * fx

                     - 6.0 * pa_xxx * pb_yy * fz * fgb

                     + 45.0 * pa_x * fz * fx * fx * pb_yy

                     + 36.0 * pa_xxx * fz * pb_yy * fx

                     - 3.0 * pa_x * fz * fga * pb_yyyy

                     + 18.0 * pa_x * fz * fx * pb_yyyy

                     + 14.0 * pa_xxx * fz * pb_yyyy);

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
    inline double fvec_xxx_yyyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xxx,
                                    double pb_yyyz,
                                    double pb_yz,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_x * fx * pb_yz * fz * fgb

                     - 4.5 * pa_x * fz * fga * pb_yz * fx

                     - 3.0 * pa_xxx * pb_yz * fz * fgb

                     + 22.5 * pa_x * fz * fx * fx * pb_yz

                     + 18.0 * pa_xxx * fz * pb_yz * fx

                     - 3.0 * pa_x * fz * fga * pb_yyyz

                     + 18.0 * pa_x * fz * fx * pb_yyyz

                     + 14.0 * pa_xxx * fz * pb_yyyz);

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
    inline double fvec_xxx_yyzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xxx,
                                    double pb_yy,
                                    double pb_yyzz,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fx * fz * fgb

                     - 0.75 * pa_x * fz * fga * fx * fx

                     - pa_xxx * fx * fz * fgb

                     + 3.0 * pa_x * fz * fx * fx * fx

                     + 2.5 * pa_xxx * fz * fx * fx

                     - 1.5 * pa_x * fx * pb_yy * fz * fgb

                     - 1.5 * pa_x * fx * fz * fgb * pb_zz

                     - 1.5 * pa_x * fz * fga * pb_yy * fx

                     - 1.5 * pa_x * fz * fga * fx * pb_zz

                     - pa_xxx * pb_yy * fz * fgb

                     - pa_xxx * fz * fgb * pb_zz

                     + 7.5 * pa_x * fz * fx * fx * pb_yy

                     + 7.5 * pa_x * fz * fx * fx * pb_zz

                     + 6.0 * pa_xxx * fz * pb_yy * fx

                     + 6.0 * pa_xxx * fz * fx * pb_zz

                     - 3.0 * pa_x * fz * fga * pb_yyzz

                     + 18.0 * pa_x * fz * fx * pb_yyzz

                     + 14.0 * pa_xxx * fz * pb_yyzz);

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
    inline double fvec_xxx_yzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xxx,
                                    double pb_yz,
                                    double pb_yzzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_x * fx * pb_yz * fz * fgb

                     - 4.5 * pa_x * fz * fga * pb_yz * fx

                     - 3.0 * pa_xxx * pb_yz * fz * fgb

                     + 22.5 * pa_x * fz * fx * fx * pb_yz

                     + 18.0 * pa_xxx * fz * pb_yz * fx

                     - 3.0 * pa_x * fz * fga * pb_yzzz

                     + 18.0 * pa_x * fz * fx * pb_yzzz

                     + 14.0 * pa_xxx * fz * pb_yzzz);

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
    inline double fvec_xxx_zzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xxx,
                                    double pb_zz,
                                    double pb_zzzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_x * fx * fx * fz * fgb

                     - 2.25 * pa_x * fz * fga * fx * fx

                     - 3.0 * pa_xxx * fx * fz * fgb

                     + 9.0 * pa_x * fz * fx * fx * fx

                     + 7.5 * pa_xxx * fz * fx * fx

                     - 9.0 * pa_x * fx * pb_zz * fz * fgb

                     - 9.0 * pa_x * fz * fga * pb_zz * fx

                     - 6.0 * pa_xxx * pb_zz * fz * fgb

                     + 45.0 * pa_x * fz * fx * fx * pb_zz

                     + 36.0 * pa_xxx * fz * pb_zz * fx

                     - 3.0 * pa_x * fz * fga * pb_zzzz

                     + 18.0 * pa_x * fz * fx * pb_zzzz

                     + 14.0 * pa_xxx * fz * pb_zzzz);

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
    inline double fvec_xxy_xxxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxy,
                                    double pa_xy,
                                    double pa_y,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double pb_xxxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * fx * fx * pa_y * fz * fgb

                     + 15.0 * fx * fx * fx * fz * pa_y

                     - 0.75 * fz * fga * pa_y * fx * fx

                     - 3.0 * pa_xxy * fx * fz * fgb

                     - 12.0 * pa_xy * fx * pb_x * fz * fgb

                     + 7.5 * pa_xxy * fz * fx * fx

                     + 60.0 * pa_xy * fx * fx * fz * pb_x

                     + 45.0 * fx * fx * fz * pa_y * pb_xx

                     - 3.0 * fx * pa_y * pb_xx * fz * fgb

                     - 3.0 * fz * fga * pa_y * pb_xx * fx

                     - 6.0 * pa_xxy * pb_xx * fz * fgb

                     + 36.0 * pa_xxy * fz * pb_xx * fx

                     + 48.0 * pa_xy * fx * fz * pb_xxx

                     - fz * fga * pa_y * pb_xxxx

                     + 6.0 * fx * fz * pa_y * pb_xxxx

                     + 14.0 * pa_xxy * fz * pb_xxxx);

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
    inline double fvec_xxy_xxxy_r_0(double fga,
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
                                    double pb_xxx,
                                    double pb_xxxy,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fx * fz * fgb

                     + 6.0 * pa_x * fx * fx * fx * fz

                     + 9.0 * fx * fx * fx * fz * pb_x

                     - 0.75 * fx * fx * pb_x * fz * fgb

                     - 0.75 * fz * fga * fx * fx * pb_x

                     - 1.5 * pa_xx * fx * pb_x * fz * fgb

                     - 3.0 * pa_xy * fx * fz * fgb * pb_y

                     + 7.5 * pa_xx * fz * fx * fx * pb_x

                     + 15.0 * pa_xy * fx * fx * fz * pb_y

                     + 15.0 * pa_x * fx * fx * fz * pb_xx

                     + 22.5 * fx * fx * fz * pa_y * pb_xy

                     - 1.5 * fx * pa_y * pb_xy * fz * fgb

                     - 1.5 * fz * fga * pa_y * pb_xy * fx

                     - 0.5 * fz * fga * fx * pb_xxx

                     - 3.0 * pa_xxy * pb_xy * fz * fgb

                     + 2.5 * fx * fx * fz * pb_xxx

                     + 18.0 * pa_xxy * fz * pb_xy * fx

                     + 6.0 * pa_xx * fz * fx * pb_xxx

                     + 36.0 * pa_xy * fx * fz * pb_xxy

                     - fz * fga * pa_y * pb_xxxy

                     + 6.0 * fx * fz * pa_y * pb_xxxy

                     + 14.0 * pa_xxy * fz * pb_xxxy);

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
    inline double fvec_xxy_xxxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxy,
                                    double pa_xy,
                                    double pa_y,
                                    double pb_xxxz,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xy * fx * fz * fgb * pb_z

                     + 15.0 * pa_xy * fx * fx * fz * pb_z

                     + 22.5 * fx * fx * fz * pa_y * pb_xz

                     - 1.5 * fx * pa_y * pb_xz * fz * fgb

                     - 1.5 * fz * fga * pa_y * pb_xz * fx

                     - 3.0 * pa_xxy * pb_xz * fz * fgb

                     + 18.0 * pa_xxy * fz * pb_xz * fx

                     + 36.0 * pa_xy * fx * fz * pb_xxz

                     - fz * fga * pa_y * pb_xxxz

                     + 6.0 * fx * fz * pa_y * pb_xxxz

                     + 14.0 * pa_xxy * fz * pb_xxxz);

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
    inline double fvec_xxy_xxyy_r_0(double fga,
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
                                    double pb_xxyy,
                                    double pb_xy,
                                    double pb_xyy,
                                    double pb_y,
                                    double pb_yy,
                                    double r_0_0)
    {
        return r_0_0 * (- fx * fx * pa_y * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_y

                     + 6.0 * fx * fx * fx * fz * pb_y

                     - 0.5 * fx * fx * fz * fgb * pb_y

                     - 0.25 * fz * fga * pa_y * fx * fx

                     - 0.5 * fz * fga * fx * fx * pb_y

                     - pa_xxy * fx * fz * fgb

                     - pa_xx * fx * fz * fgb * pb_y

                     - 2.0 * pa_xy * fx * pb_x * fz * fgb

                     + 2.5 * pa_xxy * fz * fx * fx

                     + 5.0 * pa_xx * fz * fx * fx * pb_y

                     + 10.0 * pa_xy * fx * fx * fz * pb_x

                     + 20.0 * pa_x * fx * fx * fz * pb_xy

                     + 7.5 * fx * fx * fz * pa_y * pb_yy

                     - 0.5 * fx * pa_y * pb_xx * fz * fgb

                     - 0.5 * fx * pa_y * fz * fgb * pb_yy

                     - 0.5 * fz * fga * pa_y * pb_xx * fx

                     - 0.5 * fz * fga * pa_y * fx * pb_yy

                     - fz * fga * fx * pb_xxy

                     - pa_xxy * pb_xx * fz * fgb

                     - pa_xxy * fz * fgb * pb_yy

                     + 2.5 * fx * fx * fz * pa_y * pb_xx

                     + 5.0 * fx * fx * fz * pb_xxy

                     + 6.0 * pa_xxy * fz * pb_xx * fx

                     + 6.0 * pa_xxy * fz * fx * pb_yy

                     + 12.0 * pa_xx * fz * fx * pb_xxy

                     + 24.0 * pa_xy * fx * fz * pb_xyy

                     - fz * fga * pa_y * pb_xxyy

                     + 6.0 * fx * fz * pa_y * pb_xxyy

                     + 14.0 * pa_xxy * fz * pb_xxyy);

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
    inline double fvec_xxy_xxyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
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
                                    double r_0_0)
    {
        return r_0_0 * (3.0 * fx * fx * fx * fz * pb_z

                     - 0.25 * fx * fx * fz * fgb * pb_z

                     - 0.25 * fz * fga * fx * fx * pb_z

                     - 0.5 * pa_xx * fx * fz * fgb * pb_z

                     + 2.5 * pa_xx * fz * fx * fx * pb_z

                     + 10.0 * pa_x * fx * fx * fz * pb_xz

                     + 7.5 * fx * fx * fz * pa_y * pb_yz

                     - 0.5 * fx * pa_y * fz * fgb * pb_yz

                     - 0.5 * fz * fga * pa_y * fx * pb_yz

                     - 0.5 * fz * fga * fx * pb_xxz

                     - pa_xxy * fz * fgb * pb_yz

                     + 2.5 * fx * fx * fz * pb_xxz

                     + 6.0 * pa_xxy * fz * fx * pb_yz

                     + 6.0 * pa_xx * fz * fx * pb_xxz

                     + 24.0 * pa_xy * fx * fz * pb_xyz

                     - fz * fga * pa_y * pb_xxyz

                     + 6.0 * fx * fz * pa_y * pb_xxyz

                     + 14.0 * pa_xxy * fz * pb_xxyz);

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
    inline double fvec_xxy_xxzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
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
        return r_0_0 * (- fx * fx * pa_y * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_y

                     - 0.25 * fz * fga * pa_y * fx * fx

                     - pa_xxy * fx * fz * fgb

                     - 2.0 * pa_xy * fx * pb_x * fz * fgb

                     + 2.5 * pa_xxy * fz * fx * fx

                     + 10.0 * pa_xy * fx * fx * fz * pb_x

                     + 7.5 * fx * fx * fz * pa_y * pb_zz

                     - 0.5 * fx * pa_y * pb_xx * fz * fgb

                     - 0.5 * fx * pa_y * fz * fgb * pb_zz

                     - 0.5 * fz * fga * pa_y * pb_xx * fx

                     - 0.5 * fz * fga * pa_y * fx * pb_zz

                     - pa_xxy * pb_xx * fz * fgb

                     - pa_xxy * fz * fgb * pb_zz

                     + 2.5 * fx * fx * fz * pa_y * pb_xx

                     + 6.0 * pa_xxy * fz * pb_xx * fx

                     + 6.0 * pa_xxy * fz * fx * pb_zz

                     + 24.0 * pa_xy * fx * fz * pb_xzz

                     - fz * fga * pa_y * pb_xxzz

                     + 6.0 * fx * fz * pa_y * pb_xxzz

                     + 14.0 * pa_xxy * fz * pb_xxzz);

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
    inline double fvec_xxy_xyyy_r_0(double fga,
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
                                    double pb_xyyy,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fx * fz * fgb

                     + 6.0 * pa_x * fx * fx * fx * fz

                     - 0.75 * fx * fx * pb_x * fz * fgb

                     - 0.75 * fz * fga * fx * fx * pb_x

                     - 1.5 * pa_xx * fx * pb_x * fz * fgb

                     - 3.0 * pa_xy * fx * pb_y * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_x

                     + 7.5 * pa_xx * fz * fx * fx * pb_x

                     + 15.0 * pa_xy * fx * fx * fz * pb_y

                     + 15.0 * pa_x * fx * fx * fz * pb_yy

                     - 1.5 * fx * pa_y * pb_xy * fz * fgb

                     - 1.5 * fz * fga * pa_y * pb_xy * fx

                     - 1.5 * fz * fga * fx * pb_xyy

                     - 3.0 * pa_xxy * pb_xy * fz * fgb

                     + 7.5 * fx * fx * fz * pa_y * pb_xy

                     + 7.5 * fx * fx * fz * pb_xyy

                     + 18.0 * pa_xxy * fz * pb_xy * fx

                     + 18.0 * pa_xx * fz * fx * pb_xyy

                     + 12.0 * pa_xy * fx * fz * pb_yyy

                     - fz * fga * pa_y * pb_xyyy

                     + 6.0 * fx * fz * pa_y * pb_xyyy

                     + 14.0 * pa_xxy * fz * pb_xyyy);

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
    inline double fvec_xxy_xyyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
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
                                    double r_0_0)
    {
        return r_0_0 * (- pa_xy * fx * fz * fgb * pb_z

                     + 5.0 * pa_xy * fx * fx * fz * pb_z

                     + 10.0 * pa_x * fx * fx * fz * pb_yz

                     - 0.5 * fx * pa_y * pb_xz * fz * fgb

                     - 0.5 * fz * fga * pa_y * pb_xz * fx

                     - fz * fga * fx * pb_xyz

                     - pa_xxy * pb_xz * fz * fgb

                     + 2.5 * fx * fx * fz * pa_y * pb_xz

                     + 5.0 * fx * fx * fz * pb_xyz

                     + 6.0 * pa_xxy * fz * pb_xz * fx

                     + 12.0 * pa_xx * fz * fx * pb_xyz

                     + 12.0 * pa_xy * fx * fz * pb_yyz

                     - fz * fga * pa_y * pb_xyyz

                     + 6.0 * fx * fz * pa_y * pb_xyyz

                     + 14.0 * pa_xxy * fz * pb_xyyz);

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
    inline double fvec_xxy_xyzz_r_0(double fga,
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
                                    double pb_xyzz,
                                    double pb_xzz,
                                    double pb_y,
                                    double pb_yzz,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_x * fx * fx * fz * fgb

                     + 2.0 * pa_x * fx * fx * fx * fz

                     - 0.25 * fx * fx * pb_x * fz * fgb

                     - 0.25 * fz * fga * fx * fx * pb_x

                     - 0.5 * pa_xx * fx * pb_x * fz * fgb

                     - pa_xy * fx * pb_y * fz * fgb

                     + fx * fx * fx * fz * pb_x

                     + 2.5 * pa_xx * fz * fx * fx * pb_x

                     + 5.0 * pa_xy * fx * fx * fz * pb_y

                     + 5.0 * pa_x * fx * fx * fz * pb_zz

                     - 0.5 * fx * pa_y * pb_xy * fz * fgb

                     - 0.5 * fz * fga * pa_y * pb_xy * fx

                     - 0.5 * fz * fga * fx * pb_xzz

                     - pa_xxy * pb_xy * fz * fgb

                     + 2.5 * fx * fx * fz * pa_y * pb_xy

                     + 2.5 * fx * fx * fz * pb_xzz

                     + 6.0 * pa_xxy * fz * pb_xy * fx

                     + 6.0 * pa_xx * fz * fx * pb_xzz

                     + 12.0 * pa_xy * fx * fz * pb_yzz

                     - fz * fga * pa_y * pb_xyzz

                     + 6.0 * fx * fz * pa_y * pb_xyzz

                     + 14.0 * pa_xxy * fz * pb_xyzz);

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
    inline double fvec_xxy_xzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxy,
                                    double pa_xy,
                                    double pa_y,
                                    double pb_xz,
                                    double pb_xzzz,
                                    double pb_z,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xy * fx * pb_z * fz * fgb

                     + 15.0 * pa_xy * fx * fx * fz * pb_z

                     - 1.5 * fx * pa_y * pb_xz * fz * fgb

                     - 1.5 * fz * fga * pa_y * pb_xz * fx

                     - 3.0 * pa_xxy * pb_xz * fz * fgb

                     + 7.5 * fx * fx * fz * pa_y * pb_xz

                     + 18.0 * pa_xxy * fz * pb_xz * fx

                     + 12.0 * pa_xy * fx * fz * pb_zzz

                     - fz * fga * pa_y * pb_xzzz

                     + 6.0 * fx * fz * pa_y * pb_xzzz

                     + 14.0 * pa_xxy * fz * pb_xzzz);

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
    inline double fvec_xxy_yyyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_y,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double pb_yyyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_y * fz * fgb

                     - 3.0 * fx * fx * pb_y * fz * fgb

                     - 0.75 * fz * fga * pa_y * fx * fx

                     - 3.0 * fz * fga * fx * fx * pb_y

                     - 3.0 * pa_xxy * fx * fz * fgb

                     - 6.0 * pa_xx * fx * pb_y * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_y

                     + 12.0 * fx * fx * fx * fz * pb_y

                     + 7.5 * pa_xxy * fz * fx * fx

                     + 30.0 * pa_xx * fz * fx * fx * pb_y

                     - 3.0 * fx * pa_y * pb_yy * fz * fgb

                     - 3.0 * fz * fga * pa_y * pb_yy * fx

                     - 2.0 * fz * fga * fx * pb_yyy

                     - 6.0 * pa_xxy * pb_yy * fz * fgb

                     + 15.0 * fx * fx * fz * pa_y * pb_yy

                     + 10.0 * fx * fx * fz * pb_yyy

                     + 36.0 * pa_xxy * fz * pb_yy * fx

                     + 24.0 * pa_xx * fz * fx * pb_yyy

                     - fz * fga * pa_y * pb_yyyy

                     + 6.0 * fx * fz * pa_y * pb_yyyy

                     + 14.0 * pa_xxy * fz * pb_yyyy);

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
    inline double fvec_xxy_yyyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_y,
                                    double pb_yyyz,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb * pb_z

                     - 0.75 * fz * fga * fx * fx * pb_z

                     - 1.5 * pa_xx * fx * fz * fgb * pb_z

                     + 3.0 * fx * fx * fx * fz * pb_z

                     + 7.5 * pa_xx * fz * fx * fx * pb_z

                     - 1.5 * fx * pa_y * pb_yz * fz * fgb

                     - 1.5 * fz * fga * pa_y * pb_yz * fx

                     - 1.5 * fz * fga * fx * pb_yyz

                     - 3.0 * pa_xxy * pb_yz * fz * fgb

                     + 7.5 * fx * fx * fz * pa_y * pb_yz

                     + 7.5 * fx * fx * fz * pb_yyz

                     + 18.0 * pa_xxy * fz * pb_yz * fx

                     + 18.0 * pa_xx * fz * fx * pb_yyz

                     - fz * fga * pa_y * pb_yyyz

                     + 6.0 * fx * fz * pa_y * pb_yyyz

                     + 14.0 * pa_xxy * fz * pb_yyyz);

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
    inline double fvec_xxy_yyzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_y,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyzz,
                                    double pb_yzz,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * pa_y * fz * fgb

                     - 0.5 * fx * fx * pb_y * fz * fgb

                     - 0.25 * fz * fga * pa_y * fx * fx

                     - 0.5 * fz * fga * fx * fx * pb_y

                     - pa_xxy * fx * fz * fgb

                     - pa_xx * fx * pb_y * fz * fgb

                     + fx * fx * fx * fz * pa_y

                     + 2.0 * fx * fx * fx * fz * pb_y

                     + 2.5 * pa_xxy * fz * fx * fx

                     + 5.0 * pa_xx * fz * fx * fx * pb_y

                     - 0.5 * fx * pa_y * pb_yy * fz * fgb

                     - 0.5 * fx * pa_y * fz * fgb * pb_zz

                     - 0.5 * fz * fga * pa_y * pb_yy * fx

                     - 0.5 * fz * fga * pa_y * fx * pb_zz

                     - fz * fga * fx * pb_yzz

                     - pa_xxy * pb_yy * fz * fgb

                     - pa_xxy * fz * fgb * pb_zz

                     + 2.5 * fx * fx * fz * pa_y * pb_yy

                     + 2.5 * fx * fx * fz * pa_y * pb_zz

                     + 5.0 * fx * fx * fz * pb_yzz

                     + 6.0 * pa_xxy * fz * pb_yy * fx

                     + 6.0 * pa_xxy * fz * fx * pb_zz

                     + 12.0 * pa_xx * fz * fx * pb_yzz

                     - fz * fga * pa_y * pb_yyzz

                     + 6.0 * fx * fz * pa_y * pb_yyzz

                     + 14.0 * pa_xxy * fz * pb_yyzz);

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
    inline double fvec_xxy_yzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxy,
                                    double pa_y,
                                    double pb_yz,
                                    double pb_yzzz,
                                    double pb_z,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pb_z * fz * fgb

                     - 0.75 * fz * fga * fx * fx * pb_z

                     - 1.5 * pa_xx * fx * pb_z * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_z

                     + 7.5 * pa_xx * fz * fx * fx * pb_z

                     - 1.5 * fx * pa_y * pb_yz * fz * fgb

                     - 1.5 * fz * fga * pa_y * pb_yz * fx

                     - 0.5 * fz * fga * fx * pb_zzz

                     - 3.0 * pa_xxy * pb_yz * fz * fgb

                     + 7.5 * fx * fx * fz * pa_y * pb_yz

                     + 2.5 * fx * fx * fz * pb_zzz

                     + 18.0 * pa_xxy * fz * pb_yz * fx

                     + 6.0 * pa_xx * fz * fx * pb_zzz

                     - fz * fga * pa_y * pb_yzzz

                     + 6.0 * fx * fz * pa_y * pb_yzzz

                     + 14.0 * pa_xxy * fz * pb_yzzz);

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
    inline double fvec_xxy_zzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxy,
                                    double pa_y,
                                    double pb_zz,
                                    double pb_zzzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_y * fz * fgb

                     - 0.75 * fz * fga * pa_y * fx * fx

                     - 3.0 * pa_xxy * fx * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_y

                     + 7.5 * pa_xxy * fz * fx * fx

                     - 3.0 * fx * pa_y * pb_zz * fz * fgb

                     - 3.0 * fz * fga * pa_y * pb_zz * fx

                     - 6.0 * pa_xxy * pb_zz * fz * fgb

                     + 15.0 * fx * fx * fz * pa_y * pb_zz

                     + 36.0 * pa_xxy * fz * pb_zz * fx

                     - fz * fga * pa_y * pb_zzzz

                     + 6.0 * fx * fz * pa_y * pb_zzzz

                     + 14.0 * pa_xxy * fz * pb_zzzz);

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
    inline double fvec_xxz_xxxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxz,
                                    double pa_xz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double pb_xxxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * fx * fx * pa_z * fz * fgb

                     + 15.0 * fx * fx * fx * fz * pa_z

                     - 0.75 * fz * fga * pa_z * fx * fx

                     - 3.0 * pa_xxz * fx * fz * fgb

                     - 12.0 * pa_xz * fx * pb_x * fz * fgb

                     + 7.5 * pa_xxz * fz * fx * fx

                     + 60.0 * pa_xz * fx * fx * fz * pb_x

                     + 45.0 * fx * fx * fz * pa_z * pb_xx

                     - 3.0 * fx * pa_z * pb_xx * fz * fgb

                     - 3.0 * fz * fga * pa_z * pb_xx * fx

                     - 6.0 * pa_xxz * pb_xx * fz * fgb

                     + 36.0 * pa_xxz * fz * pb_xx * fx

                     + 48.0 * pa_xz * fx * fz * pb_xxx

                     - fz * fga * pa_z * pb_xxxx

                     + 6.0 * fx * fz * pa_z * pb_xxxx

                     + 14.0 * pa_xxz * fz * pb_xxxx);

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
    inline double fvec_xxz_xxxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxz,
                                    double pa_xz,
                                    double pa_z,
                                    double pb_xxxy,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xz * fx * fz * fgb * pb_y

                     + 15.0 * pa_xz * fx * fx * fz * pb_y

                     + 22.5 * fx * fx * fz * pa_z * pb_xy

                     - 1.5 * fx * pa_z * pb_xy * fz * fgb

                     - 1.5 * fz * fga * pa_z * pb_xy * fx

                     - 3.0 * pa_xxz * pb_xy * fz * fgb

                     + 18.0 * pa_xxz * fz * pb_xy * fx

                     + 36.0 * pa_xz * fx * fz * pb_xxy

                     - fz * fga * pa_z * pb_xxxy

                     + 6.0 * fx * fz * pa_z * pb_xxxy

                     + 14.0 * pa_xxz * fz * pb_xxxy);

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
    inline double fvec_xxz_xxxz_r_0(double fga,
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
                                    double pb_xxx,
                                    double pb_xxxz,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fx * fz * fgb

                     + 6.0 * pa_x * fx * fx * fx * fz

                     + 9.0 * fx * fx * fx * fz * pb_x

                     - 0.75 * fx * fx * pb_x * fz * fgb

                     - 0.75 * fz * fga * fx * fx * pb_x

                     - 1.5 * pa_xx * fx * pb_x * fz * fgb

                     - 3.0 * pa_xz * fx * fz * fgb * pb_z

                     + 7.5 * pa_xx * fz * fx * fx * pb_x

                     + 15.0 * pa_xz * fx * fx * fz * pb_z

                     + 15.0 * pa_x * fx * fx * fz * pb_xx

                     + 22.5 * fx * fx * fz * pa_z * pb_xz

                     - 1.5 * fx * pa_z * pb_xz * fz * fgb

                     - 1.5 * fz * fga * pa_z * pb_xz * fx

                     - 0.5 * fz * fga * fx * pb_xxx

                     - 3.0 * pa_xxz * pb_xz * fz * fgb

                     + 2.5 * fx * fx * fz * pb_xxx

                     + 18.0 * pa_xxz * fz * pb_xz * fx

                     + 6.0 * pa_xx * fz * fx * pb_xxx

                     + 36.0 * pa_xz * fx * fz * pb_xxz

                     - fz * fga * pa_z * pb_xxxz

                     + 6.0 * fx * fz * pa_z * pb_xxxz

                     + 14.0 * pa_xxz * fz * pb_xxxz);

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
    inline double fvec_xxz_xxyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
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
        return r_0_0 * (- fx * fx * pa_z * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_z

                     - 0.25 * fz * fga * pa_z * fx * fx

                     - pa_xxz * fx * fz * fgb

                     - 2.0 * pa_xz * fx * pb_x * fz * fgb

                     + 2.5 * pa_xxz * fz * fx * fx

                     + 10.0 * pa_xz * fx * fx * fz * pb_x

                     + 7.5 * fx * fx * fz * pa_z * pb_yy

                     - 0.5 * fx * pa_z * pb_xx * fz * fgb

                     - 0.5 * fx * pa_z * fz * fgb * pb_yy

                     - 0.5 * fz * fga * pa_z * pb_xx * fx

                     - 0.5 * fz * fga * pa_z * fx * pb_yy

                     - pa_xxz * pb_xx * fz * fgb

                     - pa_xxz * fz * fgb * pb_yy

                     + 2.5 * fx * fx * fz * pa_z * pb_xx

                     + 6.0 * pa_xxz * fz * pb_xx * fx

                     + 6.0 * pa_xxz * fz * fx * pb_yy

                     + 24.0 * pa_xz * fx * fz * pb_xyy

                     - fz * fga * pa_z * pb_xxyy

                     + 6.0 * fx * fz * pa_z * pb_xxyy

                     + 14.0 * pa_xxz * fz * pb_xxyy);

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
    inline double fvec_xxz_xxyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
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
                                    double r_0_0)
    {
        return r_0_0 * (3.0 * fx * fx * fx * fz * pb_y

                     - 0.25 * fx * fx * fz * fgb * pb_y

                     - 0.25 * fz * fga * fx * fx * pb_y

                     - 0.5 * pa_xx * fx * fz * fgb * pb_y

                     + 2.5 * pa_xx * fz * fx * fx * pb_y

                     + 10.0 * pa_x * fx * fx * fz * pb_xy

                     + 7.5 * fx * fx * fz * pa_z * pb_yz

                     - 0.5 * fx * pa_z * fz * fgb * pb_yz

                     - 0.5 * fz * fga * pa_z * fx * pb_yz

                     - 0.5 * fz * fga * fx * pb_xxy

                     - pa_xxz * fz * fgb * pb_yz

                     + 2.5 * fx * fx * fz * pb_xxy

                     + 6.0 * pa_xxz * fz * fx * pb_yz

                     + 6.0 * pa_xx * fz * fx * pb_xxy

                     + 24.0 * pa_xz * fx * fz * pb_xyz

                     - fz * fga * pa_z * pb_xxyz

                     + 6.0 * fx * fz * pa_z * pb_xxyz

                     + 14.0 * pa_xxz * fz * pb_xxyz);

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
    inline double fvec_xxz_xxzz_r_0(double fga,
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
                                    double pb_xxzz,
                                    double pb_xz,
                                    double pb_xzz,
                                    double pb_z,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- fx * fx * pa_z * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_z

                     + 6.0 * fx * fx * fx * fz * pb_z

                     - 0.5 * fx * fx * fz * fgb * pb_z

                     - 0.25 * fz * fga * pa_z * fx * fx

                     - 0.5 * fz * fga * fx * fx * pb_z

                     - pa_xxz * fx * fz * fgb

                     - pa_xx * fx * fz * fgb * pb_z

                     - 2.0 * pa_xz * fx * pb_x * fz * fgb

                     + 2.5 * pa_xxz * fz * fx * fx

                     + 5.0 * pa_xx * fz * fx * fx * pb_z

                     + 10.0 * pa_xz * fx * fx * fz * pb_x

                     + 20.0 * pa_x * fx * fx * fz * pb_xz

                     + 7.5 * fx * fx * fz * pa_z * pb_zz

                     - 0.5 * fx * pa_z * pb_xx * fz * fgb

                     - 0.5 * fx * pa_z * fz * fgb * pb_zz

                     - 0.5 * fz * fga * pa_z * pb_xx * fx

                     - 0.5 * fz * fga * pa_z * fx * pb_zz

                     - fz * fga * fx * pb_xxz

                     - pa_xxz * pb_xx * fz * fgb

                     - pa_xxz * fz * fgb * pb_zz

                     + 2.5 * fx * fx * fz * pa_z * pb_xx

                     + 5.0 * fx * fx * fz * pb_xxz

                     + 6.0 * pa_xxz * fz * pb_xx * fx

                     + 6.0 * pa_xxz * fz * fx * pb_zz

                     + 12.0 * pa_xx * fz * fx * pb_xxz

                     + 24.0 * pa_xz * fx * fz * pb_xzz

                     - fz * fga * pa_z * pb_xxzz

                     + 6.0 * fx * fz * pa_z * pb_xxzz

                     + 14.0 * pa_xxz * fz * pb_xxzz);

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
    inline double fvec_xxz_xyyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxz,
                                    double pa_xz,
                                    double pa_z,
                                    double pb_xy,
                                    double pb_xyyy,
                                    double pb_y,
                                    double pb_yyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xz * fx * pb_y * fz * fgb

                     + 15.0 * pa_xz * fx * fx * fz * pb_y

                     - 1.5 * fx * pa_z * pb_xy * fz * fgb

                     - 1.5 * fz * fga * pa_z * pb_xy * fx

                     - 3.0 * pa_xxz * pb_xy * fz * fgb

                     + 7.5 * fx * fx * fz * pa_z * pb_xy

                     + 18.0 * pa_xxz * fz * pb_xy * fx

                     + 12.0 * pa_xz * fx * fz * pb_yyy

                     - fz * fga * pa_z * pb_xyyy

                     + 6.0 * fx * fz * pa_z * pb_xyyy

                     + 14.0 * pa_xxz * fz * pb_xyyy);

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
    inline double fvec_xxz_xyyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
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
                                    double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_x * fx * fx * fz * fgb

                     + 2.0 * pa_x * fx * fx * fx * fz

                     - 0.25 * fx * fx * pb_x * fz * fgb

                     - 0.25 * fz * fga * fx * fx * pb_x

                     - 0.5 * pa_xx * fx * pb_x * fz * fgb

                     - pa_xz * fx * fz * fgb * pb_z

                     + fx * fx * fx * fz * pb_x

                     + 2.5 * pa_xx * fz * fx * fx * pb_x

                     + 5.0 * pa_xz * fx * fx * fz * pb_z

                     + 5.0 * pa_x * fx * fx * fz * pb_yy

                     - 0.5 * fx * pa_z * pb_xz * fz * fgb

                     - 0.5 * fz * fga * pa_z * pb_xz * fx

                     - 0.5 * fz * fga * fx * pb_xyy

                     - pa_xxz * pb_xz * fz * fgb

                     + 2.5 * fx * fx * fz * pa_z * pb_xz

                     + 2.5 * fx * fx * fz * pb_xyy

                     + 6.0 * pa_xxz * fz * pb_xz * fx

                     + 6.0 * pa_xx * fz * fx * pb_xyy

                     + 12.0 * pa_xz * fx * fz * pb_yyz

                     - fz * fga * pa_z * pb_xyyz

                     + 6.0 * fx * fz * pa_z * pb_xyyz

                     + 14.0 * pa_xxz * fz * pb_xyyz);

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
    inline double fvec_xxz_xyzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
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
                                    double r_0_0)
    {
        return r_0_0 * (- pa_xz * fx * pb_y * fz * fgb

                     + 5.0 * pa_xz * fx * fx * fz * pb_y

                     + 10.0 * pa_x * fx * fx * fz * pb_yz

                     - 0.5 * fx * pa_z * pb_xy * fz * fgb

                     - 0.5 * fz * fga * pa_z * pb_xy * fx

                     - fz * fga * fx * pb_xyz

                     - pa_xxz * pb_xy * fz * fgb

                     + 2.5 * fx * fx * fz * pa_z * pb_xy

                     + 5.0 * fx * fx * fz * pb_xyz

                     + 6.0 * pa_xxz * fz * pb_xy * fx

                     + 12.0 * pa_xx * fz * fx * pb_xyz

                     + 12.0 * pa_xz * fx * fz * pb_yzz

                     - fz * fga * pa_z * pb_xyzz

                     + 6.0 * fx * fz * pa_z * pb_xyzz

                     + 14.0 * pa_xxz * fz * pb_xyzz);

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
    inline double fvec_xxz_xzzz_r_0(double fga,
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
                                    double pb_xzzz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fx * fz * fgb

                     + 6.0 * pa_x * fx * fx * fx * fz

                     - 0.75 * fx * fx * pb_x * fz * fgb

                     - 0.75 * fz * fga * fx * fx * pb_x

                     - 1.5 * pa_xx * fx * pb_x * fz * fgb

                     - 3.0 * pa_xz * fx * pb_z * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_x

                     + 7.5 * pa_xx * fz * fx * fx * pb_x

                     + 15.0 * pa_xz * fx * fx * fz * pb_z

                     + 15.0 * pa_x * fx * fx * fz * pb_zz

                     - 1.5 * fx * pa_z * pb_xz * fz * fgb

                     - 1.5 * fz * fga * pa_z * pb_xz * fx

                     - 1.5 * fz * fga * fx * pb_xzz

                     - 3.0 * pa_xxz * pb_xz * fz * fgb

                     + 7.5 * fx * fx * fz * pa_z * pb_xz

                     + 7.5 * fx * fx * fz * pb_xzz

                     + 18.0 * pa_xxz * fz * pb_xz * fx

                     + 18.0 * pa_xx * fz * fx * pb_xzz

                     + 12.0 * pa_xz * fx * fz * pb_zzz

                     - fz * fga * pa_z * pb_xzzz

                     + 6.0 * fx * fz * pa_z * pb_xzzz

                     + 14.0 * pa_xxz * fz * pb_xzzz);

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
    inline double fvec_xxz_yyyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xxz,
                                    double pa_z,
                                    double pb_yy,
                                    double pb_yyyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_z * fz * fgb

                     - 0.75 * fz * fga * pa_z * fx * fx

                     - 3.0 * pa_xxz * fx * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_z

                     + 7.5 * pa_xxz * fz * fx * fx

                     - 3.0 * fx * pa_z * pb_yy * fz * fgb

                     - 3.0 * fz * fga * pa_z * pb_yy * fx

                     - 6.0 * pa_xxz * pb_yy * fz * fgb

                     + 15.0 * fx * fx * fz * pa_z * pb_yy

                     + 36.0 * pa_xxz * fz * pb_yy * fx

                     - fz * fga * pa_z * pb_yyyy

                     + 6.0 * fx * fz * pa_z * pb_yyyy

                     + 14.0 * pa_xxz * fz * pb_yyyy);

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
    inline double fvec_xxz_yyyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxz,
                                    double pa_z,
                                    double pb_y,
                                    double pb_yyy,
                                    double pb_yyyz,
                                    double pb_yz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pb_y * fz * fgb

                     - 0.75 * fz * fga * fx * fx * pb_y

                     - 1.5 * pa_xx * fx * pb_y * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_y

                     + 7.5 * pa_xx * fz * fx * fx * pb_y

                     - 1.5 * fx * pa_z * pb_yz * fz * fgb

                     - 1.5 * fz * fga * pa_z * pb_yz * fx

                     - 0.5 * fz * fga * fx * pb_yyy

                     - 3.0 * pa_xxz * pb_yz * fz * fgb

                     + 7.5 * fx * fx * fz * pa_z * pb_yz

                     + 2.5 * fx * fx * fz * pb_yyy

                     + 18.0 * pa_xxz * fz * pb_yz * fx

                     + 6.0 * pa_xx * fz * fx * pb_yyy

                     - fz * fga * pa_z * pb_yyyz

                     + 6.0 * fx * fz * pa_z * pb_yyyz

                     + 14.0 * pa_xxz * fz * pb_yyyz);

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
    inline double fvec_xxz_yyzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxz,
                                    double pa_z,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_yyzz,
                                    double pb_z,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * pa_z * fz * fgb

                     - 0.5 * fx * fx * fz * fgb * pb_z

                     - 0.25 * fz * fga * pa_z * fx * fx

                     - 0.5 * fz * fga * fx * fx * pb_z

                     - pa_xxz * fx * fz * fgb

                     - pa_xx * fx * fz * fgb * pb_z

                     + fx * fx * fx * fz * pa_z

                     + 2.0 * fx * fx * fx * fz * pb_z

                     + 2.5 * pa_xxz * fz * fx * fx

                     + 5.0 * pa_xx * fz * fx * fx * pb_z

                     - 0.5 * fx * pa_z * pb_yy * fz * fgb

                     - 0.5 * fx * pa_z * fz * fgb * pb_zz

                     - 0.5 * fz * fga * pa_z * pb_yy * fx

                     - 0.5 * fz * fga * pa_z * fx * pb_zz

                     - fz * fga * fx * pb_yyz

                     - pa_xxz * pb_yy * fz * fgb

                     - pa_xxz * fz * fgb * pb_zz

                     + 2.5 * fx * fx * fz * pa_z * pb_yy

                     + 2.5 * fx * fx * fz * pa_z * pb_zz

                     + 5.0 * fx * fx * fz * pb_yyz

                     + 6.0 * pa_xxz * fz * pb_yy * fx

                     + 6.0 * pa_xxz * fz * fx * pb_zz

                     + 12.0 * pa_xx * fz * fx * pb_yyz

                     - fz * fga * pa_z * pb_yyzz

                     + 6.0 * fx * fz * pa_z * pb_yyzz

                     + 14.0 * pa_xxz * fz * pb_yyzz);

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
    inline double fvec_xxz_yzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxz,
                                    double pa_z,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double pb_yzzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pb_y * fz * fgb

                     - 0.75 * fz * fga * fx * fx * pb_y

                     - 1.5 * pa_xx * fx * pb_y * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_y

                     + 7.5 * pa_xx * fz * fx * fx * pb_y

                     - 1.5 * fx * pa_z * pb_yz * fz * fgb

                     - 1.5 * fz * fga * pa_z * pb_yz * fx

                     - 1.5 * fz * fga * fx * pb_yzz

                     - 3.0 * pa_xxz * pb_yz * fz * fgb

                     + 7.5 * fx * fx * fz * pa_z * pb_yz

                     + 7.5 * fx * fx * fz * pb_yzz

                     + 18.0 * pa_xxz * fz * pb_yz * fx

                     + 18.0 * pa_xx * fz * fx * pb_yzz

                     - fz * fga * pa_z * pb_yzzz

                     + 6.0 * fx * fz * pa_z * pb_yzzz

                     + 14.0 * pa_xxz * fz * pb_yzzz);

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
    inline double fvec_xxz_zzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xx,
                                    double pa_xxz,
                                    double pa_z,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double pb_zzzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_z * fz * fgb

                     - 3.0 * fx * fx * pb_z * fz * fgb

                     - 0.75 * fz * fga * pa_z * fx * fx

                     - 3.0 * fz * fga * fx * fx * pb_z

                     - 3.0 * pa_xxz * fx * fz * fgb

                     - 6.0 * pa_xx * fx * pb_z * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_z

                     + 12.0 * fx * fx * fx * fz * pb_z

                     + 7.5 * pa_xxz * fz * fx * fx

                     + 30.0 * pa_xx * fz * fx * fx * pb_z

                     - 3.0 * fx * pa_z * pb_zz * fz * fgb

                     - 3.0 * fz * fga * pa_z * pb_zz * fx

                     - 2.0 * fz * fga * fx * pb_zzz

                     - 6.0 * pa_xxz * pb_zz * fz * fgb

                     + 15.0 * fx * fx * fz * pa_z * pb_zz

                     + 10.0 * fx * fx * fz * pb_zzz

                     + 36.0 * pa_xxz * fz * pb_zz * fx

                     + 24.0 * pa_xx * fz * fx * pb_zzz

                     - fz * fga * pa_z * pb_zzzz

                     + 6.0 * fx * fz * pa_z * pb_zzzz

                     + 14.0 * pa_xxz * fz * pb_zzzz);

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
    inline double fvec_xyy_xxxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xyy,
                                    double pa_yy,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double pb_xxxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fx * fz * fgb

                     - 0.75 * pa_x * fz * fga * fx * fx

                     - 3.0 * fx * fx * pb_x * fz * fgb

                     - 3.0 * fx * fx * fz * fga * pb_x

                     - 3.0 * pa_xyy * fx * fz * fgb

                     + 3.0 * pa_x * fz * fx * fx * fx

                     - 6.0 * fx * pa_yy * pb_x * fz * fgb

                     + 12.0 * fx * fx * fx * fz * pb_x

                     + 7.5 * pa_xyy * fz * fx * fx

                     + 30.0 * fx * fx * pa_yy * fz * pb_x

                     - 3.0 * pa_x * fx * pb_xx * fz * fgb

                     - 3.0 * pa_x * fz * fga * pb_xx * fx

                     - 2.0 * fx * fz * fga * pb_xxx

                     - 6.0 * pa_xyy * pb_xx * fz * fgb

                     + 15.0 * pa_x * fz * fx * fx * pb_xx

                     + 10.0 * fx * fx * fz * pb_xxx

                     + 36.0 * pa_xyy * fz * pb_xx * fx

                     + 24.0 * fx * pa_yy * fz * pb_xxx

                     - pa_x * fz * fga * pb_xxxx

                     + 6.0 * pa_x * fz * fx * pb_xxxx

                     + 14.0 * pa_xyy * fz * pb_xxxx);

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
    inline double fvec_xyy_xxxy_r_0(double fga,
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
                                    double pb_xxx,
                                    double pb_xxxy,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_y * fz * fgb

                     + 6.0 * fx * fx * fx * pa_y * fz

                     - 0.75 * fx * fx * fz * fgb * pb_y

                     - 0.75 * fx * fx * fz * fga * pb_y

                     - 3.0 * pa_xy * fx * pb_x * fz * fgb

                     - 1.5 * fx * pa_yy * fz * fgb * pb_y

                     + 3.0 * fx * fx * fx * fz * pb_y

                     + 15.0 * pa_xy * fz * fx * fx * pb_x

                     + 7.5 * fx * fx * pa_yy * fz * pb_y

                     + 15.0 * fx * fx * pa_y * fz * pb_xx

                     - 1.5 * pa_x * fx * pb_xy * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_xy * fx

                     - 1.5 * fx * fz * fga * pb_xxy

                     - 3.0 * pa_xyy * pb_xy * fz * fgb

                     + 7.5 * pa_x * fz * fx * fx * pb_xy

                     + 7.5 * fx * fx * fz * pb_xxy

                     + 18.0 * pa_xyy * fz * pb_xy * fx

                     + 12.0 * pa_xy * fz * fx * pb_xxx

                     + 18.0 * fx * pa_yy * fz * pb_xxy

                     - pa_x * fz * fga * pb_xxxy

                     + 6.0 * pa_x * fz * fx * pb_xxxy

                     + 14.0 * pa_xyy * fz * pb_xxxy);

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
    inline double fvec_xyy_xxxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xyy,
                                    double pa_yy,
                                    double pb_xxxz,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb * pb_z

                     - 0.75 * fx * fx * fz * fga * pb_z

                     - 1.5 * fx * pa_yy * fz * fgb * pb_z

                     + 3.0 * fx * fx * fx * fz * pb_z

                     + 7.5 * fx * fx * pa_yy * fz * pb_z

                     - 1.5 * pa_x * fx * pb_xz * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_xz * fx

                     - 1.5 * fx * fz * fga * pb_xxz

                     - 3.0 * pa_xyy * pb_xz * fz * fgb

                     + 7.5 * pa_x * fz * fx * fx * pb_xz

                     + 7.5 * fx * fx * fz * pb_xxz

                     + 18.0 * pa_xyy * fz * pb_xz * fx

                     + 18.0 * fx * pa_yy * fz * pb_xxz

                     - pa_x * fz * fga * pb_xxxz

                     + 6.0 * pa_x * fz * fx * pb_xxxz

                     + 14.0 * pa_xyy * fz * pb_xxxz);

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
    inline double fvec_xyy_xxyy_r_0(double fga,
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
                                    double pb_xxyy,
                                    double pb_xy,
                                    double pb_xyy,
                                    double pb_y,
                                    double pb_yy,
                                    double r_0_0)
    {
        return r_0_0 * (- pa_x * fx * fx * fz * fgb

                     + 3.0 * pa_x * fx * fx * fx * fz

                     + 6.0 * fx * fx * fx * fz * pb_x

                     - 0.25 * pa_x * fz * fga * fx * fx

                     - 0.5 * fx * fx * pb_x * fz * fgb

                     - 0.5 * fx * fx * fz * fga * pb_x

                     - pa_xyy * fx * fz * fgb

                     - 2.0 * pa_xy * fx * fz * fgb * pb_y

                     - fx * pa_yy * pb_x * fz * fgb

                     + 2.5 * pa_xyy * fz * fx * fx

                     + 10.0 * pa_xy * fz * fx * fx * pb_y

                     + 7.5 * pa_x * fx * fx * fz * pb_xx

                     + 5.0 * fx * fx * pa_yy * fz * pb_x

                     + 20.0 * fx * fx * pa_y * fz * pb_xy

                     - 0.5 * pa_x * fx * pb_xx * fz * fgb

                     - 0.5 * pa_x * fx * fz * fgb * pb_yy

                     - 0.5 * pa_x * fz * fga * pb_xx * fx

                     - 0.5 * pa_x * fz * fga * fx * pb_yy

                     - fx * fz * fga * pb_xyy

                     - pa_xyy * pb_xx * fz * fgb

                     - pa_xyy * fz * fgb * pb_yy

                     + 2.5 * pa_x * fz * fx * fx * pb_yy

                     + 5.0 * fx * fx * fz * pb_xyy

                     + 6.0 * pa_xyy * fz * pb_xx * fx

                     + 6.0 * pa_xyy * fz * fx * pb_yy

                     + 24.0 * pa_xy * fz * fx * pb_xxy

                     + 12.0 * fx * pa_yy * fz * pb_xyy

                     - pa_x * fz * fga * pb_xxyy

                     + 6.0 * pa_x * fz * fx * pb_xxyy

                     + 14.0 * pa_xyy * fz * pb_xxyy);

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
    inline double fvec_xyy_xxyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
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
                                    double r_0_0)
    {
        return r_0_0 * (- pa_xy * fx * fz * fgb * pb_z

                     + 5.0 * pa_xy * fz * fx * fx * pb_z

                     + 10.0 * fx * fx * pa_y * fz * pb_xz

                     - 0.5 * pa_x * fx * fz * fgb * pb_yz

                     - 0.5 * pa_x * fz * fga * fx * pb_yz

                     - fx * fz * fga * pb_xyz

                     - pa_xyy * fz * fgb * pb_yz

                     + 2.5 * pa_x * fz * fx * fx * pb_yz

                     + 5.0 * fx * fx * fz * pb_xyz

                     + 6.0 * pa_xyy * fz * fx * pb_yz

                     + 12.0 * pa_xy * fz * fx * pb_xxz

                     + 12.0 * fx * pa_yy * fz * pb_xyz

                     - pa_x * fz * fga * pb_xxyz

                     + 6.0 * pa_x * fz * fx * pb_xxyz

                     + 14.0 * pa_xyy * fz * pb_xxyz);

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
    inline double fvec_xyy_xxzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xyy,
                                    double pa_yy,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxzz,
                                    double pb_xzz,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_x * fx * fx * fz * fgb

                     - 0.25 * pa_x * fz * fga * fx * fx

                     - 0.5 * fx * fx * pb_x * fz * fgb

                     - 0.5 * fx * fx * fz * fga * pb_x

                     - pa_xyy * fx * fz * fgb

                     + pa_x * fz * fx * fx * fx

                     - fx * pa_yy * pb_x * fz * fgb

                     + 2.0 * fx * fx * fx * fz * pb_x

                     + 2.5 * pa_xyy * fz * fx * fx

                     + 5.0 * fx * fx * pa_yy * fz * pb_x

                     - 0.5 * pa_x * fx * pb_xx * fz * fgb

                     - 0.5 * pa_x * fx * fz * fgb * pb_zz

                     - 0.5 * pa_x * fz * fga * pb_xx * fx

                     - 0.5 * pa_x * fz * fga * fx * pb_zz

                     - fx * fz * fga * pb_xzz

                     - pa_xyy * pb_xx * fz * fgb

                     - pa_xyy * fz * fgb * pb_zz

                     + 2.5 * pa_x * fz * fx * fx * pb_xx

                     + 2.5 * pa_x * fz * fx * fx * pb_zz

                     + 5.0 * fx * fx * fz * pb_xzz

                     + 6.0 * pa_xyy * fz * pb_xx * fx

                     + 6.0 * pa_xyy * fz * fx * pb_zz

                     + 12.0 * fx * pa_yy * fz * pb_xzz

                     - pa_x * fz * fga * pb_xxzz

                     + 6.0 * pa_x * fz * fx * pb_xxzz

                     + 14.0 * pa_xyy * fz * pb_xxzz);

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
    inline double fvec_xyy_xyyy_r_0(double fga,
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
                                    double pb_xyyy,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_y * fz * fgb

                     + 6.0 * fx * fx * fx * pa_y * fz

                     + 9.0 * fx * fx * fx * fz * pb_y

                     - 0.75 * fx * fx * pb_y * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pb_y

                     - 3.0 * pa_xy * fx * pb_x * fz * fgb

                     - 1.5 * fx * pa_yy * pb_y * fz * fgb

                     + 15.0 * pa_xy * fz * fx * fx * pb_x

                     + 22.5 * pa_x * fx * fx * fz * pb_xy

                     + 7.5 * fx * fx * pa_yy * fz * pb_y

                     + 15.0 * fx * fx * pa_y * fz * pb_yy

                     - 1.5 * pa_x * fx * pb_xy * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_xy * fx

                     - 0.5 * fx * fz * fga * pb_yyy

                     - 3.0 * pa_xyy * pb_xy * fz * fgb

                     + 2.5 * fx * fx * fz * pb_yyy

                     + 18.0 * pa_xyy * fz * pb_xy * fx

                     + 36.0 * pa_xy * fz * fx * pb_xyy

                     + 6.0 * fx * pa_yy * fz * pb_yyy

                     - pa_x * fz * fga * pb_xyyy

                     + 6.0 * pa_x * fz * fx * pb_xyyy

                     + 14.0 * pa_xyy * fz * pb_xyyy);

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
    inline double fvec_xyy_xyyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
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
                                    double r_0_0)
    {
        return r_0_0 * (3.0 * fx * fx * fx * fz * pb_z

                     - 0.25 * fx * fx * fz * fgb * pb_z

                     - 0.25 * fx * fx * fz * fga * pb_z

                     - 0.5 * fx * pa_yy * fz * fgb * pb_z

                     + 7.5 * pa_x * fx * fx * fz * pb_xz

                     + 2.5 * fx * fx * pa_yy * fz * pb_z

                     + 10.0 * fx * fx * pa_y * fz * pb_yz

                     - 0.5 * pa_x * fx * pb_xz * fz * fgb

                     - 0.5 * pa_x * fz * fga * pb_xz * fx

                     - 0.5 * fx * fz * fga * pb_yyz

                     - pa_xyy * pb_xz * fz * fgb

                     + 2.5 * fx * fx * fz * pb_yyz

                     + 6.0 * pa_xyy * fz * pb_xz * fx

                     + 24.0 * pa_xy * fz * fx * pb_xyz

                     + 6.0 * fx * pa_yy * fz * pb_yyz

                     - pa_x * fz * fga * pb_xyyz

                     + 6.0 * pa_x * fz * fx * pb_xyyz

                     + 14.0 * pa_xyy * fz * pb_xyyz);

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
    inline double fvec_xyy_xyzz_r_0(double fga,
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
                                    double pb_xyzz,
                                    double pb_xzz,
                                    double pb_y,
                                    double pb_yzz,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * pa_y * fz * fgb

                     + 2.0 * fx * fx * fx * pa_y * fz

                     - 0.25 * fx * fx * pb_y * fz * fgb

                     - 0.25 * fx * fx * fz * fga * pb_y

                     - pa_xy * fx * pb_x * fz * fgb

                     - 0.5 * fx * pa_yy * pb_y * fz * fgb

                     + fx * fx * fx * fz * pb_y

                     + 5.0 * pa_xy * fz * fx * fx * pb_x

                     + 2.5 * fx * fx * pa_yy * fz * pb_y

                     + 5.0 * fx * fx * pa_y * fz * pb_zz

                     - 0.5 * pa_x * fx * pb_xy * fz * fgb

                     - 0.5 * pa_x * fz * fga * pb_xy * fx

                     - 0.5 * fx * fz * fga * pb_yzz

                     - pa_xyy * pb_xy * fz * fgb

                     + 2.5 * pa_x * fz * fx * fx * pb_xy

                     + 2.5 * fx * fx * fz * pb_yzz

                     + 6.0 * pa_xyy * fz * pb_xy * fx

                     + 12.0 * pa_xy * fz * fx * pb_xzz

                     + 6.0 * fx * pa_yy * fz * pb_yzz

                     - pa_x * fz * fga * pb_xyzz

                     + 6.0 * pa_x * fz * fx * pb_xyzz

                     + 14.0 * pa_xyy * fz * pb_xyzz);

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
    inline double fvec_xyy_xzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xyy,
                                    double pa_yy,
                                    double pb_xz,
                                    double pb_xzzz,
                                    double pb_z,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pb_z * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pb_z

                     - 1.5 * fx * pa_yy * pb_z * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_z

                     + 7.5 * fx * fx * pa_yy * fz * pb_z

                     - 1.5 * pa_x * fx * pb_xz * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_xz * fx

                     - 0.5 * fx * fz * fga * pb_zzz

                     - 3.0 * pa_xyy * pb_xz * fz * fgb

                     + 7.5 * pa_x * fz * fx * fx * pb_xz

                     + 2.5 * fx * fx * fz * pb_zzz

                     + 18.0 * pa_xyy * fz * pb_xz * fx

                     + 6.0 * fx * pa_yy * fz * pb_zzz

                     - pa_x * fz * fga * pb_xzzz

                     + 6.0 * pa_x * fz * fx * pb_xzzz

                     + 14.0 * pa_xyy * fz * pb_xzzz);

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
    inline double fvec_xyy_yyyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double pb_yyyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_x * fx * fx * fz * fgb

                     + 15.0 * pa_x * fx * fx * fx * fz

                     - 0.75 * pa_x * fz * fga * fx * fx

                     - 3.0 * pa_xyy * fx * fz * fgb

                     - 12.0 * pa_xy * fx * pb_y * fz * fgb

                     + 7.5 * pa_xyy * fz * fx * fx

                     + 60.0 * pa_xy * fz * fx * fx * pb_y

                     + 45.0 * pa_x * fx * fx * fz * pb_yy

                     - 3.0 * pa_x * fx * pb_yy * fz * fgb

                     - 3.0 * pa_x * fz * fga * pb_yy * fx

                     - 6.0 * pa_xyy * pb_yy * fz * fgb

                     + 36.0 * pa_xyy * fz * pb_yy * fx

                     + 48.0 * pa_xy * fz * fx * pb_yyy

                     - pa_x * fz * fga * pb_yyyy

                     + 6.0 * pa_x * fz * fx * pb_yyyy

                     + 14.0 * pa_xyy * fz * pb_yyyy);

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
    inline double fvec_xyy_yyyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pb_yyyz,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xy * fx * fz * fgb * pb_z

                     + 15.0 * pa_xy * fz * fx * fx * pb_z

                     + 22.5 * pa_x * fx * fx * fz * pb_yz

                     - 1.5 * pa_x * fx * pb_yz * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_yz * fx

                     - 3.0 * pa_xyy * pb_yz * fz * fgb

                     + 18.0 * pa_xyy * fz * pb_yz * fx

                     + 36.0 * pa_xy * fz * fx * pb_yyz

                     - pa_x * fz * fga * pb_yyyz

                     + 6.0 * pa_x * fz * fx * pb_yyyz

                     + 14.0 * pa_xyy * fz * pb_yyyz);

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
    inline double fvec_xyy_yyzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyzz,
                                    double pb_yzz,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- pa_x * fx * fx * fz * fgb

                     + 3.0 * pa_x * fx * fx * fx * fz

                     - 0.25 * pa_x * fz * fga * fx * fx

                     - pa_xyy * fx * fz * fgb

                     - 2.0 * pa_xy * fx * pb_y * fz * fgb

                     + 2.5 * pa_xyy * fz * fx * fx

                     + 10.0 * pa_xy * fz * fx * fx * pb_y

                     + 7.5 * pa_x * fx * fx * fz * pb_zz

                     - 0.5 * pa_x * fx * pb_yy * fz * fgb

                     - 0.5 * pa_x * fx * fz * fgb * pb_zz

                     - 0.5 * pa_x * fz * fga * pb_yy * fx

                     - 0.5 * pa_x * fz * fga * fx * pb_zz

                     - pa_xyy * pb_yy * fz * fgb

                     - pa_xyy * fz * fgb * pb_zz

                     + 2.5 * pa_x * fz * fx * fx * pb_yy

                     + 6.0 * pa_xyy * fz * pb_yy * fx

                     + 6.0 * pa_xyy * fz * fx * pb_zz

                     + 24.0 * pa_xy * fz * fx * pb_yzz

                     - pa_x * fz * fga * pb_yyzz

                     + 6.0 * pa_x * fz * fx * pb_yyzz

                     + 14.0 * pa_xyy * fz * pb_yyzz);

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
    inline double fvec_xyy_yzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xy,
                                    double pa_xyy,
                                    double pb_yz,
                                    double pb_yzzz,
                                    double pb_z,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xy * fx * pb_z * fz * fgb

                     + 15.0 * pa_xy * fz * fx * fx * pb_z

                     - 1.5 * pa_x * fx * pb_yz * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_yz * fx

                     - 3.0 * pa_xyy * pb_yz * fz * fgb

                     + 7.5 * pa_x * fz * fx * fx * pb_yz

                     + 18.0 * pa_xyy * fz * pb_yz * fx

                     + 12.0 * pa_xy * fz * fx * pb_zzz

                     - pa_x * fz * fga * pb_yzzz

                     + 6.0 * pa_x * fz * fx * pb_yzzz

                     + 14.0 * pa_xyy * fz * pb_yzzz);

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
    inline double fvec_xyy_zzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xyy,
                                    double pb_zz,
                                    double pb_zzzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fx * fz * fgb

                     - 0.75 * pa_x * fz * fga * fx * fx

                     - 3.0 * pa_xyy * fx * fz * fgb

                     + 3.0 * pa_x * fz * fx * fx * fx

                     + 7.5 * pa_xyy * fz * fx * fx

                     - 3.0 * pa_x * fx * pb_zz * fz * fgb

                     - 3.0 * pa_x * fz * fga * pb_zz * fx

                     - 6.0 * pa_xyy * pb_zz * fz * fgb

                     + 15.0 * pa_x * fz * fx * fx * pb_zz

                     + 36.0 * pa_xyy * fz * pb_zz * fx

                     - pa_x * fz * fga * pb_zzzz

                     + 6.0 * pa_x * fz * fx * pb_zzzz

                     + 14.0 * pa_xyy * fz * pb_zzzz);

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
    inline double fvec_xyz_xxxx_r_0(double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xyz,
                                    double pa_yz,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double pb_xxxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xyz * fx * fz * fgb

                     - 6.0 * fx * pa_yz * pb_x * fz * fgb

                     + 7.5 * pa_xyz * fz * fx * fx

                     + 30.0 * fx * fx * pa_yz * fz * pb_x

                     - 6.0 * pa_xyz * pb_xx * fz * fgb

                     + 36.0 * pa_xyz * fz * pb_xx * fx

                     + 24.0 * fx * pa_yz * fz * pb_xxx

                     + 14.0 * pa_xyz * fz * pb_xxxx);

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
    inline double fvec_xyz_xxxy_r_0(double fgb,
                                    double fx,
                                    double fz,
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
        return r_0_0 * (- 0.75 * fx * fx * pa_z * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_z

                     - 1.5 * pa_xz * fx * pb_x * fz * fgb

                     - 1.5 * fx * pa_yz * fz * fgb * pb_y

                     + 7.5 * pa_xz * fx * fx * fz * pb_x

                     + 7.5 * fx * fx * pa_yz * fz * pb_y

                     + 7.5 * fx * fx * fz * pa_z * pb_xx

                     - 3.0 * pa_xyz * pb_xy * fz * fgb

                     + 18.0 * pa_xyz * fz * pb_xy * fx

                     + 6.0 * pa_xz * fx * fz * pb_xxx

                     + 18.0 * fx * pa_yz * fz * pb_xxy

                     + 14.0 * pa_xyz * fz * pb_xxxy);

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
    inline double fvec_xyz_xxxz_r_0(double fgb,
                                    double fx,
                                    double fz,
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
        return r_0_0 * (- 0.75 * fx * fx * pa_y * fz * fgb

                     + 3.0 * fx * fx * fx * pa_y * fz

                     - 1.5 * pa_xy * fx * pb_x * fz * fgb

                     - 1.5 * fx * pa_yz * fz * fgb * pb_z

                     + 7.5 * pa_xy * fz * fx * fx * pb_x

                     + 7.5 * fx * fx * pa_yz * fz * pb_z

                     + 7.5 * fx * fx * pa_y * fz * pb_xx

                     - 3.0 * pa_xyz * pb_xz * fz * fgb

                     + 18.0 * pa_xyz * fz * pb_xz * fx

                     + 6.0 * pa_xy * fz * fx * pb_xxx

                     + 18.0 * fx * pa_yz * fz * pb_xxz

                     + 14.0 * pa_xyz * fz * pb_xxxz);

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
    inline double fvec_xyz_xxyy_r_0(double fgb,
                                    double fx,
                                    double fz,
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
        return r_0_0 * (- pa_xyz * fx * fz * fgb

                     - pa_xz * fx * fz * fgb * pb_y

                     - fx * pa_yz * pb_x * fz * fgb

                     + 2.5 * pa_xyz * fz * fx * fx

                     + 5.0 * pa_xz * fx * fx * fz * pb_y

                     + 5.0 * fx * fx * pa_yz * fz * pb_x

                     + 10.0 * fx * fx * fz * pa_z * pb_xy

                     - pa_xyz * pb_xx * fz * fgb

                     - pa_xyz * fz * fgb * pb_yy

                     + 6.0 * pa_xyz * fz * pb_xx * fx

                     + 6.0 * pa_xyz * fz * fx * pb_yy

                     + 12.0 * pa_xz * fx * fz * pb_xxy

                     + 12.0 * fx * pa_yz * fz * pb_xyy

                     + 14.0 * pa_xyz * fz * pb_xxyy);

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
    inline double fvec_xyz_xxyz_r_0(double fgb,
                                    double fx,
                                    double fz,
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
                                    double r_0_0)
    {
        return r_0_0 * (- 0.25 * pa_x * fx * fx * fz * fgb

                     + pa_x * fx * fx * fx * fz

                     + 2.0 * fx * fx * fx * fz * pb_x

                     - 0.5 * pa_xy * fx * fz * fgb * pb_y

                     - 0.5 * pa_xz * fx * fz * fgb * pb_z

                     + 2.5 * pa_xy * fz * fx * fx * pb_y

                     + 2.5 * pa_xz * fx * fx * fz * pb_z

                     + 2.5 * pa_x * fx * fx * fz * pb_xx

                     + 5.0 * fx * fx * pa_y * fz * pb_xy

                     + 5.0 * fx * fx * fz * pa_z * pb_xz

                     - pa_xyz * fz * fgb * pb_yz

                     + 6.0 * pa_xyz * fz * fx * pb_yz

                     + 6.0 * pa_xy * fz * fx * pb_xxy

                     + 6.0 * pa_xz * fx * fz * pb_xxz

                     + 12.0 * fx * pa_yz * fz * pb_xyz

                     + 14.0 * pa_xyz * fz * pb_xxyz);

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
    inline double fvec_xyz_xxzz_r_0(double fgb,
                                    double fx,
                                    double fz,
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
        return r_0_0 * (- pa_xyz * fx * fz * fgb

                     - pa_xy * fx * fz * fgb * pb_z

                     - fx * pa_yz * pb_x * fz * fgb

                     + 2.5 * pa_xyz * fz * fx * fx

                     + 5.0 * pa_xy * fz * fx * fx * pb_z

                     + 5.0 * fx * fx * pa_yz * fz * pb_x

                     + 10.0 * fx * fx * pa_y * fz * pb_xz

                     - pa_xyz * pb_xx * fz * fgb

                     - pa_xyz * fz * fgb * pb_zz

                     + 6.0 * pa_xyz * fz * pb_xx * fx

                     + 6.0 * pa_xyz * fz * fx * pb_zz

                     + 12.0 * pa_xy * fz * fx * pb_xxz

                     + 12.0 * fx * pa_yz * fz * pb_xzz

                     + 14.0 * pa_xyz * fz * pb_xxzz);

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
    inline double fvec_xyz_xyyy_r_0(double fgb,
                                    double fx,
                                    double fz,
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
        return r_0_0 * (- 0.75 * fx * fx * pa_z * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_z

                     - 1.5 * pa_xz * fx * pb_x * fz * fgb

                     - 1.5 * fx * pa_yz * pb_y * fz * fgb

                     + 7.5 * pa_xz * fx * fx * fz * pb_x

                     + 7.5 * fx * fx * pa_yz * fz * pb_y

                     + 7.5 * fx * fx * fz * pa_z * pb_yy

                     - 3.0 * pa_xyz * pb_xy * fz * fgb

                     + 18.0 * pa_xyz * fz * pb_xy * fx

                     + 18.0 * pa_xz * fx * fz * pb_xyy

                     + 6.0 * fx * pa_yz * fz * pb_yyy

                     + 14.0 * pa_xyz * fz * pb_xyyy);

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
    inline double fvec_xyz_xyyz_r_0(double fgb,
                                    double fx,
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
        return r_0_0 * (- 0.25 * fx * fx * pa_y * fz * fgb

                     + fx * fx * fx * pa_y * fz

                     + 2.0 * fx * fx * fx * fz * pb_y

                     - 0.5 * pa_xy * fx * pb_x * fz * fgb

                     - 0.5 * fx * pa_yz * fz * fgb * pb_z

                     + 2.5 * pa_xy * fz * fx * fx * pb_x

                     + 5.0 * pa_x * fx * fx * fz * pb_xy

                     + 2.5 * fx * fx * pa_yz * fz * pb_z

                     + 2.5 * fx * fx * pa_y * fz * pb_yy

                     + 5.0 * fx * fx * fz * pa_z * pb_yz

                     - pa_xyz * pb_xz * fz * fgb

                     + 6.0 * pa_xyz * fz * pb_xz * fx

                     + 6.0 * pa_xy * fz * fx * pb_xyy

                     + 12.0 * pa_xz * fx * fz * pb_xyz

                     + 6.0 * fx * pa_yz * fz * pb_yyz

                     + 14.0 * pa_xyz * fz * pb_xyyz);

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
    inline double fvec_xyz_xyzz_r_0(double fgb,
                                    double fx,
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
        return r_0_0 * (- 0.25 * fx * fx * pa_z * fz * fgb

                     + fx * fx * fx * fz * pa_z

                     + 2.0 * fx * fx * fx * fz * pb_z

                     - 0.5 * pa_xz * fx * pb_x * fz * fgb

                     - 0.5 * fx * pa_yz * pb_y * fz * fgb

                     + 2.5 * pa_xz * fx * fx * fz * pb_x

                     + 5.0 * pa_x * fx * fx * fz * pb_xz

                     + 2.5 * fx * fx * pa_yz * fz * pb_y

                     + 5.0 * fx * fx * pa_y * fz * pb_yz

                     + 2.5 * fx * fx * fz * pa_z * pb_zz

                     - pa_xyz * pb_xy * fz * fgb

                     + 6.0 * pa_xyz * fz * pb_xy * fx

                     + 12.0 * pa_xy * fz * fx * pb_xyz

                     + 6.0 * pa_xz * fx * fz * pb_xzz

                     + 6.0 * fx * pa_yz * fz * pb_yzz

                     + 14.0 * pa_xyz * fz * pb_xyzz);

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
    inline double fvec_xyz_xzzz_r_0(double fgb,
                                    double fx,
                                    double fz,
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
        return r_0_0 * (- 0.75 * fx * fx * pa_y * fz * fgb

                     + 3.0 * fx * fx * fx * pa_y * fz

                     - 1.5 * pa_xy * fx * pb_x * fz * fgb

                     - 1.5 * fx * pa_yz * pb_z * fz * fgb

                     + 7.5 * pa_xy * fz * fx * fx * pb_x

                     + 7.5 * fx * fx * pa_yz * fz * pb_z

                     + 7.5 * fx * fx * pa_y * fz * pb_zz

                     - 3.0 * pa_xyz * pb_xz * fz * fgb

                     + 18.0 * pa_xyz * fz * pb_xz * fx

                     + 18.0 * pa_xy * fz * fx * pb_xzz

                     + 6.0 * fx * pa_yz * fz * pb_zzz

                     + 14.0 * pa_xyz * fz * pb_xzzz);

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
    inline double fvec_xyz_yyyy_r_0(double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xyz,
                                    double pa_xz,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double pb_yyyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xyz * fx * fz * fgb

                     - 6.0 * pa_xz * fx * pb_y * fz * fgb

                     + 7.5 * pa_xyz * fz * fx * fx

                     + 30.0 * pa_xz * fx * fx * fz * pb_y

                     - 6.0 * pa_xyz * pb_yy * fz * fgb

                     + 36.0 * pa_xyz * fz * pb_yy * fx

                     + 24.0 * pa_xz * fx * fz * pb_yyy

                     + 14.0 * pa_xyz * fz * pb_yyyy);

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
    inline double fvec_xyz_yyyz_r_0(double fgb,
                                    double fx,
                                    double fz,
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
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * pa_x * fx * fx * fz * fgb

                     + 3.0 * pa_x * fx * fx * fx * fz

                     - 1.5 * pa_xy * fx * pb_y * fz * fgb

                     - 1.5 * pa_xz * fx * fz * fgb * pb_z

                     + 7.5 * pa_xy * fz * fx * fx * pb_y

                     + 7.5 * pa_xz * fx * fx * fz * pb_z

                     + 7.5 * pa_x * fx * fx * fz * pb_yy

                     - 3.0 * pa_xyz * pb_yz * fz * fgb

                     + 18.0 * pa_xyz * fz * pb_yz * fx

                     + 6.0 * pa_xy * fz * fx * pb_yyy

                     + 18.0 * pa_xz * fx * fz * pb_yyz

                     + 14.0 * pa_xyz * fz * pb_yyyz);

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
    inline double fvec_xyz_yyzz_r_0(double fgb,
                                    double fx,
                                    double fz,
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
                                    double r_0_0)
    {
        return r_0_0 * (- pa_xyz * fx * fz * fgb

                     - pa_xy * fx * fz * fgb * pb_z

                     - pa_xz * fx * pb_y * fz * fgb

                     + 2.5 * pa_xyz * fz * fx * fx

                     + 5.0 * pa_xy * fz * fx * fx * pb_z

                     + 5.0 * pa_xz * fx * fx * fz * pb_y

                     + 10.0 * pa_x * fx * fx * fz * pb_yz

                     - pa_xyz * pb_yy * fz * fgb

                     - pa_xyz * fz * fgb * pb_zz

                     + 6.0 * pa_xyz * fz * pb_yy * fx

                     + 6.0 * pa_xyz * fz * fx * pb_zz

                     + 12.0 * pa_xy * fz * fx * pb_yyz

                     + 12.0 * pa_xz * fx * fz * pb_yzz

                     + 14.0 * pa_xyz * fz * pb_yyzz);

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
    inline double fvec_xyz_yzzz_r_0(double fgb,
                                    double fx,
                                    double fz,
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
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * pa_x * fx * fx * fz * fgb

                     + 3.0 * pa_x * fx * fx * fx * fz

                     - 1.5 * pa_xy * fx * pb_y * fz * fgb

                     - 1.5 * pa_xz * fx * pb_z * fz * fgb

                     + 7.5 * pa_xy * fz * fx * fx * pb_y

                     + 7.5 * pa_xz * fx * fx * fz * pb_z

                     + 7.5 * pa_x * fx * fx * fz * pb_zz

                     - 3.0 * pa_xyz * pb_yz * fz * fgb

                     + 18.0 * pa_xyz * fz * pb_yz * fx

                     + 18.0 * pa_xy * fz * fx * pb_yzz

                     + 6.0 * pa_xz * fx * fz * pb_zzz

                     + 14.0 * pa_xyz * fz * pb_yzzz);

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
    inline double fvec_xyz_zzzz_r_0(double fgb,
                                    double fx,
                                    double fz,
                                    double pa_xy,
                                    double pa_xyz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double pb_zzzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xyz * fx * fz * fgb

                     - 6.0 * pa_xy * fx * pb_z * fz * fgb

                     + 7.5 * pa_xyz * fz * fx * fx

                     + 30.0 * pa_xy * fz * fx * fx * pb_z

                     - 6.0 * pa_xyz * pb_zz * fz * fgb

                     + 36.0 * pa_xyz * fz * pb_zz * fx

                     + 24.0 * pa_xy * fz * fx * pb_zzz

                     + 14.0 * pa_xyz * fz * pb_zzzz);

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
    inline double fvec_xzz_xxxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xzz,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxx,
                                    double pb_xxxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fx * fz * fgb

                     - 0.75 * pa_x * fz * fga * fx * fx

                     - 3.0 * fx * fx * pb_x * fz * fgb

                     - 3.0 * fx * fx * fz * fga * pb_x

                     - 3.0 * pa_xzz * fx * fz * fgb

                     + 3.0 * pa_x * fz * fx * fx * fx

                     - 6.0 * fx * pa_zz * pb_x * fz * fgb

                     + 12.0 * fx * fx * fx * fz * pb_x

                     + 7.5 * pa_xzz * fz * fx * fx

                     + 30.0 * fx * fx * pa_zz * fz * pb_x

                     - 3.0 * pa_x * fx * pb_xx * fz * fgb

                     - 3.0 * pa_x * fz * fga * pb_xx * fx

                     - 2.0 * fx * fz * fga * pb_xxx

                     - 6.0 * pa_xzz * pb_xx * fz * fgb

                     + 15.0 * pa_x * fz * fx * fx * pb_xx

                     + 10.0 * fx * fx * fz * pb_xxx

                     + 36.0 * pa_xzz * fz * pb_xx * fx

                     + 24.0 * fx * pa_zz * fz * pb_xxx

                     - pa_x * fz * fga * pb_xxxx

                     + 6.0 * pa_x * fz * fx * pb_xxxx

                     + 14.0 * pa_xzz * fz * pb_xxxx);

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
    inline double fvec_xzz_xxxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xzz,
                                    double pa_zz,
                                    double pb_xxxy,
                                    double pb_xxy,
                                    double pb_xy,
                                    double pb_y,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb * pb_y

                     - 0.75 * fx * fx * fz * fga * pb_y

                     - 1.5 * fx * pa_zz * fz * fgb * pb_y

                     + 3.0 * fx * fx * fx * fz * pb_y

                     + 7.5 * fx * fx * pa_zz * fz * pb_y

                     - 1.5 * pa_x * fx * pb_xy * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_xy * fx

                     - 1.5 * fx * fz * fga * pb_xxy

                     - 3.0 * pa_xzz * pb_xy * fz * fgb

                     + 7.5 * pa_x * fz * fx * fx * pb_xy

                     + 7.5 * fx * fx * fz * pb_xxy

                     + 18.0 * pa_xzz * fz * pb_xy * fx

                     + 18.0 * fx * pa_zz * fz * pb_xxy

                     - pa_x * fz * fga * pb_xxxy

                     + 6.0 * pa_x * fz * fx * pb_xxxy

                     + 14.0 * pa_xzz * fz * pb_xxxy);

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
    inline double fvec_xzz_xxxz_r_0(double fga,
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
                                    double pb_xxx,
                                    double pb_xxxz,
                                    double pb_xxz,
                                    double pb_xz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_z * fz * fgb

                     + 6.0 * fx * fx * fx * pa_z * fz

                     - 0.75 * fx * fx * fz * fgb * pb_z

                     - 0.75 * fx * fx * fz * fga * pb_z

                     - 3.0 * pa_xz * fx * pb_x * fz * fgb

                     - 1.5 * fx * pa_zz * fz * fgb * pb_z

                     + 3.0 * fx * fx * fx * fz * pb_z

                     + 15.0 * pa_xz * fz * fx * fx * pb_x

                     + 7.5 * fx * fx * pa_zz * fz * pb_z

                     + 15.0 * fx * fx * pa_z * fz * pb_xx

                     - 1.5 * pa_x * fx * pb_xz * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_xz * fx

                     - 1.5 * fx * fz * fga * pb_xxz

                     - 3.0 * pa_xzz * pb_xz * fz * fgb

                     + 7.5 * pa_x * fz * fx * fx * pb_xz

                     + 7.5 * fx * fx * fz * pb_xxz

                     + 18.0 * pa_xzz * fz * pb_xz * fx

                     + 12.0 * pa_xz * fz * fx * pb_xxx

                     + 18.0 * fx * pa_zz * fz * pb_xxz

                     - pa_x * fz * fga * pb_xxxz

                     + 6.0 * pa_x * fz * fx * pb_xxxz

                     + 14.0 * pa_xzz * fz * pb_xxxz);

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
    inline double fvec_xzz_xxyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xzz,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xx,
                                    double pb_xxyy,
                                    double pb_xyy,
                                    double pb_yy,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_x * fx * fx * fz * fgb

                     - 0.25 * pa_x * fz * fga * fx * fx

                     - 0.5 * fx * fx * pb_x * fz * fgb

                     - 0.5 * fx * fx * fz * fga * pb_x

                     - pa_xzz * fx * fz * fgb

                     + pa_x * fz * fx * fx * fx

                     - fx * pa_zz * pb_x * fz * fgb

                     + 2.0 * fx * fx * fx * fz * pb_x

                     + 2.5 * pa_xzz * fz * fx * fx

                     + 5.0 * fx * fx * pa_zz * fz * pb_x

                     - 0.5 * pa_x * fx * pb_xx * fz * fgb

                     - 0.5 * pa_x * fx * fz * fgb * pb_yy

                     - 0.5 * pa_x * fz * fga * pb_xx * fx

                     - 0.5 * pa_x * fz * fga * fx * pb_yy

                     - fx * fz * fga * pb_xyy

                     - pa_xzz * pb_xx * fz * fgb

                     - pa_xzz * fz * fgb * pb_yy

                     + 2.5 * pa_x * fz * fx * fx * pb_xx

                     + 2.5 * pa_x * fz * fx * fx * pb_yy

                     + 5.0 * fx * fx * fz * pb_xyy

                     + 6.0 * pa_xzz * fz * pb_xx * fx

                     + 6.0 * pa_xzz * fz * fx * pb_yy

                     + 12.0 * fx * pa_zz * fz * pb_xyy

                     - pa_x * fz * fga * pb_xxyy

                     + 6.0 * pa_x * fz * fx * pb_xxyy

                     + 14.0 * pa_xzz * fz * pb_xxyy);

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
    inline double fvec_xzz_xxyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
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
                                    double r_0_0)
    {
        return r_0_0 * (- pa_xz * fx * fz * fgb * pb_y

                     + 5.0 * pa_xz * fz * fx * fx * pb_y

                     + 10.0 * fx * fx * pa_z * fz * pb_xy

                     - 0.5 * pa_x * fx * fz * fgb * pb_yz

                     - 0.5 * pa_x * fz * fga * fx * pb_yz

                     - fx * fz * fga * pb_xyz

                     - pa_xzz * fz * fgb * pb_yz

                     + 2.5 * pa_x * fz * fx * fx * pb_yz

                     + 5.0 * fx * fx * fz * pb_xyz

                     + 6.0 * pa_xzz * fz * fx * pb_yz

                     + 12.0 * pa_xz * fz * fx * pb_xxy

                     + 12.0 * fx * pa_zz * fz * pb_xyz

                     - pa_x * fz * fga * pb_xxyz

                     + 6.0 * pa_x * fz * fx * pb_xxyz

                     + 14.0 * pa_xzz * fz * pb_xxyz);

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
    inline double fvec_xzz_xxzz_r_0(double fga,
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
                                    double pb_xxzz,
                                    double pb_xz,
                                    double pb_xzz,
                                    double pb_z,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- pa_x * fx * fx * fz * fgb

                     + 3.0 * pa_x * fx * fx * fx * fz

                     + 6.0 * fx * fx * fx * fz * pb_x

                     - 0.25 * pa_x * fz * fga * fx * fx

                     - 0.5 * fx * fx * pb_x * fz * fgb

                     - 0.5 * fx * fx * fz * fga * pb_x

                     - pa_xzz * fx * fz * fgb

                     - 2.0 * pa_xz * fx * fz * fgb * pb_z

                     - fx * pa_zz * pb_x * fz * fgb

                     + 2.5 * pa_xzz * fz * fx * fx

                     + 10.0 * pa_xz * fz * fx * fx * pb_z

                     + 7.5 * pa_x * fx * fx * fz * pb_xx

                     + 5.0 * fx * fx * pa_zz * fz * pb_x

                     + 20.0 * fx * fx * pa_z * fz * pb_xz

                     - 0.5 * pa_x * fx * pb_xx * fz * fgb

                     - 0.5 * pa_x * fx * fz * fgb * pb_zz

                     - 0.5 * pa_x * fz * fga * pb_xx * fx

                     - 0.5 * pa_x * fz * fga * fx * pb_zz

                     - fx * fz * fga * pb_xzz

                     - pa_xzz * pb_xx * fz * fgb

                     - pa_xzz * fz * fgb * pb_zz

                     + 2.5 * pa_x * fz * fx * fx * pb_zz

                     + 5.0 * fx * fx * fz * pb_xzz

                     + 6.0 * pa_xzz * fz * pb_xx * fx

                     + 6.0 * pa_xzz * fz * fx * pb_zz

                     + 24.0 * pa_xz * fz * fx * pb_xxz

                     + 12.0 * fx * pa_zz * fz * pb_xzz

                     - pa_x * fz * fga * pb_xxzz

                     + 6.0 * pa_x * fz * fx * pb_xxzz

                     + 14.0 * pa_xzz * fz * pb_xxzz);

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
    inline double fvec_xzz_xyyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xzz,
                                    double pa_zz,
                                    double pb_xy,
                                    double pb_xyyy,
                                    double pb_y,
                                    double pb_yyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pb_y * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pb_y

                     - 1.5 * fx * pa_zz * pb_y * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_y

                     + 7.5 * fx * fx * pa_zz * fz * pb_y

                     - 1.5 * pa_x * fx * pb_xy * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_xy * fx

                     - 0.5 * fx * fz * fga * pb_yyy

                     - 3.0 * pa_xzz * pb_xy * fz * fgb

                     + 7.5 * pa_x * fz * fx * fx * pb_xy

                     + 2.5 * fx * fx * fz * pb_yyy

                     + 18.0 * pa_xzz * fz * pb_xy * fx

                     + 6.0 * fx * pa_zz * fz * pb_yyy

                     - pa_x * fz * fga * pb_xyyy

                     + 6.0 * pa_x * fz * fx * pb_xyyy

                     + 14.0 * pa_xzz * fz * pb_xyyy);

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
    inline double fvec_xzz_xyyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
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
                                    double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * pa_z * fz * fgb

                     + 2.0 * fx * fx * fx * pa_z * fz

                     - 0.25 * fx * fx * fz * fgb * pb_z

                     - 0.25 * fx * fx * fz * fga * pb_z

                     - pa_xz * fx * pb_x * fz * fgb

                     - 0.5 * fx * pa_zz * fz * fgb * pb_z

                     + fx * fx * fx * fz * pb_z

                     + 5.0 * pa_xz * fz * fx * fx * pb_x

                     + 2.5 * fx * fx * pa_zz * fz * pb_z

                     + 5.0 * fx * fx * pa_z * fz * pb_yy

                     - 0.5 * pa_x * fx * pb_xz * fz * fgb

                     - 0.5 * pa_x * fz * fga * pb_xz * fx

                     - 0.5 * fx * fz * fga * pb_yyz

                     - pa_xzz * pb_xz * fz * fgb

                     + 2.5 * pa_x * fz * fx * fx * pb_xz

                     + 2.5 * fx * fx * fz * pb_yyz

                     + 6.0 * pa_xzz * fz * pb_xz * fx

                     + 12.0 * pa_xz * fz * fx * pb_xyy

                     + 6.0 * fx * pa_zz * fz * pb_yyz

                     - pa_x * fz * fga * pb_xyyz

                     + 6.0 * pa_x * fz * fx * pb_xyyz

                     + 14.0 * pa_xzz * fz * pb_xyyz);

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
    inline double fvec_xzz_xyzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
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
                                    double r_0_0)
    {
        return r_0_0 * (3.0 * fx * fx * fx * fz * pb_y

                     - 0.25 * fx * fx * pb_y * fz * fgb

                     - 0.25 * fx * fx * fz * fga * pb_y

                     - 0.5 * fx * pa_zz * pb_y * fz * fgb

                     + 7.5 * pa_x * fx * fx * fz * pb_xy

                     + 2.5 * fx * fx * pa_zz * fz * pb_y

                     + 10.0 * fx * fx * pa_z * fz * pb_yz

                     - 0.5 * pa_x * fx * pb_xy * fz * fgb

                     - 0.5 * pa_x * fz * fga * pb_xy * fx

                     - 0.5 * fx * fz * fga * pb_yzz

                     - pa_xzz * pb_xy * fz * fgb

                     + 2.5 * fx * fx * fz * pb_yzz

                     + 6.0 * pa_xzz * fz * pb_xy * fx

                     + 24.0 * pa_xz * fz * fx * pb_xyz

                     + 6.0 * fx * pa_zz * fz * pb_yzz

                     - pa_x * fz * fga * pb_xyzz

                     + 6.0 * pa_x * fz * fx * pb_xyzz

                     + 14.0 * pa_xzz * fz * pb_xyzz);

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
    inline double fvec_xzz_xzzz_r_0(double fga,
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
                                    double pb_xzzz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_z * fz * fgb

                     + 6.0 * fx * fx * fx * pa_z * fz

                     + 9.0 * fx * fx * fx * fz * pb_z

                     - 0.75 * fx * fx * pb_z * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pb_z

                     - 3.0 * pa_xz * fx * pb_x * fz * fgb

                     - 1.5 * fx * pa_zz * pb_z * fz * fgb

                     + 15.0 * pa_xz * fz * fx * fx * pb_x

                     + 22.5 * pa_x * fx * fx * fz * pb_xz

                     + 7.5 * fx * fx * pa_zz * fz * pb_z

                     + 15.0 * fx * fx * pa_z * fz * pb_zz

                     - 1.5 * pa_x * fx * pb_xz * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_xz * fx

                     - 0.5 * fx * fz * fga * pb_zzz

                     - 3.0 * pa_xzz * pb_xz * fz * fgb

                     + 2.5 * fx * fx * fz * pb_zzz

                     + 18.0 * pa_xzz * fz * pb_xz * fx

                     + 36.0 * pa_xz * fz * fx * pb_xzz

                     + 6.0 * fx * pa_zz * fz * pb_zzz

                     - pa_x * fz * fga * pb_xzzz

                     + 6.0 * pa_x * fz * fx * pb_xzzz

                     + 14.0 * pa_xzz * fz * pb_xzzz);

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
    inline double fvec_xzz_yyyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xzz,
                                    double pb_yy,
                                    double pb_yyyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fx * fz * fgb

                     - 0.75 * pa_x * fz * fga * fx * fx

                     - 3.0 * pa_xzz * fx * fz * fgb

                     + 3.0 * pa_x * fz * fx * fx * fx

                     + 7.5 * pa_xzz * fz * fx * fx

                     - 3.0 * pa_x * fx * pb_yy * fz * fgb

                     - 3.0 * pa_x * fz * fga * pb_yy * fx

                     - 6.0 * pa_xzz * pb_yy * fz * fgb

                     + 15.0 * pa_x * fz * fx * fx * pb_yy

                     + 36.0 * pa_xzz * fz * pb_yy * fx

                     - pa_x * fz * fga * pb_yyyy

                     + 6.0 * pa_x * fz * fx * pb_yyyy

                     + 14.0 * pa_xzz * fz * pb_yyyy);

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
    inline double fvec_xzz_yyyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pb_y,
                                    double pb_yyy,
                                    double pb_yyyz,
                                    double pb_yz,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xz * fx * pb_y * fz * fgb

                     + 15.0 * pa_xz * fz * fx * fx * pb_y

                     - 1.5 * pa_x * fx * pb_yz * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_yz * fx

                     - 3.0 * pa_xzz * pb_yz * fz * fgb

                     + 7.5 * pa_x * fz * fx * fx * pb_yz

                     + 18.0 * pa_xzz * fz * pb_yz * fx

                     + 12.0 * pa_xz * fz * fx * pb_yyy

                     - pa_x * fz * fga * pb_yyyz

                     + 6.0 * pa_x * fz * fx * pb_yyyz

                     + 14.0 * pa_xzz * fz * pb_yyyz);

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
    inline double fvec_xzz_yyzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_yyzz,
                                    double pb_z,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- pa_x * fx * fx * fz * fgb

                     + 3.0 * pa_x * fx * fx * fx * fz

                     - 0.25 * pa_x * fz * fga * fx * fx

                     - pa_xzz * fx * fz * fgb

                     - 2.0 * pa_xz * fx * fz * fgb * pb_z

                     + 2.5 * pa_xzz * fz * fx * fx

                     + 10.0 * pa_xz * fz * fx * fx * pb_z

                     + 7.5 * pa_x * fx * fx * fz * pb_yy

                     - 0.5 * pa_x * fx * pb_yy * fz * fgb

                     - 0.5 * pa_x * fx * fz * fgb * pb_zz

                     - 0.5 * pa_x * fz * fga * pb_yy * fx

                     - 0.5 * pa_x * fz * fga * fx * pb_zz

                     - pa_xzz * pb_yy * fz * fgb

                     - pa_xzz * fz * fgb * pb_zz

                     + 2.5 * pa_x * fz * fx * fx * pb_zz

                     + 6.0 * pa_xzz * fz * pb_yy * fx

                     + 6.0 * pa_xzz * fz * fx * pb_zz

                     + 24.0 * pa_xz * fz * fx * pb_yyz

                     - pa_x * fz * fga * pb_yyzz

                     + 6.0 * pa_x * fz * fx * pb_yyzz

                     + 14.0 * pa_xzz * fz * pb_yyzz);

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
    inline double fvec_xzz_yzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double pb_yzzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xz * fx * pb_y * fz * fgb

                     + 15.0 * pa_xz * fz * fx * fx * pb_y

                     + 22.5 * pa_x * fx * fx * fz * pb_yz

                     - 1.5 * pa_x * fx * pb_yz * fz * fgb

                     - 1.5 * pa_x * fz * fga * pb_yz * fx

                     - 3.0 * pa_xzz * pb_yz * fz * fgb

                     + 18.0 * pa_xzz * fz * pb_yz * fx

                     + 36.0 * pa_xz * fz * fx * pb_yzz

                     - pa_x * fz * fga * pb_yzzz

                     + 6.0 * pa_x * fz * fx * pb_yzzz

                     + 14.0 * pa_xzz * fz * pb_yzzz);

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
    inline double fvec_xzz_zzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_x,
                                    double pa_xz,
                                    double pa_xzz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double pb_zzzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_x * fx * fx * fz * fgb

                     + 15.0 * pa_x * fx * fx * fx * fz

                     - 0.75 * pa_x * fz * fga * fx * fx

                     - 3.0 * pa_xzz * fx * fz * fgb

                     - 12.0 * pa_xz * fx * pb_z * fz * fgb

                     + 7.5 * pa_xzz * fz * fx * fx

                     + 60.0 * pa_xz * fz * fx * fx * pb_z

                     + 45.0 * pa_x * fx * fx * fz * pb_zz

                     - 3.0 * pa_x * fx * pb_zz * fz * fgb

                     - 3.0 * pa_x * fz * fga * pb_zz * fx

                     - 6.0 * pa_xzz * pb_zz * fz * fgb

                     + 36.0 * pa_xzz * fz * pb_zz * fx

                     + 48.0 * pa_xz * fz * fx * pb_zzz

                     - pa_x * fz * fga * pb_zzzz

                     + 6.0 * pa_x * fz * fx * pb_zzzz

                     + 14.0 * pa_xzz * fz * pb_zzzz);

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
    inline double fvec_yyy_xxxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yyy,
                                    double pb_xx,
                                    double pb_xxxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_y * fx * fx * fz * fgb

                     - 2.25 * pa_y * fz * fga * fx * fx

                     - 3.0 * pa_yyy * fx * fz * fgb

                     + 9.0 * pa_y * fz * fx * fx * fx

                     + 7.5 * pa_yyy * fz * fx * fx

                     - 9.0 * pa_y * fx * pb_xx * fz * fgb

                     - 9.0 * pa_y * fz * fga * pb_xx * fx

                     - 6.0 * pa_yyy * pb_xx * fz * fgb

                     + 45.0 * pa_y * fz * fx * fx * pb_xx

                     + 36.0 * pa_yyy * fz * pb_xx * fx

                     - 3.0 * pa_y * fz * fga * pb_xxxx

                     + 18.0 * pa_y * fz * fx * pb_xxxx

                     + 14.0 * pa_yyy * fz * pb_xxxx);

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
    inline double fvec_yyy_xxxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_x,
                                    double pb_xxx,
                                    double pb_xxxy,
                                    double pb_xy,
                                    double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pb_x * fz * fgb

                     - 2.25 * fx * fx * fz * fga * pb_x

                     - 4.5 * pa_yy * fx * pb_x * fz * fgb

                     + 9.0 * fx * fx * fx * fz * pb_x

                     + 22.5 * pa_yy * fz * fx * fx * pb_x

                     - 4.5 * pa_y * fx * pb_xy * fz * fgb

                     - 4.5 * pa_y * fz * fga * pb_xy * fx

                     - 1.5 * fx * fz * fga * pb_xxx

                     - 3.0 * pa_yyy * pb_xy * fz * fgb

                     + 22.5 * pa_y * fz * fx * fx * pb_xy

                     + 7.5 * fx * fx * fz * pb_xxx

                     + 18.0 * pa_yyy * fz * pb_xy * fx

                     + 18.0 * pa_yy * fz * fx * pb_xxx

                     - 3.0 * pa_y * fz * fga * pb_xxxy

                     + 18.0 * pa_y * fz * fx * pb_xxxy

                     + 14.0 * pa_yyy * fz * pb_xxxy);

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
    inline double fvec_yyy_xxxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yyy,
                                    double pb_xxxz,
                                    double pb_xz,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_y * fx * pb_xz * fz * fgb

                     - 4.5 * pa_y * fz * fga * pb_xz * fx

                     - 3.0 * pa_yyy * pb_xz * fz * fgb

                     + 22.5 * pa_y * fz * fx * fx * pb_xz

                     + 18.0 * pa_yyy * fz * pb_xz * fx

                     - 3.0 * pa_y * fz * fga * pb_xxxz

                     + 18.0 * pa_y * fz * fx * pb_xxxz

                     + 14.0 * pa_yyy * fz * pb_xxxz);

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
    inline double fvec_yyy_xxyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_xx,
                                    double pb_xxy,
                                    double pb_xxyy,
                                    double pb_y,
                                    double pb_yy,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_y * fx * fx * fz * fgb

                     + 9.0 * pa_y * fx * fx * fx * fz

                     - 0.75 * pa_y * fz * fga * fx * fx

                     - 1.5 * fx * fx * fz * fgb * pb_y

                     - 1.5 * fx * fx * fz * fga * pb_y

                     - pa_yyy * fx * fz * fgb

                     - 3.0 * pa_yy * fx * fz * fgb * pb_y

                     + 6.0 * fx * fx * fx * fz * pb_y

                     + 2.5 * pa_yyy * fz * fx * fx

                     + 15.0 * pa_yy * fz * fx * fx * pb_y

                     + 22.5 * pa_y * fx * fx * fz * pb_xx

                     - 1.5 * pa_y * fx * pb_xx * fz * fgb

                     - 1.5 * pa_y * fx * fz * fgb * pb_yy

                     - 1.5 * pa_y * fz * fga * pb_xx * fx

                     - 1.5 * pa_y * fz * fga * fx * pb_yy

                     - 3.0 * fx * fz * fga * pb_xxy

                     - pa_yyy * pb_xx * fz * fgb

                     - pa_yyy * fz * fgb * pb_yy

                     + 7.5 * pa_y * fz * fx * fx * pb_yy

                     + 15.0 * fx * fx * fz * pb_xxy

                     + 6.0 * pa_yyy * fz * pb_xx * fx

                     + 6.0 * pa_yyy * fz * fx * pb_yy

                     + 36.0 * pa_yy * fz * fx * pb_xxy

                     - 3.0 * pa_y * fz * fga * pb_xxyy

                     + 18.0 * pa_y * fz * fx * pb_xxyy

                     + 14.0 * pa_yyy * fz * pb_xxyy);

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
    inline double fvec_yyy_xxyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_xxyz,
                                    double pb_xxz,
                                    double pb_yz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb * pb_z

                     - 0.75 * fx * fx * fz * fga * pb_z

                     - 1.5 * pa_yy * fx * fz * fgb * pb_z

                     + 3.0 * fx * fx * fx * fz * pb_z

                     + 7.5 * pa_yy * fz * fx * fx * pb_z

                     - 1.5 * pa_y * fx * fz * fgb * pb_yz

                     - 1.5 * pa_y * fz * fga * fx * pb_yz

                     - 1.5 * fx * fz * fga * pb_xxz

                     - pa_yyy * fz * fgb * pb_yz

                     + 7.5 * pa_y * fz * fx * fx * pb_yz

                     + 7.5 * fx * fx * fz * pb_xxz

                     + 6.0 * pa_yyy * fz * fx * pb_yz

                     + 18.0 * pa_yy * fz * fx * pb_xxz

                     - 3.0 * pa_y * fz * fga * pb_xxyz

                     + 18.0 * pa_y * fz * fx * pb_xxyz

                     + 14.0 * pa_yyy * fz * pb_xxyz);

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
    inline double fvec_yyy_xxzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yyy,
                                    double pb_xx,
                                    double pb_xxzz,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_y * fx * fx * fz * fgb

                     - 0.75 * pa_y * fz * fga * fx * fx

                     - pa_yyy * fx * fz * fgb

                     + 3.0 * pa_y * fz * fx * fx * fx

                     + 2.5 * pa_yyy * fz * fx * fx

                     - 1.5 * pa_y * fx * pb_xx * fz * fgb

                     - 1.5 * pa_y * fx * fz * fgb * pb_zz

                     - 1.5 * pa_y * fz * fga * pb_xx * fx

                     - 1.5 * pa_y * fz * fga * fx * pb_zz

                     - pa_yyy * pb_xx * fz * fgb

                     - pa_yyy * fz * fgb * pb_zz

                     + 7.5 * pa_y * fz * fx * fx * pb_xx

                     + 7.5 * pa_y * fz * fx * fx * pb_zz

                     + 6.0 * pa_yyy * fz * pb_xx * fx

                     + 6.0 * pa_yyy * fz * fx * pb_zz

                     - 3.0 * pa_y * fz * fga * pb_xxzz

                     + 18.0 * pa_y * fz * fx * pb_xxzz

                     + 14.0 * pa_yyy * fz * pb_xxzz);

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
    inline double fvec_yyy_xyyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyy,
                                    double pb_xyyy,
                                    double r_0_0)
    {
        return r_0_0 * (15.0 * fx * fx * fx * fz * pb_x

                     - 2.25 * fx * fx * pb_x * fz * fgb

                     - 2.25 * fx * fx * fz * fga * pb_x

                     - 4.5 * pa_yy * fx * pb_x * fz * fgb

                     + 22.5 * pa_yy * fz * fx * fx * pb_x

                     + 67.5 * pa_y * fx * fx * fz * pb_xy

                     - 4.5 * pa_y * fx * pb_xy * fz * fgb

                     - 4.5 * pa_y * fz * fga * pb_xy * fx

                     - 4.5 * fx * fz * fga * pb_xyy

                     - 3.0 * pa_yyy * pb_xy * fz * fgb

                     + 22.5 * fx * fx * fz * pb_xyy

                     + 18.0 * pa_yyy * fz * pb_xy * fx

                     + 54.0 * pa_yy * fz * fx * pb_xyy

                     - 3.0 * pa_y * fz * fga * pb_xyyy

                     + 18.0 * pa_y * fz * fx * pb_xyyy

                     + 14.0 * pa_yyy * fz * pb_xyyy);

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
    inline double fvec_yyy_xyyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_xyyz,
                                    double pb_xyz,
                                    double pb_xz,
                                    double r_0_0)
    {
        return r_0_0 * (22.5 * pa_y * fx * fx * fz * pb_xz

                     - 1.5 * pa_y * fx * pb_xz * fz * fgb

                     - 1.5 * pa_y * fz * fga * pb_xz * fx

                     - 3.0 * fx * fz * fga * pb_xyz

                     - pa_yyy * pb_xz * fz * fgb

                     + 15.0 * fx * fx * fz * pb_xyz

                     + 6.0 * pa_yyy * fz * pb_xz * fx

                     + 36.0 * pa_yy * fz * fx * pb_xyz

                     - 3.0 * pa_y * fz * fga * pb_xyyz

                     + 18.0 * pa_y * fz * fx * pb_xyyz

                     + 14.0 * pa_yyy * fz * pb_xyyz);

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
    inline double fvec_yyy_xyzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyzz,
                                    double pb_xzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pb_x * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pb_x

                     - 1.5 * pa_yy * fx * pb_x * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_x

                     + 7.5 * pa_yy * fz * fx * fx * pb_x

                     - 1.5 * pa_y * fx * pb_xy * fz * fgb

                     - 1.5 * pa_y * fz * fga * pb_xy * fx

                     - 1.5 * fx * fz * fga * pb_xzz

                     - pa_yyy * pb_xy * fz * fgb

                     + 7.5 * pa_y * fz * fx * fx * pb_xy

                     + 7.5 * fx * fx * fz * pb_xzz

                     + 6.0 * pa_yyy * fz * pb_xy * fx

                     + 18.0 * pa_yy * fz * fx * pb_xzz

                     - 3.0 * pa_y * fz * fga * pb_xyzz

                     + 18.0 * pa_y * fz * fx * pb_xyzz

                     + 14.0 * pa_yyy * fz * pb_xyzz);

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
    inline double fvec_yyy_xzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yyy,
                                    double pb_xz,
                                    double pb_xzzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_y * fx * pb_xz * fz * fgb

                     - 4.5 * pa_y * fz * fga * pb_xz * fx

                     - 3.0 * pa_yyy * pb_xz * fz * fgb

                     + 22.5 * pa_y * fz * fx * fx * pb_xz

                     + 18.0 * pa_yyy * fz * pb_xz * fx

                     - 3.0 * pa_y * fz * fga * pb_xzzz

                     + 18.0 * pa_y * fz * fx * pb_xzzz

                     + 14.0 * pa_yyy * fz * pb_xzzz);

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
    inline double fvec_yyy_yyyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double pb_yyyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 13.5 * pa_y * fx * fx * fz * fgb

                     + 45.0 * pa_y * fx * fx * fx * fz

                     + 60.0 * fx * fx * fx * fz * pb_y

                     - 2.25 * pa_y * fz * fga * fx * fx

                     - 9.0 * fx * fx * pb_y * fz * fgb

                     - 9.0 * fx * fx * fz * fga * pb_y

                     - 3.0 * pa_yyy * fx * fz * fgb

                     - 18.0 * pa_yy * fx * pb_y * fz * fgb

                     + 7.5 * pa_yyy * fz * fx * fx

                     + 90.0 * pa_yy * fz * fx * fx * pb_y

                     + 135.0 * pa_y * fx * fx * fz * pb_yy

                     - 9.0 * pa_y * fx * pb_yy * fz * fgb

                     - 9.0 * pa_y * fz * fga * pb_yy * fx

                     - 6.0 * fx * fz * fga * pb_yyy

                     - 6.0 * pa_yyy * pb_yy * fz * fgb

                     + 30.0 * fx * fx * fz * pb_yyy

                     + 36.0 * pa_yyy * fz * pb_yy * fx

                     + 72.0 * pa_yy * fz * fx * pb_yyy

                     - 3.0 * pa_y * fz * fga * pb_yyyy

                     + 18.0 * pa_y * fz * fx * pb_yyyy

                     + 14.0 * pa_yyy * fz * pb_yyyy);

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
    inline double fvec_yyy_yyyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_yyyz,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (15.0 * fx * fx * fx * fz * pb_z

                     - 2.25 * fx * fx * fz * fgb * pb_z

                     - 2.25 * fx * fx * fz * fga * pb_z

                     - 4.5 * pa_yy * fx * fz * fgb * pb_z

                     + 22.5 * pa_yy * fz * fx * fx * pb_z

                     + 67.5 * pa_y * fx * fx * fz * pb_yz

                     - 4.5 * pa_y * fx * pb_yz * fz * fgb

                     - 4.5 * pa_y * fz * fga * pb_yz * fx

                     - 4.5 * fx * fz * fga * pb_yyz

                     - 3.0 * pa_yyy * pb_yz * fz * fgb

                     + 22.5 * fx * fx * fz * pb_yyz

                     + 18.0 * pa_yyy * fz * pb_yz * fx

                     + 54.0 * pa_yy * fz * fx * pb_yyz

                     - 3.0 * pa_y * fz * fga * pb_yyyz

                     + 18.0 * pa_y * fz * fx * pb_yyyz

                     + 14.0 * pa_yyy * fz * pb_yyyz);

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
    inline double fvec_yyy_yyzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyzz,
                                    double pb_yzz,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_y * fx * fx * fz * fgb

                     + 9.0 * pa_y * fx * fx * fx * fz

                     - 0.75 * pa_y * fz * fga * fx * fx

                     - 1.5 * fx * fx * pb_y * fz * fgb

                     - 1.5 * fx * fx * fz * fga * pb_y

                     - pa_yyy * fx * fz * fgb

                     - 3.0 * pa_yy * fx * pb_y * fz * fgb

                     + 6.0 * fx * fx * fx * fz * pb_y

                     + 2.5 * pa_yyy * fz * fx * fx

                     + 15.0 * pa_yy * fz * fx * fx * pb_y

                     + 22.5 * pa_y * fx * fx * fz * pb_zz

                     - 1.5 * pa_y * fx * pb_yy * fz * fgb

                     - 1.5 * pa_y * fx * fz * fgb * pb_zz

                     - 1.5 * pa_y * fz * fga * pb_yy * fx

                     - 1.5 * pa_y * fz * fga * fx * pb_zz

                     - 3.0 * fx * fz * fga * pb_yzz

                     - pa_yyy * pb_yy * fz * fgb

                     - pa_yyy * fz * fgb * pb_zz

                     + 7.5 * pa_y * fz * fx * fx * pb_yy

                     + 15.0 * fx * fx * fz * pb_yzz

                     + 6.0 * pa_yyy * fz * pb_yy * fx

                     + 6.0 * pa_yyy * fz * fx * pb_zz

                     + 36.0 * pa_yy * fz * fx * pb_yzz

                     - 3.0 * pa_y * fz * fga * pb_yyzz

                     + 18.0 * pa_y * fz * fx * pb_yyzz

                     + 14.0 * pa_yyy * fz * pb_yyzz);

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
    inline double fvec_yyy_yzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yy,
                                    double pa_yyy,
                                    double pb_yz,
                                    double pb_yzzz,
                                    double pb_z,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pb_z * fz * fgb

                     - 2.25 * fx * fx * fz * fga * pb_z

                     - 4.5 * pa_yy * fx * pb_z * fz * fgb

                     + 9.0 * fx * fx * fx * fz * pb_z

                     + 22.5 * pa_yy * fz * fx * fx * pb_z

                     - 4.5 * pa_y * fx * pb_yz * fz * fgb

                     - 4.5 * pa_y * fz * fga * pb_yz * fx

                     - 1.5 * fx * fz * fga * pb_zzz

                     - 3.0 * pa_yyy * pb_yz * fz * fgb

                     + 22.5 * pa_y * fz * fx * fx * pb_yz

                     + 7.5 * fx * fx * fz * pb_zzz

                     + 18.0 * pa_yyy * fz * pb_yz * fx

                     + 18.0 * pa_yy * fz * fx * pb_zzz

                     - 3.0 * pa_y * fz * fga * pb_yzzz

                     + 18.0 * pa_y * fz * fx * pb_yzzz

                     + 14.0 * pa_yyy * fz * pb_yzzz);

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
    inline double fvec_yyy_zzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yyy,
                                    double pb_zz,
                                    double pb_zzzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_y * fx * fx * fz * fgb

                     - 2.25 * pa_y * fz * fga * fx * fx

                     - 3.0 * pa_yyy * fx * fz * fgb

                     + 9.0 * pa_y * fz * fx * fx * fx

                     + 7.5 * pa_yyy * fz * fx * fx

                     - 9.0 * pa_y * fx * pb_zz * fz * fgb

                     - 9.0 * pa_y * fz * fga * pb_zz * fx

                     - 6.0 * pa_yyy * pb_zz * fz * fgb

                     + 45.0 * pa_y * fz * fx * fx * pb_zz

                     + 36.0 * pa_yyy * fz * pb_zz * fx

                     - 3.0 * pa_y * fz * fga * pb_zzzz

                     + 18.0 * pa_y * fz * fx * pb_zzzz

                     + 14.0 * pa_yyy * fz * pb_zzzz);

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
    inline double fvec_yyz_xxxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yyz,
                                    double pa_z,
                                    double pb_xx,
                                    double pb_xxxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_z * fz * fgb

                     - 0.75 * fz * fga * pa_z * fx * fx

                     - 3.0 * pa_yyz * fx * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_z

                     + 7.5 * pa_yyz * fz * fx * fx

                     - 3.0 * fx * pa_z * pb_xx * fz * fgb

                     - 3.0 * fz * fga * pa_z * pb_xx * fx

                     - 6.0 * pa_yyz * pb_xx * fz * fgb

                     + 15.0 * fx * fx * fz * pa_z * pb_xx

                     + 36.0 * pa_yyz * fz * pb_xx * fx

                     - fz * fga * pa_z * pb_xxxx

                     + 6.0 * fx * fz * pa_z * pb_xxxx

                     + 14.0 * pa_yyz * fz * pb_xxxx);

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
    inline double fvec_yyz_xxxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yyz,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xxx,
                                    double pb_xxxy,
                                    double pb_xy,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_yz * fx * pb_x * fz * fgb

                     + 15.0 * pa_yz * fx * fx * fz * pb_x

                     - 1.5 * fx * pa_z * pb_xy * fz * fgb

                     - 1.5 * fz * fga * pa_z * pb_xy * fx

                     - 3.0 * pa_yyz * pb_xy * fz * fgb

                     + 7.5 * fx * fx * fz * pa_z * pb_xy

                     + 18.0 * pa_yyz * fz * pb_xy * fx

                     + 12.0 * pa_yz * fx * fz * pb_xxx

                     - fz * fga * pa_z * pb_xxxy

                     + 6.0 * fx * fz * pa_z * pb_xxxy

                     + 14.0 * pa_yyz * fz * pb_xxxy);

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
    inline double fvec_yyz_xxxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xxx,
                                    double pb_xxxz,
                                    double pb_xz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pb_x * fz * fgb

                     - 0.75 * fz * fga * fx * fx * pb_x

                     - 1.5 * pa_yy * fx * pb_x * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_x

                     + 7.5 * pa_yy * fz * fx * fx * pb_x

                     - 1.5 * fx * pa_z * pb_xz * fz * fgb

                     - 1.5 * fz * fga * pa_z * pb_xz * fx

                     - 0.5 * fz * fga * fx * pb_xxx

                     - 3.0 * pa_yyz * pb_xz * fz * fgb

                     + 7.5 * fx * fx * fz * pa_z * pb_xz

                     + 2.5 * fx * fx * fz * pb_xxx

                     + 18.0 * pa_yyz * fz * pb_xz * fx

                     + 6.0 * pa_yy * fz * fx * pb_xxx

                     - fz * fga * pa_z * pb_xxxz

                     + 6.0 * fx * fz * pa_z * pb_xxxz

                     + 14.0 * pa_yyz * fz * pb_xxxz);

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
    inline double fvec_yyz_xxyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
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
        return r_0_0 * (- fx * fx * pa_z * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_z

                     - 0.25 * fz * fga * pa_z * fx * fx

                     - pa_yyz * fx * fz * fgb

                     - 2.0 * pa_yz * fx * fz * fgb * pb_y

                     + 2.5 * pa_yyz * fz * fx * fx

                     + 10.0 * pa_yz * fx * fx * fz * pb_y

                     + 7.5 * fx * fx * fz * pa_z * pb_xx

                     - 0.5 * fx * pa_z * pb_xx * fz * fgb

                     - 0.5 * fx * pa_z * fz * fgb * pb_yy

                     - 0.5 * fz * fga * pa_z * pb_xx * fx

                     - 0.5 * fz * fga * pa_z * fx * pb_yy

                     - pa_yyz * pb_xx * fz * fgb

                     - pa_yyz * fz * fgb * pb_yy

                     + 2.5 * fx * fx * fz * pa_z * pb_yy

                     + 6.0 * pa_yyz * fz * pb_xx * fx

                     + 6.0 * pa_yyz * fz * fx * pb_yy

                     + 24.0 * pa_yz * fx * fz * pb_xxy

                     - fz * fga * pa_z * pb_xxyy

                     + 6.0 * fx * fz * pa_z * pb_xxyy

                     + 14.0 * pa_yyz * fz * pb_xxyy);

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
    inline double fvec_yyz_xxyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
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
                                    double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_y * fx * fx * fz * fgb

                     + 2.0 * pa_y * fx * fx * fx * fz

                     - 0.25 * fx * fx * fz * fgb * pb_y

                     - 0.25 * fz * fga * fx * fx * pb_y

                     - 0.5 * pa_yy * fx * fz * fgb * pb_y

                     - pa_yz * fx * fz * fgb * pb_z

                     + fx * fx * fx * fz * pb_y

                     + 2.5 * pa_yy * fz * fx * fx * pb_y

                     + 5.0 * pa_yz * fx * fx * fz * pb_z

                     + 5.0 * pa_y * fx * fx * fz * pb_xx

                     - 0.5 * fx * pa_z * fz * fgb * pb_yz

                     - 0.5 * fz * fga * pa_z * fx * pb_yz

                     - 0.5 * fz * fga * fx * pb_xxy

                     - pa_yyz * fz * fgb * pb_yz

                     + 2.5 * fx * fx * fz * pa_z * pb_yz

                     + 2.5 * fx * fx * fz * pb_xxy

                     + 6.0 * pa_yyz * fz * fx * pb_yz

                     + 6.0 * pa_yy * fz * fx * pb_xxy

                     + 12.0 * pa_yz * fx * fz * pb_xxz

                     - fz * fga * pa_z * pb_xxyz

                     + 6.0 * fx * fz * pa_z * pb_xxyz

                     + 14.0 * pa_yyz * fz * pb_xxyz);

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
    inline double fvec_yyz_xxzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_z,
                                    double pb_xx,
                                    double pb_xxz,
                                    double pb_xxzz,
                                    double pb_z,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * pa_z * fz * fgb

                     - 0.5 * fx * fx * fz * fgb * pb_z

                     - 0.25 * fz * fga * pa_z * fx * fx

                     - 0.5 * fz * fga * fx * fx * pb_z

                     - pa_yyz * fx * fz * fgb

                     - pa_yy * fx * fz * fgb * pb_z

                     + fx * fx * fx * fz * pa_z

                     + 2.0 * fx * fx * fx * fz * pb_z

                     + 2.5 * pa_yyz * fz * fx * fx

                     + 5.0 * pa_yy * fz * fx * fx * pb_z

                     - 0.5 * fx * pa_z * pb_xx * fz * fgb

                     - 0.5 * fx * pa_z * fz * fgb * pb_zz

                     - 0.5 * fz * fga * pa_z * pb_xx * fx

                     - 0.5 * fz * fga * pa_z * fx * pb_zz

                     - fz * fga * fx * pb_xxz

                     - pa_yyz * pb_xx * fz * fgb

                     - pa_yyz * fz * fgb * pb_zz

                     + 2.5 * fx * fx * fz * pa_z * pb_xx

                     + 2.5 * fx * fx * fz * pa_z * pb_zz

                     + 5.0 * fx * fx * fz * pb_xxz

                     + 6.0 * pa_yyz * fz * pb_xx * fx

                     + 6.0 * pa_yyz * fz * fx * pb_zz

                     + 12.0 * pa_yy * fz * fx * pb_xxz

                     - fz * fga * pa_z * pb_xxzz

                     + 6.0 * fx * fz * pa_z * pb_xxzz

                     + 14.0 * pa_yyz * fz * pb_xxzz);

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
    inline double fvec_yyz_xyyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yyz,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyy,
                                    double pb_xyyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_yz * fx * pb_x * fz * fgb

                     + 15.0 * pa_yz * fx * fx * fz * pb_x

                     + 22.5 * fx * fx * fz * pa_z * pb_xy

                     - 1.5 * fx * pa_z * pb_xy * fz * fgb

                     - 1.5 * fz * fga * pa_z * pb_xy * fx

                     - 3.0 * pa_yyz * pb_xy * fz * fgb

                     + 18.0 * pa_yyz * fz * pb_xy * fx

                     + 36.0 * pa_yz * fx * fz * pb_xyy

                     - fz * fga * pa_z * pb_xyyy

                     + 6.0 * fx * fz * pa_z * pb_xyyy

                     + 14.0 * pa_yyz * fz * pb_xyyy);

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
    inline double fvec_yyz_xyyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
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
                                    double r_0_0)
    {
        return r_0_0 * (3.0 * fx * fx * fx * fz * pb_x

                     - 0.25 * fx * fx * pb_x * fz * fgb

                     - 0.25 * fz * fga * fx * fx * pb_x

                     - 0.5 * pa_yy * fx * pb_x * fz * fgb

                     + 2.5 * pa_yy * fz * fx * fx * pb_x

                     + 10.0 * pa_y * fx * fx * fz * pb_xy

                     + 7.5 * fx * fx * fz * pa_z * pb_xz

                     - 0.5 * fx * pa_z * pb_xz * fz * fgb

                     - 0.5 * fz * fga * pa_z * pb_xz * fx

                     - 0.5 * fz * fga * fx * pb_xyy

                     - pa_yyz * pb_xz * fz * fgb

                     + 2.5 * fx * fx * fz * pb_xyy

                     + 6.0 * pa_yyz * fz * pb_xz * fx

                     + 6.0 * pa_yy * fz * fx * pb_xyy

                     + 24.0 * pa_yz * fx * fz * pb_xyz

                     - fz * fga * pa_z * pb_xyyz

                     + 6.0 * fx * fz * pa_z * pb_xyyz

                     + 14.0 * pa_yyz * fz * pb_xyyz);

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
    inline double fvec_yyz_xyzz_r_0(double fga,
                                    double fgb,
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
                                    double pb_xyzz,
                                    double pb_xz,
                                    double pb_xzz,
                                    double r_0_0)
    {
        return r_0_0 * (- pa_yz * fx * pb_x * fz * fgb

                     + 5.0 * pa_yz * fx * fx * fz * pb_x

                     + 10.0 * pa_y * fx * fx * fz * pb_xz

                     - 0.5 * fx * pa_z * pb_xy * fz * fgb

                     - 0.5 * fz * fga * pa_z * pb_xy * fx

                     - fz * fga * fx * pb_xyz

                     - pa_yyz * pb_xy * fz * fgb

                     + 2.5 * fx * fx * fz * pa_z * pb_xy

                     + 5.0 * fx * fx * fz * pb_xyz

                     + 6.0 * pa_yyz * fz * pb_xy * fx

                     + 12.0 * pa_yy * fz * fx * pb_xyz

                     + 12.0 * pa_yz * fx * fz * pb_xzz

                     - fz * fga * pa_z * pb_xyzz

                     + 6.0 * fx * fz * pa_z * pb_xyzz

                     + 14.0 * pa_yyz * fz * pb_xyzz);

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
    inline double fvec_yyz_xzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_z,
                                    double pb_x,
                                    double pb_xz,
                                    double pb_xzz,
                                    double pb_xzzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pb_x * fz * fgb

                     - 0.75 * fz * fga * fx * fx * pb_x

                     - 1.5 * pa_yy * fx * pb_x * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_x

                     + 7.5 * pa_yy * fz * fx * fx * pb_x

                     - 1.5 * fx * pa_z * pb_xz * fz * fgb

                     - 1.5 * fz * fga * pa_z * pb_xz * fx

                     - 1.5 * fz * fga * fx * pb_xzz

                     - 3.0 * pa_yyz * pb_xz * fz * fgb

                     + 7.5 * fx * fx * fz * pa_z * pb_xz

                     + 7.5 * fx * fx * fz * pb_xzz

                     + 18.0 * pa_yyz * fz * pb_xz * fx

                     + 18.0 * pa_yy * fz * fx * pb_xzz

                     - fz * fga * pa_z * pb_xzzz

                     + 6.0 * fx * fz * pa_z * pb_xzzz

                     + 14.0 * pa_yyz * fz * pb_xzzz);

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
    inline double fvec_yyz_yyyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yyz,
                                    double pa_yz,
                                    double pa_z,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double pb_yyyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * fx * fx * pa_z * fz * fgb

                     + 15.0 * fx * fx * fx * fz * pa_z

                     - 0.75 * fz * fga * pa_z * fx * fx

                     - 3.0 * pa_yyz * fx * fz * fgb

                     - 12.0 * pa_yz * fx * pb_y * fz * fgb

                     + 7.5 * pa_yyz * fz * fx * fx

                     + 60.0 * pa_yz * fx * fx * fz * pb_y

                     + 45.0 * fx * fx * fz * pa_z * pb_yy

                     - 3.0 * fx * pa_z * pb_yy * fz * fgb

                     - 3.0 * fz * fga * pa_z * pb_yy * fx

                     - 6.0 * pa_yyz * pb_yy * fz * fgb

                     + 36.0 * pa_yyz * fz * pb_yy * fx

                     + 48.0 * pa_yz * fx * fz * pb_yyy

                     - fz * fga * pa_z * pb_yyyy

                     + 6.0 * fx * fz * pa_z * pb_yyyy

                     + 14.0 * pa_yyz * fz * pb_yyyy);

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
    inline double fvec_yyz_yyyz_r_0(double fga,
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
                                    double pb_yyy,
                                    double pb_yyyz,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_y * fx * fx * fz * fgb

                     + 6.0 * pa_y * fx * fx * fx * fz

                     + 9.0 * fx * fx * fx * fz * pb_y

                     - 0.75 * fx * fx * pb_y * fz * fgb

                     - 0.75 * fz * fga * fx * fx * pb_y

                     - 1.5 * pa_yy * fx * pb_y * fz * fgb

                     - 3.0 * pa_yz * fx * fz * fgb * pb_z

                     + 7.5 * pa_yy * fz * fx * fx * pb_y

                     + 15.0 * pa_yz * fx * fx * fz * pb_z

                     + 15.0 * pa_y * fx * fx * fz * pb_yy

                     + 22.5 * fx * fx * fz * pa_z * pb_yz

                     - 1.5 * fx * pa_z * pb_yz * fz * fgb

                     - 1.5 * fz * fga * pa_z * pb_yz * fx

                     - 0.5 * fz * fga * fx * pb_yyy

                     - 3.0 * pa_yyz * pb_yz * fz * fgb

                     + 2.5 * fx * fx * fz * pb_yyy

                     + 18.0 * pa_yyz * fz * pb_yz * fx

                     + 6.0 * pa_yy * fz * fx * pb_yyy

                     + 36.0 * pa_yz * fx * fz * pb_yyz

                     - fz * fga * pa_z * pb_yyyz

                     + 6.0 * fx * fz * pa_z * pb_yyyz

                     + 14.0 * pa_yyz * fz * pb_yyyz);

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
    inline double fvec_yyz_yyzz_r_0(double fga,
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
                                    double pb_yyzz,
                                    double pb_yz,
                                    double pb_yzz,
                                    double pb_z,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- fx * fx * pa_z * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_z

                     + 6.0 * fx * fx * fx * fz * pb_z

                     - 0.5 * fx * fx * fz * fgb * pb_z

                     - 0.25 * fz * fga * pa_z * fx * fx

                     - 0.5 * fz * fga * fx * fx * pb_z

                     - pa_yyz * fx * fz * fgb

                     - pa_yy * fx * fz * fgb * pb_z

                     - 2.0 * pa_yz * fx * pb_y * fz * fgb

                     + 2.5 * pa_yyz * fz * fx * fx

                     + 5.0 * pa_yy * fz * fx * fx * pb_z

                     + 10.0 * pa_yz * fx * fx * fz * pb_y

                     + 20.0 * pa_y * fx * fx * fz * pb_yz

                     + 7.5 * fx * fx * fz * pa_z * pb_zz

                     - 0.5 * fx * pa_z * pb_yy * fz * fgb

                     - 0.5 * fx * pa_z * fz * fgb * pb_zz

                     - 0.5 * fz * fga * pa_z * pb_yy * fx

                     - 0.5 * fz * fga * pa_z * fx * pb_zz

                     - fz * fga * fx * pb_yyz

                     - pa_yyz * pb_yy * fz * fgb

                     - pa_yyz * fz * fgb * pb_zz

                     + 2.5 * fx * fx * fz * pa_z * pb_yy

                     + 5.0 * fx * fx * fz * pb_yyz

                     + 6.0 * pa_yyz * fz * pb_yy * fx

                     + 6.0 * pa_yyz * fz * fx * pb_zz

                     + 12.0 * pa_yy * fz * fx * pb_yyz

                     + 24.0 * pa_yz * fx * fz * pb_yzz

                     - fz * fga * pa_z * pb_yyzz

                     + 6.0 * fx * fz * pa_z * pb_yyzz

                     + 14.0 * pa_yyz * fz * pb_yyzz);

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
    inline double fvec_yyz_yzzz_r_0(double fga,
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
                                    double pb_yzzz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_y * fx * fx * fz * fgb

                     + 6.0 * pa_y * fx * fx * fx * fz

                     - 0.75 * fx * fx * pb_y * fz * fgb

                     - 0.75 * fz * fga * fx * fx * pb_y

                     - 1.5 * pa_yy * fx * pb_y * fz * fgb

                     - 3.0 * pa_yz * fx * pb_z * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_y

                     + 7.5 * pa_yy * fz * fx * fx * pb_y

                     + 15.0 * pa_yz * fx * fx * fz * pb_z

                     + 15.0 * pa_y * fx * fx * fz * pb_zz

                     - 1.5 * fx * pa_z * pb_yz * fz * fgb

                     - 1.5 * fz * fga * pa_z * pb_yz * fx

                     - 1.5 * fz * fga * fx * pb_yzz

                     - 3.0 * pa_yyz * pb_yz * fz * fgb

                     + 7.5 * fx * fx * fz * pa_z * pb_yz

                     + 7.5 * fx * fx * fz * pb_yzz

                     + 18.0 * pa_yyz * fz * pb_yz * fx

                     + 18.0 * pa_yy * fz * fx * pb_yzz

                     + 12.0 * pa_yz * fx * fz * pb_zzz

                     - fz * fga * pa_z * pb_yzzz

                     + 6.0 * fx * fz * pa_z * pb_yzzz

                     + 14.0 * pa_yyz * fz * pb_yzzz);

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
    inline double fvec_yyz_zzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_yy,
                                    double pa_yyz,
                                    double pa_z,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double pb_zzzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_z * fz * fgb

                     - 3.0 * fx * fx * pb_z * fz * fgb

                     - 0.75 * fz * fga * pa_z * fx * fx

                     - 3.0 * fz * fga * fx * fx * pb_z

                     - 3.0 * pa_yyz * fx * fz * fgb

                     - 6.0 * pa_yy * fx * pb_z * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pa_z

                     + 12.0 * fx * fx * fx * fz * pb_z

                     + 7.5 * pa_yyz * fz * fx * fx

                     + 30.0 * pa_yy * fz * fx * fx * pb_z

                     - 3.0 * fx * pa_z * pb_zz * fz * fgb

                     - 3.0 * fz * fga * pa_z * pb_zz * fx

                     - 2.0 * fz * fga * fx * pb_zzz

                     - 6.0 * pa_yyz * pb_zz * fz * fgb

                     + 15.0 * fx * fx * fz * pa_z * pb_zz

                     + 10.0 * fx * fx * fz * pb_zzz

                     + 36.0 * pa_yyz * fz * pb_zz * fx

                     + 24.0 * pa_yy * fz * fx * pb_zzz

                     - fz * fga * pa_z * pb_zzzz

                     + 6.0 * fx * fz * pa_z * pb_zzzz

                     + 14.0 * pa_yyz * fz * pb_zzzz);

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
    inline double fvec_yzz_xxxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yzz,
                                    double pb_xx,
                                    double pb_xxxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_y * fx * fx * fz * fgb

                     - 0.75 * pa_y * fz * fga * fx * fx

                     - 3.0 * pa_yzz * fx * fz * fgb

                     + 3.0 * pa_y * fz * fx * fx * fx

                     + 7.5 * pa_yzz * fz * fx * fx

                     - 3.0 * pa_y * fx * pb_xx * fz * fgb

                     - 3.0 * pa_y * fz * fga * pb_xx * fx

                     - 6.0 * pa_yzz * pb_xx * fz * fgb

                     + 15.0 * pa_y * fz * fx * fx * pb_xx

                     + 36.0 * pa_yzz * fz * pb_xx * fx

                     - pa_y * fz * fga * pb_xxxx

                     + 6.0 * pa_y * fz * fx * pb_xxxx

                     + 14.0 * pa_yzz * fz * pb_xxxx);

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
    inline double fvec_yzz_xxxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yzz,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xxx,
                                    double pb_xxxy,
                                    double pb_xy,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pb_x * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pb_x

                     - 1.5 * fx * pa_zz * pb_x * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_x

                     + 7.5 * fx * fx * pa_zz * fz * pb_x

                     - 1.5 * pa_y * fx * pb_xy * fz * fgb

                     - 1.5 * pa_y * fz * fga * pb_xy * fx

                     - 0.5 * fx * fz * fga * pb_xxx

                     - 3.0 * pa_yzz * pb_xy * fz * fgb

                     + 7.5 * pa_y * fz * fx * fx * pb_xy

                     + 2.5 * fx * fx * fz * pb_xxx

                     + 18.0 * pa_yzz * fz * pb_xy * fx

                     + 6.0 * fx * pa_zz * fz * pb_xxx

                     - pa_y * fz * fga * pb_xxxy

                     + 6.0 * pa_y * fz * fx * pb_xxxy

                     + 14.0 * pa_yzz * fz * pb_xxxy);

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
    inline double fvec_yzz_xxxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pb_x,
                                    double pb_xxx,
                                    double pb_xxxz,
                                    double pb_xz,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_yz * fx * pb_x * fz * fgb

                     + 15.0 * pa_yz * fz * fx * fx * pb_x

                     - 1.5 * pa_y * fx * pb_xz * fz * fgb

                     - 1.5 * pa_y * fz * fga * pb_xz * fx

                     - 3.0 * pa_yzz * pb_xz * fz * fgb

                     + 7.5 * pa_y * fz * fx * fx * pb_xz

                     + 18.0 * pa_yzz * fz * pb_xz * fx

                     + 12.0 * pa_yz * fz * fx * pb_xxx

                     - pa_y * fz * fga * pb_xxxz

                     + 6.0 * pa_y * fz * fx * pb_xxxz

                     + 14.0 * pa_yzz * fz * pb_xxxz);

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
    inline double fvec_yzz_xxyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yzz,
                                    double pa_zz,
                                    double pb_xx,
                                    double pb_xxy,
                                    double pb_xxyy,
                                    double pb_y,
                                    double pb_yy,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_y * fx * fx * fz * fgb

                     - 0.25 * pa_y * fz * fga * fx * fx

                     - 0.5 * fx * fx * fz * fgb * pb_y

                     - 0.5 * fx * fx * fz * fga * pb_y

                     - pa_yzz * fx * fz * fgb

                     + pa_y * fz * fx * fx * fx

                     - fx * pa_zz * fz * fgb * pb_y

                     + 2.0 * fx * fx * fx * fz * pb_y

                     + 2.5 * pa_yzz * fz * fx * fx

                     + 5.0 * fx * fx * pa_zz * fz * pb_y

                     - 0.5 * pa_y * fx * pb_xx * fz * fgb

                     - 0.5 * pa_y * fx * fz * fgb * pb_yy

                     - 0.5 * pa_y * fz * fga * pb_xx * fx

                     - 0.5 * pa_y * fz * fga * fx * pb_yy

                     - fx * fz * fga * pb_xxy

                     - pa_yzz * pb_xx * fz * fgb

                     - pa_yzz * fz * fgb * pb_yy

                     + 2.5 * pa_y * fz * fx * fx * pb_xx

                     + 2.5 * pa_y * fz * fx * fx * pb_yy

                     + 5.0 * fx * fx * fz * pb_xxy

                     + 6.0 * pa_yzz * fz * pb_xx * fx

                     + 6.0 * pa_yzz * fz * fx * pb_yy

                     + 12.0 * fx * pa_zz * fz * pb_xxy

                     - pa_y * fz * fga * pb_xxyy

                     + 6.0 * pa_y * fz * fx * pb_xxyy

                     + 14.0 * pa_yzz * fz * pb_xxyy);

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
    inline double fvec_yzz_xxyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
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
                                    double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * pa_z * fz * fgb

                     + 2.0 * fx * fx * fx * pa_z * fz

                     - 0.25 * fx * fx * fz * fgb * pb_z

                     - 0.25 * fx * fx * fz * fga * pb_z

                     - pa_yz * fx * fz * fgb * pb_y

                     - 0.5 * fx * pa_zz * fz * fgb * pb_z

                     + fx * fx * fx * fz * pb_z

                     + 5.0 * pa_yz * fz * fx * fx * pb_y

                     + 2.5 * fx * fx * pa_zz * fz * pb_z

                     + 5.0 * fx * fx * pa_z * fz * pb_xx

                     - 0.5 * pa_y * fx * fz * fgb * pb_yz

                     - 0.5 * pa_y * fz * fga * fx * pb_yz

                     - 0.5 * fx * fz * fga * pb_xxz

                     - pa_yzz * fz * fgb * pb_yz

                     + 2.5 * pa_y * fz * fx * fx * pb_yz

                     + 2.5 * fx * fx * fz * pb_xxz

                     + 6.0 * pa_yzz * fz * fx * pb_yz

                     + 12.0 * pa_yz * fz * fx * pb_xxy

                     + 6.0 * fx * pa_zz * fz * pb_xxz

                     - pa_y * fz * fga * pb_xxyz

                     + 6.0 * pa_y * fz * fx * pb_xxyz

                     + 14.0 * pa_yzz * fz * pb_xxyz);

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
    inline double fvec_yzz_xxzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pb_xx,
                                    double pb_xxz,
                                    double pb_xxzz,
                                    double pb_z,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- pa_y * fx * fx * fz * fgb

                     + 3.0 * pa_y * fx * fx * fx * fz

                     - 0.25 * pa_y * fz * fga * fx * fx

                     - pa_yzz * fx * fz * fgb

                     - 2.0 * pa_yz * fx * fz * fgb * pb_z

                     + 2.5 * pa_yzz * fz * fx * fx

                     + 10.0 * pa_yz * fz * fx * fx * pb_z

                     + 7.5 * pa_y * fx * fx * fz * pb_xx

                     - 0.5 * pa_y * fx * pb_xx * fz * fgb

                     - 0.5 * pa_y * fx * fz * fgb * pb_zz

                     - 0.5 * pa_y * fz * fga * pb_xx * fx

                     - 0.5 * pa_y * fz * fga * fx * pb_zz

                     - pa_yzz * pb_xx * fz * fgb

                     - pa_yzz * fz * fgb * pb_zz

                     + 2.5 * pa_y * fz * fx * fx * pb_zz

                     + 6.0 * pa_yzz * fz * pb_xx * fx

                     + 6.0 * pa_yzz * fz * fx * pb_zz

                     + 24.0 * pa_yz * fz * fx * pb_xxz

                     - pa_y * fz * fga * pb_xxzz

                     + 6.0 * pa_y * fz * fx * pb_xxzz

                     + 14.0 * pa_yzz * fz * pb_xxzz);

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
    inline double fvec_yzz_xyyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yzz,
                                    double pa_zz,
                                    double pb_x,
                                    double pb_xy,
                                    double pb_xyy,
                                    double pb_xyyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pb_x * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pb_x

                     - 1.5 * fx * pa_zz * pb_x * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_x

                     + 7.5 * fx * fx * pa_zz * fz * pb_x

                     - 1.5 * pa_y * fx * pb_xy * fz * fgb

                     - 1.5 * pa_y * fz * fga * pb_xy * fx

                     - 1.5 * fx * fz * fga * pb_xyy

                     - 3.0 * pa_yzz * pb_xy * fz * fgb

                     + 7.5 * pa_y * fz * fx * fx * pb_xy

                     + 7.5 * fx * fx * fz * pb_xyy

                     + 18.0 * pa_yzz * fz * pb_xy * fx

                     + 18.0 * fx * pa_zz * fz * pb_xyy

                     - pa_y * fz * fga * pb_xyyy

                     + 6.0 * pa_y * fz * fx * pb_xyyy

                     + 14.0 * pa_yzz * fz * pb_xyyy);

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
    inline double fvec_yzz_xyyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
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
                                    double r_0_0)
    {
        return r_0_0 * (- pa_yz * fx * pb_x * fz * fgb

                     + 5.0 * pa_yz * fz * fx * fx * pb_x

                     + 10.0 * fx * fx * pa_z * fz * pb_xy

                     - 0.5 * pa_y * fx * pb_xz * fz * fgb

                     - 0.5 * pa_y * fz * fga * pb_xz * fx

                     - fx * fz * fga * pb_xyz

                     - pa_yzz * pb_xz * fz * fgb

                     + 2.5 * pa_y * fz * fx * fx * pb_xz

                     + 5.0 * fx * fx * fz * pb_xyz

                     + 6.0 * pa_yzz * fz * pb_xz * fx

                     + 12.0 * pa_yz * fz * fx * pb_xyy

                     + 12.0 * fx * pa_zz * fz * pb_xyz

                     - pa_y * fz * fga * pb_xyyz

                     + 6.0 * pa_y * fz * fx * pb_xyyz

                     + 14.0 * pa_yzz * fz * pb_xyyz);

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
    inline double fvec_yzz_xyzz_r_0(double fga,
                                    double fgb,
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
                                    double pb_xyzz,
                                    double pb_xz,
                                    double pb_xzz,
                                    double r_0_0)
    {
        return r_0_0 * (3.0 * fx * fx * fx * fz * pb_x

                     - 0.25 * fx * fx * pb_x * fz * fgb

                     - 0.25 * fx * fx * fz * fga * pb_x

                     - 0.5 * fx * pa_zz * pb_x * fz * fgb

                     + 7.5 * pa_y * fx * fx * fz * pb_xy

                     + 2.5 * fx * fx * pa_zz * fz * pb_x

                     + 10.0 * fx * fx * pa_z * fz * pb_xz

                     - 0.5 * pa_y * fx * pb_xy * fz * fgb

                     - 0.5 * pa_y * fz * fga * pb_xy * fx

                     - 0.5 * fx * fz * fga * pb_xzz

                     - pa_yzz * pb_xy * fz * fgb

                     + 2.5 * fx * fx * fz * pb_xzz

                     + 6.0 * pa_yzz * fz * pb_xy * fx

                     + 24.0 * pa_yz * fz * fx * pb_xyz

                     + 6.0 * fx * pa_zz * fz * pb_xzz

                     - pa_y * fz * fga * pb_xyzz

                     + 6.0 * pa_y * fz * fx * pb_xyzz

                     + 14.0 * pa_yzz * fz * pb_xyzz);

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
    inline double fvec_yzz_xzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pb_x,
                                    double pb_xz,
                                    double pb_xzz,
                                    double pb_xzzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_yz * fx * pb_x * fz * fgb

                     + 15.0 * pa_yz * fz * fx * fx * pb_x

                     + 22.5 * pa_y * fx * fx * fz * pb_xz

                     - 1.5 * pa_y * fx * pb_xz * fz * fgb

                     - 1.5 * pa_y * fz * fga * pb_xz * fx

                     - 3.0 * pa_yzz * pb_xz * fz * fgb

                     + 18.0 * pa_yzz * fz * pb_xz * fx

                     + 36.0 * pa_yz * fz * fx * pb_xzz

                     - pa_y * fz * fga * pb_xzzz

                     + 6.0 * pa_y * fz * fx * pb_xzzz

                     + 14.0 * pa_yzz * fz * pb_xzzz);

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
    inline double fvec_yzz_yyyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yzz,
                                    double pa_zz,
                                    double pb_y,
                                    double pb_yy,
                                    double pb_yyy,
                                    double pb_yyyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_y * fx * fx * fz * fgb

                     - 0.75 * pa_y * fz * fga * fx * fx

                     - 3.0 * fx * fx * pb_y * fz * fgb

                     - 3.0 * fx * fx * fz * fga * pb_y

                     - 3.0 * pa_yzz * fx * fz * fgb

                     + 3.0 * pa_y * fz * fx * fx * fx

                     - 6.0 * fx * pa_zz * pb_y * fz * fgb

                     + 12.0 * fx * fx * fx * fz * pb_y

                     + 7.5 * pa_yzz * fz * fx * fx

                     + 30.0 * fx * fx * pa_zz * fz * pb_y

                     - 3.0 * pa_y * fx * pb_yy * fz * fgb

                     - 3.0 * pa_y * fz * fga * pb_yy * fx

                     - 2.0 * fx * fz * fga * pb_yyy

                     - 6.0 * pa_yzz * pb_yy * fz * fgb

                     + 15.0 * pa_y * fz * fx * fx * pb_yy

                     + 10.0 * fx * fx * fz * pb_yyy

                     + 36.0 * pa_yzz * fz * pb_yy * fx

                     + 24.0 * fx * pa_zz * fz * pb_yyy

                     - pa_y * fz * fga * pb_yyyy

                     + 6.0 * pa_y * fz * fx * pb_yyyy

                     + 14.0 * pa_yzz * fz * pb_yyyy);

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
    inline double fvec_yzz_yyyz_r_0(double fga,
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
                                    double pb_yyy,
                                    double pb_yyyz,
                                    double pb_yyz,
                                    double pb_yz,
                                    double pb_z,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_z * fz * fgb

                     + 6.0 * fx * fx * fx * pa_z * fz

                     - 0.75 * fx * fx * fz * fgb * pb_z

                     - 0.75 * fx * fx * fz * fga * pb_z

                     - 3.0 * pa_yz * fx * pb_y * fz * fgb

                     - 1.5 * fx * pa_zz * fz * fgb * pb_z

                     + 3.0 * fx * fx * fx * fz * pb_z

                     + 15.0 * pa_yz * fz * fx * fx * pb_y

                     + 7.5 * fx * fx * pa_zz * fz * pb_z

                     + 15.0 * fx * fx * pa_z * fz * pb_yy

                     - 1.5 * pa_y * fx * pb_yz * fz * fgb

                     - 1.5 * pa_y * fz * fga * pb_yz * fx

                     - 1.5 * fx * fz * fga * pb_yyz

                     - 3.0 * pa_yzz * pb_yz * fz * fgb

                     + 7.5 * pa_y * fz * fx * fx * pb_yz

                     + 7.5 * fx * fx * fz * pb_yyz

                     + 18.0 * pa_yzz * fz * pb_yz * fx

                     + 12.0 * pa_yz * fz * fx * pb_yyy

                     + 18.0 * fx * pa_zz * fz * pb_yyz

                     - pa_y * fz * fga * pb_yyyz

                     + 6.0 * pa_y * fz * fx * pb_yyyz

                     + 14.0 * pa_yzz * fz * pb_yyyz);

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
    inline double fvec_yzz_yyzz_r_0(double fga,
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
                                    double pb_yyzz,
                                    double pb_yz,
                                    double pb_yzz,
                                    double pb_z,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- pa_y * fx * fx * fz * fgb

                     + 3.0 * pa_y * fx * fx * fx * fz

                     + 6.0 * fx * fx * fx * fz * pb_y

                     - 0.25 * pa_y * fz * fga * fx * fx

                     - 0.5 * fx * fx * pb_y * fz * fgb

                     - 0.5 * fx * fx * fz * fga * pb_y

                     - pa_yzz * fx * fz * fgb

                     - 2.0 * pa_yz * fx * fz * fgb * pb_z

                     - fx * pa_zz * pb_y * fz * fgb

                     + 2.5 * pa_yzz * fz * fx * fx

                     + 10.0 * pa_yz * fz * fx * fx * pb_z

                     + 7.5 * pa_y * fx * fx * fz * pb_yy

                     + 5.0 * fx * fx * pa_zz * fz * pb_y

                     + 20.0 * fx * fx * pa_z * fz * pb_yz

                     - 0.5 * pa_y * fx * pb_yy * fz * fgb

                     - 0.5 * pa_y * fx * fz * fgb * pb_zz

                     - 0.5 * pa_y * fz * fga * pb_yy * fx

                     - 0.5 * pa_y * fz * fga * fx * pb_zz

                     - fx * fz * fga * pb_yzz

                     - pa_yzz * pb_yy * fz * fgb

                     - pa_yzz * fz * fgb * pb_zz

                     + 2.5 * pa_y * fz * fx * fx * pb_zz

                     + 5.0 * fx * fx * fz * pb_yzz

                     + 6.0 * pa_yzz * fz * pb_yy * fx

                     + 6.0 * pa_yzz * fz * fx * pb_zz

                     + 24.0 * pa_yz * fz * fx * pb_yyz

                     + 12.0 * fx * pa_zz * fz * pb_yzz

                     - pa_y * fz * fga * pb_yyzz

                     + 6.0 * pa_y * fz * fx * pb_yyzz

                     + 14.0 * pa_yzz * fz * pb_yyzz);

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
    inline double fvec_yzz_yzzz_r_0(double fga,
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
                                    double pb_yzzz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * pa_z * fz * fgb

                     + 6.0 * fx * fx * fx * pa_z * fz

                     + 9.0 * fx * fx * fx * fz * pb_z

                     - 0.75 * fx * fx * pb_z * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pb_z

                     - 3.0 * pa_yz * fx * pb_y * fz * fgb

                     - 1.5 * fx * pa_zz * pb_z * fz * fgb

                     + 15.0 * pa_yz * fz * fx * fx * pb_y

                     + 22.5 * pa_y * fx * fx * fz * pb_yz

                     + 7.5 * fx * fx * pa_zz * fz * pb_z

                     + 15.0 * fx * fx * pa_z * fz * pb_zz

                     - 1.5 * pa_y * fx * pb_yz * fz * fgb

                     - 1.5 * pa_y * fz * fga * pb_yz * fx

                     - 0.5 * fx * fz * fga * pb_zzz

                     - 3.0 * pa_yzz * pb_yz * fz * fgb

                     + 2.5 * fx * fx * fz * pb_zzz

                     + 18.0 * pa_yzz * fz * pb_yz * fx

                     + 36.0 * pa_yz * fz * fx * pb_yzz

                     + 6.0 * fx * pa_zz * fz * pb_zzz

                     - pa_y * fz * fga * pb_yzzz

                     + 6.0 * pa_y * fz * fx * pb_yzzz

                     + 14.0 * pa_yzz * fz * pb_yzzz);

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
    inline double fvec_yzz_zzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_y,
                                    double pa_yz,
                                    double pa_yzz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double pb_zzzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_y * fx * fx * fz * fgb

                     + 15.0 * pa_y * fx * fx * fx * fz

                     - 0.75 * pa_y * fz * fga * fx * fx

                     - 3.0 * pa_yzz * fx * fz * fgb

                     - 12.0 * pa_yz * fx * pb_z * fz * fgb

                     + 7.5 * pa_yzz * fz * fx * fx

                     + 60.0 * pa_yz * fz * fx * fx * pb_z

                     + 45.0 * pa_y * fx * fx * fz * pb_zz

                     - 3.0 * pa_y * fx * pb_zz * fz * fgb

                     - 3.0 * pa_y * fz * fga * pb_zz * fx

                     - 6.0 * pa_yzz * pb_zz * fz * fgb

                     + 36.0 * pa_yzz * fz * pb_zz * fx

                     + 48.0 * pa_yz * fz * fx * pb_zzz

                     - pa_y * fz * fga * pb_zzzz

                     + 6.0 * pa_y * fz * fx * pb_zzzz

                     + 14.0 * pa_yzz * fz * pb_zzzz);

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
    inline double fvec_zzz_xxxx_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_xx,
                                    double pb_xxxx,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_z * fx * fx * fz * fgb

                     - 2.25 * pa_z * fz * fga * fx * fx

                     - 3.0 * pa_zzz * fx * fz * fgb

                     + 9.0 * pa_z * fz * fx * fx * fx

                     + 7.5 * pa_zzz * fz * fx * fx

                     - 9.0 * pa_z * fx * pb_xx * fz * fgb

                     - 9.0 * pa_z * fz * fga * pb_xx * fx

                     - 6.0 * pa_zzz * pb_xx * fz * fgb

                     + 45.0 * pa_z * fz * fx * fx * pb_xx

                     + 36.0 * pa_zzz * fz * pb_xx * fx

                     - 3.0 * pa_z * fz * fga * pb_xxxx

                     + 18.0 * pa_z * fz * fx * pb_xxxx

                     + 14.0 * pa_zzz * fz * pb_xxxx);

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
    inline double fvec_zzz_xxxy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_xxxy,
                                    double pb_xy,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_z * fx * pb_xy * fz * fgb

                     - 4.5 * pa_z * fz * fga * pb_xy * fx

                     - 3.0 * pa_zzz * pb_xy * fz * fgb

                     + 22.5 * pa_z * fz * fx * fx * pb_xy

                     + 18.0 * pa_zzz * fz * pb_xy * fx

                     - 3.0 * pa_z * fz * fga * pb_xxxy

                     + 18.0 * pa_z * fz * fx * pb_xxxy

                     + 14.0 * pa_zzz * fz * pb_xxxy);

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
    inline double fvec_zzz_xxxz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_x,
                                    double pb_xxx,
                                    double pb_xxxz,
                                    double pb_xz,
                                    double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pb_x * fz * fgb

                     - 2.25 * fx * fx * fz * fga * pb_x

                     - 4.5 * pa_zz * fx * pb_x * fz * fgb

                     + 9.0 * fx * fx * fx * fz * pb_x

                     + 22.5 * pa_zz * fz * fx * fx * pb_x

                     - 4.5 * pa_z * fx * pb_xz * fz * fgb

                     - 4.5 * pa_z * fz * fga * pb_xz * fx

                     - 1.5 * fx * fz * fga * pb_xxx

                     - 3.0 * pa_zzz * pb_xz * fz * fgb

                     + 22.5 * pa_z * fz * fx * fx * pb_xz

                     + 7.5 * fx * fx * fz * pb_xxx

                     + 18.0 * pa_zzz * fz * pb_xz * fx

                     + 18.0 * pa_zz * fz * fx * pb_xxx

                     - 3.0 * pa_z * fz * fga * pb_xxxz

                     + 18.0 * pa_z * fz * fx * pb_xxxz

                     + 14.0 * pa_zzz * fz * pb_xxxz);

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
    inline double fvec_zzz_xxyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_xx,
                                    double pb_xxyy,
                                    double pb_yy,
                                    double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_z * fx * fx * fz * fgb

                     - 0.75 * pa_z * fz * fga * fx * fx

                     - pa_zzz * fx * fz * fgb

                     + 3.0 * pa_z * fz * fx * fx * fx

                     + 2.5 * pa_zzz * fz * fx * fx

                     - 1.5 * pa_z * fx * pb_xx * fz * fgb

                     - 1.5 * pa_z * fx * fz * fgb * pb_yy

                     - 1.5 * pa_z * fz * fga * pb_xx * fx

                     - 1.5 * pa_z * fz * fga * fx * pb_yy

                     - pa_zzz * pb_xx * fz * fgb

                     - pa_zzz * fz * fgb * pb_yy

                     + 7.5 * pa_z * fz * fx * fx * pb_xx

                     + 7.5 * pa_z * fz * fx * fx * pb_yy

                     + 6.0 * pa_zzz * fz * pb_xx * fx

                     + 6.0 * pa_zzz * fz * fx * pb_yy

                     - 3.0 * pa_z * fz * fga * pb_xxyy

                     + 18.0 * pa_z * fz * fx * pb_xxyy

                     + 14.0 * pa_zzz * fz * pb_xxyy);

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
    inline double fvec_zzz_xxyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_xxy,
                                    double pb_xxyz,
                                    double pb_y,
                                    double pb_yz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb * pb_y

                     - 0.75 * fx * fx * fz * fga * pb_y

                     - 1.5 * pa_zz * fx * fz * fgb * pb_y

                     + 3.0 * fx * fx * fx * fz * pb_y

                     + 7.5 * pa_zz * fz * fx * fx * pb_y

                     - 1.5 * pa_z * fx * fz * fgb * pb_yz

                     - 1.5 * pa_z * fz * fga * fx * pb_yz

                     - 1.5 * fx * fz * fga * pb_xxy

                     - pa_zzz * fz * fgb * pb_yz

                     + 7.5 * pa_z * fz * fx * fx * pb_yz

                     + 7.5 * fx * fx * fz * pb_xxy

                     + 6.0 * pa_zzz * fz * fx * pb_yz

                     + 18.0 * pa_zz * fz * fx * pb_xxy

                     - 3.0 * pa_z * fz * fga * pb_xxyz

                     + 18.0 * pa_z * fz * fx * pb_xxyz

                     + 14.0 * pa_zzz * fz * pb_xxyz);

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
    inline double fvec_zzz_xxzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_xx,
                                    double pb_xxz,
                                    double pb_xxzz,
                                    double pb_z,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_z * fx * fx * fz * fgb

                     + 9.0 * pa_z * fx * fx * fx * fz

                     - 0.75 * pa_z * fz * fga * fx * fx

                     - 1.5 * fx * fx * fz * fgb * pb_z

                     - 1.5 * fx * fx * fz * fga * pb_z

                     - pa_zzz * fx * fz * fgb

                     - 3.0 * pa_zz * fx * fz * fgb * pb_z

                     + 6.0 * fx * fx * fx * fz * pb_z

                     + 2.5 * pa_zzz * fz * fx * fx

                     + 15.0 * pa_zz * fz * fx * fx * pb_z

                     + 22.5 * pa_z * fx * fx * fz * pb_xx

                     - 1.5 * pa_z * fx * pb_xx * fz * fgb

                     - 1.5 * pa_z * fx * fz * fgb * pb_zz

                     - 1.5 * pa_z * fz * fga * pb_xx * fx

                     - 1.5 * pa_z * fz * fga * fx * pb_zz

                     - 3.0 * fx * fz * fga * pb_xxz

                     - pa_zzz * pb_xx * fz * fgb

                     - pa_zzz * fz * fgb * pb_zz

                     + 7.5 * pa_z * fz * fx * fx * pb_zz

                     + 15.0 * fx * fx * fz * pb_xxz

                     + 6.0 * pa_zzz * fz * pb_xx * fx

                     + 6.0 * pa_zzz * fz * fx * pb_zz

                     + 36.0 * pa_zz * fz * fx * pb_xxz

                     - 3.0 * pa_z * fz * fga * pb_xxzz

                     + 18.0 * pa_z * fz * fx * pb_xxzz

                     + 14.0 * pa_zzz * fz * pb_xxzz);

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
    inline double fvec_zzz_xyyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_xy,
                                    double pb_xyyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_z * fx * pb_xy * fz * fgb

                     - 4.5 * pa_z * fz * fga * pb_xy * fx

                     - 3.0 * pa_zzz * pb_xy * fz * fgb

                     + 22.5 * pa_z * fz * fx * fx * pb_xy

                     + 18.0 * pa_zzz * fz * pb_xy * fx

                     - 3.0 * pa_z * fz * fga * pb_xyyy

                     + 18.0 * pa_z * fz * fx * pb_xyyy

                     + 14.0 * pa_zzz * fz * pb_xyyy);

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
    inline double fvec_zzz_xyyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_x,
                                    double pb_xyy,
                                    double pb_xyyz,
                                    double pb_xz,
                                    double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * pb_x * fz * fgb

                     - 0.75 * fx * fx * fz * fga * pb_x

                     - 1.5 * pa_zz * fx * pb_x * fz * fgb

                     + 3.0 * fx * fx * fx * fz * pb_x

                     + 7.5 * pa_zz * fz * fx * fx * pb_x

                     - 1.5 * pa_z * fx * pb_xz * fz * fgb

                     - 1.5 * pa_z * fz * fga * pb_xz * fx

                     - 1.5 * fx * fz * fga * pb_xyy

                     - pa_zzz * pb_xz * fz * fgb

                     + 7.5 * pa_z * fz * fx * fx * pb_xz

                     + 7.5 * fx * fx * fz * pb_xyy

                     + 6.0 * pa_zzz * fz * pb_xz * fx

                     + 18.0 * pa_zz * fz * fx * pb_xyy

                     - 3.0 * pa_z * fz * fga * pb_xyyz

                     + 18.0 * pa_z * fz * fx * pb_xyyz

                     + 14.0 * pa_zzz * fz * pb_xyyz);

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
    inline double fvec_zzz_xyzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_xy,
                                    double pb_xyz,
                                    double pb_xyzz,
                                    double r_0_0)
    {
        return r_0_0 * (22.5 * pa_z * fx * fx * fz * pb_xy

                     - 1.5 * pa_z * fx * pb_xy * fz * fgb

                     - 1.5 * pa_z * fz * fga * pb_xy * fx

                     - 3.0 * fx * fz * fga * pb_xyz

                     - pa_zzz * pb_xy * fz * fgb

                     + 15.0 * fx * fx * fz * pb_xyz

                     + 6.0 * pa_zzz * fz * pb_xy * fx

                     + 36.0 * pa_zz * fz * fx * pb_xyz

                     - 3.0 * pa_z * fz * fga * pb_xyzz

                     + 18.0 * pa_z * fz * fx * pb_xyzz

                     + 14.0 * pa_zzz * fz * pb_xyzz);

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
    inline double fvec_zzz_xzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_x,
                                    double pb_xz,
                                    double pb_xzz,
                                    double pb_xzzz,
                                    double r_0_0)
    {
        return r_0_0 * (15.0 * fx * fx * fx * fz * pb_x

                     - 2.25 * fx * fx * pb_x * fz * fgb

                     - 2.25 * fx * fx * fz * fga * pb_x

                     - 4.5 * pa_zz * fx * pb_x * fz * fgb

                     + 22.5 * pa_zz * fz * fx * fx * pb_x

                     + 67.5 * pa_z * fx * fx * fz * pb_xz

                     - 4.5 * pa_z * fx * pb_xz * fz * fgb

                     - 4.5 * pa_z * fz * fga * pb_xz * fx

                     - 4.5 * fx * fz * fga * pb_xzz

                     - 3.0 * pa_zzz * pb_xz * fz * fgb

                     + 22.5 * fx * fx * fz * pb_xzz

                     + 18.0 * pa_zzz * fz * pb_xz * fx

                     + 54.0 * pa_zz * fz * fx * pb_xzz

                     - 3.0 * pa_z * fz * fga * pb_xzzz

                     + 18.0 * pa_z * fz * fx * pb_xzzz

                     + 14.0 * pa_zzz * fz * pb_xzzz);

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
    inline double fvec_zzz_yyyy_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zzz,
                                    double pb_yy,
                                    double pb_yyyy,
                                    double r_0_0)
    {
        return r_0_0 * (- 4.5 * pa_z * fx * fx * fz * fgb

                     - 2.25 * pa_z * fz * fga * fx * fx

                     - 3.0 * pa_zzz * fx * fz * fgb

                     + 9.0 * pa_z * fz * fx * fx * fx

                     + 7.5 * pa_zzz * fz * fx * fx

                     - 9.0 * pa_z * fx * pb_yy * fz * fgb

                     - 9.0 * pa_z * fz * fga * pb_yy * fx

                     - 6.0 * pa_zzz * pb_yy * fz * fgb

                     + 45.0 * pa_z * fz * fx * fx * pb_yy

                     + 36.0 * pa_zzz * fz * pb_yy * fx

                     - 3.0 * pa_z * fz * fga * pb_yyyy

                     + 18.0 * pa_z * fz * fx * pb_yyyy

                     + 14.0 * pa_zzz * fz * pb_yyyy);

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
    inline double fvec_zzz_yyyz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_y,
                                    double pb_yyy,
                                    double pb_yyyz,
                                    double pb_yz,
                                    double r_0_0)
    {
        return r_0_0 * (- 2.25 * fx * fx * pb_y * fz * fgb

                     - 2.25 * fx * fx * fz * fga * pb_y

                     - 4.5 * pa_zz * fx * pb_y * fz * fgb

                     + 9.0 * fx * fx * fx * fz * pb_y

                     + 22.5 * pa_zz * fz * fx * fx * pb_y

                     - 4.5 * pa_z * fx * pb_yz * fz * fgb

                     - 4.5 * pa_z * fz * fga * pb_yz * fx

                     - 1.5 * fx * fz * fga * pb_yyy

                     - 3.0 * pa_zzz * pb_yz * fz * fgb

                     + 22.5 * pa_z * fz * fx * fx * pb_yz

                     + 7.5 * fx * fx * fz * pb_yyy

                     + 18.0 * pa_zzz * fz * pb_yz * fx

                     + 18.0 * pa_zz * fz * fx * pb_yyy

                     - 3.0 * pa_z * fz * fga * pb_yyyz

                     + 18.0 * pa_z * fz * fx * pb_yyyz

                     + 14.0 * pa_zzz * fz * pb_yyyz);

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
    inline double fvec_zzz_yyzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_yy,
                                    double pb_yyz,
                                    double pb_yyzz,
                                    double pb_z,
                                    double pb_zz,
                                    double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_z * fx * fx * fz * fgb

                     + 9.0 * pa_z * fx * fx * fx * fz

                     - 0.75 * pa_z * fz * fga * fx * fx

                     - 1.5 * fx * fx * fz * fgb * pb_z

                     - 1.5 * fx * fx * fz * fga * pb_z

                     - pa_zzz * fx * fz * fgb

                     - 3.0 * pa_zz * fx * fz * fgb * pb_z

                     + 6.0 * fx * fx * fx * fz * pb_z

                     + 2.5 * pa_zzz * fz * fx * fx

                     + 15.0 * pa_zz * fz * fx * fx * pb_z

                     + 22.5 * pa_z * fx * fx * fz * pb_yy

                     - 1.5 * pa_z * fx * pb_yy * fz * fgb

                     - 1.5 * pa_z * fx * fz * fgb * pb_zz

                     - 1.5 * pa_z * fz * fga * pb_yy * fx

                     - 1.5 * pa_z * fz * fga * fx * pb_zz

                     - 3.0 * fx * fz * fga * pb_yyz

                     - pa_zzz * pb_yy * fz * fgb

                     - pa_zzz * fz * fgb * pb_zz

                     + 7.5 * pa_z * fz * fx * fx * pb_zz

                     + 15.0 * fx * fx * fz * pb_yyz

                     + 6.0 * pa_zzz * fz * pb_yy * fx

                     + 6.0 * pa_zzz * fz * fx * pb_zz

                     + 36.0 * pa_zz * fz * fx * pb_yyz

                     - 3.0 * pa_z * fz * fga * pb_yyzz

                     + 18.0 * pa_z * fz * fx * pb_yyzz

                     + 14.0 * pa_zzz * fz * pb_yyzz);

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
    inline double fvec_zzz_yzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_y,
                                    double pb_yz,
                                    double pb_yzz,
                                    double pb_yzzz,
                                    double r_0_0)
    {
        return r_0_0 * (15.0 * fx * fx * fx * fz * pb_y

                     - 2.25 * fx * fx * pb_y * fz * fgb

                     - 2.25 * fx * fx * fz * fga * pb_y

                     - 4.5 * pa_zz * fx * pb_y * fz * fgb

                     + 22.5 * pa_zz * fz * fx * fx * pb_y

                     + 67.5 * pa_z * fx * fx * fz * pb_yz

                     - 4.5 * pa_z * fx * pb_yz * fz * fgb

                     - 4.5 * pa_z * fz * fga * pb_yz * fx

                     - 4.5 * fx * fz * fga * pb_yzz

                     - 3.0 * pa_zzz * pb_yz * fz * fgb

                     + 22.5 * fx * fx * fz * pb_yzz

                     + 18.0 * pa_zzz * fz * pb_yz * fx

                     + 54.0 * pa_zz * fz * fx * pb_yzz

                     - 3.0 * pa_z * fz * fga * pb_yzzz

                     + 18.0 * pa_z * fz * fx * pb_yzzz

                     + 14.0 * pa_zzz * fz * pb_yzzz);

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

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_zzzz_r_0(double fga,
                                    double fgb,
                                    double fx,
                                    double fz,
                                    double pa_z,
                                    double pa_zz,
                                    double pa_zzz,
                                    double pb_z,
                                    double pb_zz,
                                    double pb_zzz,
                                    double pb_zzzz,
                                    double r_0_0)
    {
        return r_0_0 * (- 13.5 * pa_z * fx * fx * fz * fgb

                     + 45.0 * pa_z * fx * fx * fx * fz

                     + 60.0 * fx * fx * fx * fz * pb_z

                     - 2.25 * pa_z * fz * fga * fx * fx

                     - 9.0 * fx * fx * pb_z * fz * fgb

                     - 9.0 * fx * fx * fz * fga * pb_z

                     - 3.0 * pa_zzz * fx * fz * fgb

                     - 18.0 * pa_zz * fx * pb_z * fz * fgb

                     + 7.5 * pa_zzz * fz * fx * fx

                     + 90.0 * pa_zz * fz * fx * fx * pb_z

                     + 135.0 * pa_z * fx * fx * fz * pb_zz

                     - 9.0 * pa_z * fx * pb_zz * fz * fgb

                     - 9.0 * pa_z * fz * fga * pb_zz * fx

                     - 6.0 * fx * fz * fga * pb_zzz

                     - 6.0 * pa_zzz * pb_zz * fz * fgb

                     + 30.0 * fx * fx * fz * pb_zzz

                     + 36.0 * pa_zzz * fz * pb_zz * fx

                     + 72.0 * pa_zz * fz * fx * pb_zzz

                     - 3.0 * pa_z * fz * fga * pb_zzzz

                     + 18.0 * pa_z * fz * fx * pb_zzzz

                     + 14.0 * pa_zzz * fz * pb_zzzz);

    }


} // kinrecfunc namespace

