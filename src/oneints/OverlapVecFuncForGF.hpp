//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

namespace ovlvecfunc { // ovlvecfunc namespace

    // SIMD elementary functions for (G||F) integrals

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


} // ovlrecfunc namespace

