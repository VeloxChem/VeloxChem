//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

namespace ovlvecfunc { // ovlvecfunc namespace

    // SIMD elementary functions for (S||P) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_0_x_s_0(double pb_x,
                               double s_0_0)
    {
        return s_0_0 * (pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_y_s_0(double pb_y,
                               double s_0_0)
    {
        return s_0_0 * (pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_z_s_0(double pb_z,
                               double s_0_0)
    {
        return s_0_0 * (pb_z);

    }

    // SIMD elementary functions for (P||S) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_x_0_s_0(double pa_x,
                               double s_0_0)
    {
        return s_0_0 * (pa_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_0_s_0(double pa_y,
                               double s_0_0)
    {
        return s_0_0 * (pa_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_0_s_0(double pa_z,
                               double s_0_0)
    {
        return s_0_0 * (pa_z);

    }

    // SIMD elementary functions for (S||D) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_0_xx_s_0(double fx,
                                double pb_xx,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx

                     + pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xy_s_0(double pb_xy,
                                double s_0_0)
    {
        return s_0_0 * (pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xz_s_0(double pb_xz,
                                double s_0_0)
    {
        return s_0_0 * (pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_yy_s_0(double fx,
                                double pb_yy,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx

                     + pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_yz_s_0(double pb_yz,
                                double s_0_0)
    {
        return s_0_0 * (pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_zz_s_0(double fx,
                                double pb_zz,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx

                     + pb_zz);

    }

    // SIMD elementary functions for (D||S) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_xx_0_s_0(double fx,
                                double pa_xx,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx

                     + pa_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_0_s_0(double pa_xy,
                                double s_0_0)
    {
        return s_0_0 * (pa_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_0_s_0(double pa_xz,
                                double s_0_0)
    {
        return s_0_0 * (pa_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_0_s_0(double fx,
                                double pa_yy,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx

                     + pa_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_0_s_0(double pa_yz,
                                double s_0_0)
    {
        return s_0_0 * (pa_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_0_s_0(double fx,
                                double pa_zz,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx

                     + pa_zz);

    }

    // SIMD elementary functions for (S||F) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_0_xxx_s_0(double fx,
                                 double pb_x,
                                 double pb_xxx,
                                 double s_0_0)
    {
        return s_0_0 * (1.5 * pb_x * fx

                     + pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xxy_s_0(double fx,
                                 double pb_xxy,
                                 double pb_y,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_y

                     + pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xxz_s_0(double fx,
                                 double pb_xxz,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_z

                     + pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xyy_s_0(double fx,
                                 double pb_x,
                                 double pb_xyy,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pb_x * fx

                     + pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xyz_s_0(double pb_xyz,
                                 double s_0_0)
    {
        return s_0_0 * (pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xzz_s_0(double fx,
                                 double pb_x,
                                 double pb_xzz,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pb_x * fx

                     + pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_yyy_s_0(double fx,
                                 double pb_y,
                                 double pb_yyy,
                                 double s_0_0)
    {
        return s_0_0 * (1.5 * pb_y * fx

                     + pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_yyz_s_0(double fx,
                                 double pb_yyz,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_z

                     + pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_yzz_s_0(double fx,
                                 double pb_y,
                                 double pb_yzz,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pb_y * fx

                     + pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_zzz_s_0(double fx,
                                 double pb_z,
                                 double pb_zzz,
                                 double s_0_0)
    {
        return s_0_0 * (1.5 * pb_z * fx

                     + pb_zzz);

    }

    // SIMD elementary functions for (F||S) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_0_s_0(double fx,
                                 double pa_x,
                                 double pa_xxx,
                                 double s_0_0)
    {
        return s_0_0 * (1.5 * pa_x * fx

                     + pa_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_0_s_0(double fx,
                                 double pa_xxy,
                                 double pa_y,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pa_y

                     + pa_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_0_s_0(double fx,
                                 double pa_xxz,
                                 double pa_z,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pa_z

                     + pa_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_0_s_0(double fx,
                                 double pa_x,
                                 double pa_xyy,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx

                     + pa_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_0_s_0(double pa_xyz,
                                 double s_0_0)
    {
        return s_0_0 * (pa_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_0_s_0(double fx,
                                 double pa_x,
                                 double pa_xzz,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx

                     + pa_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_0_s_0(double fx,
                                 double pa_y,
                                 double pa_yyy,
                                 double s_0_0)
    {
        return s_0_0 * (1.5 * pa_y * fx

                     + pa_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_0_s_0(double fx,
                                 double pa_yyz,
                                 double pa_z,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pa_z

                     + pa_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_0_s_0(double fx,
                                 double pa_y,
                                 double pa_yzz,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_y * fx

                     + pa_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_0_s_0(double fx,
                                 double pa_z,
                                 double pa_zzz,
                                 double s_0_0)
    {
        return s_0_0 * (1.5 * pa_z * fx

                     + pa_zzz);

    }

    // SIMD elementary functions for (S||G) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_0_xxxx_s_0(double fx,
                                  double pb_xx,
                                  double pb_xxxx,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx

                     + 3.0 * pb_xx * fx

                     + pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xxxy_s_0(double fx,
                                  double pb_xxxy,
                                  double pb_xy,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pb_xy * fx

                     + pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xxxz_s_0(double fx,
                                  double pb_xxxz,
                                  double pb_xz,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pb_xz * fx

                     + pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xxyy_s_0(double fx,
                                  double pb_xx,
                                  double pb_xxyy,
                                  double pb_yy,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pb_xx * fx

                     + 0.5 * fx * pb_yy

                     + pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xxyz_s_0(double fx,
                                  double pb_xxyz,
                                  double pb_yz,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_yz

                     + pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xxzz_s_0(double fx,
                                  double pb_xx,
                                  double pb_xxzz,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pb_xx * fx

                     + 0.5 * fx * pb_zz

                     + pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xyyy_s_0(double fx,
                                  double pb_xy,
                                  double pb_xyyy,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pb_xy * fx

                     + pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xyyz_s_0(double fx,
                                  double pb_xyyz,
                                  double pb_xz,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pb_xz * fx

                     + pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xyzz_s_0(double fx,
                                  double pb_xy,
                                  double pb_xyzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pb_xy * fx

                     + pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xzzz_s_0(double fx,
                                  double pb_xz,
                                  double pb_xzzz,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pb_xz * fx

                     + pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_yyyy_s_0(double fx,
                                  double pb_yy,
                                  double pb_yyyy,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx

                     + 3.0 * pb_yy * fx

                     + pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_yyyz_s_0(double fx,
                                  double pb_yyyz,
                                  double pb_yz,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pb_yz * fx

                     + pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_yyzz_s_0(double fx,
                                  double pb_yy,
                                  double pb_yyzz,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pb_yy * fx

                     + 0.5 * fx * pb_zz

                     + pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_yzzz_s_0(double fx,
                                  double pb_yz,
                                  double pb_yzzz,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pb_yz * fx

                     + pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_zzzz_s_0(double fx,
                                  double pb_zz,
                                  double pb_zzzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx

                     + 3.0 * pb_zz * fx

                     + pb_zzzz);

    }

    // SIMD elementary functions for (G||S) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_0_s_0(double fx,
                                  double pa_xx,
                                  double pa_xxxx,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx

                     + 3.0 * pa_xx * fx

                     + pa_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_0_s_0(double fx,
                                  double pa_xxxy,
                                  double pa_xy,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xy * fx

                     + pa_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_0_s_0(double fx,
                                  double pa_xxxz,
                                  double pa_xz,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xz * fx

                     + pa_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_0_s_0(double fx,
                                  double pa_xx,
                                  double pa_xxyy,
                                  double pa_yy,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_xx * fx

                     + 0.5 * fx * pa_yy

                     + pa_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_0_s_0(double fx,
                                  double pa_xxyz,
                                  double pa_yz,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pa_yz

                     + pa_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_0_s_0(double fx,
                                  double pa_xx,
                                  double pa_xxzz,
                                  double pa_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_xx * fx

                     + 0.5 * fx * pa_zz

                     + pa_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_0_s_0(double fx,
                                  double pa_xy,
                                  double pa_xyyy,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xy * fx

                     + pa_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_0_s_0(double fx,
                                  double pa_xyyz,
                                  double pa_xz,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xz * fx

                     + pa_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_0_s_0(double fx,
                                  double pa_xy,
                                  double pa_xyzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xy * fx

                     + pa_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_0_s_0(double fx,
                                  double pa_xz,
                                  double pa_xzzz,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xz * fx

                     + pa_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_0_s_0(double fx,
                                  double pa_yy,
                                  double pa_yyyy,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx

                     + 3.0 * pa_yy * fx

                     + pa_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_0_s_0(double fx,
                                  double pa_yyyz,
                                  double pa_yz,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_yz * fx

                     + pa_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_0_s_0(double fx,
                                  double pa_yy,
                                  double pa_yyzz,
                                  double pa_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_yy * fx

                     + 0.5 * fx * pa_zz

                     + pa_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_0_s_0(double fx,
                                  double pa_yz,
                                  double pa_yzzz,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_yz * fx

                     + pa_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_0_s_0(double fx,
                                  double pa_zz,
                                  double pa_zzzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx

                     + 3.0 * pa_zz * fx

                     + pa_zzzz);

    }


} // ovlrecfunc namespace

