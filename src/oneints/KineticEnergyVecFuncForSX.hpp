//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

namespace kinvecfunc { // kinvecfunc namespace

    // SIMD elementary functions for (S|T|P) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_0_x_s_0(double pb_x,
                               double s_0_0)
    {
        return s_0_0 * (pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_x_r_0(double fz,
                               double pb_x,
                               double r_0_0)
    {
        return r_0_0 * (2.0 * fz * pb_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_y_s_0(double pb_y,
                               double s_0_0)
    {
        return s_0_0 * (pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_y_r_0(double fz,
                               double pb_y,
                               double r_0_0)
    {
        return r_0_0 * (2.0 * fz * pb_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_z_s_0(double pb_z,
                               double s_0_0)
    {
        return s_0_0 * (pb_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_z_r_0(double fz,
                               double pb_z,
                               double r_0_0)
    {
        return r_0_0 * (2.0 * fz * pb_z);

    }

    // SIMD elementary functions for (P|T|S) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_x_0_s_0(double pa_x,
                               double s_0_0)
    {
        return s_0_0 * (pa_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_x_0_r_0(double fz,
                               double pa_x,
                               double r_0_0)
    {
        return r_0_0 * (2.0 * fz * pa_x);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_0_s_0(double pa_y,
                               double s_0_0)
    {
        return s_0_0 * (pa_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_y_0_r_0(double fz,
                               double pa_y,
                               double r_0_0)
    {
        return r_0_0 * (2.0 * fz * pa_y);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_0_s_0(double pa_z,
                               double s_0_0)
    {
        return s_0_0 * (pa_z);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_z_0_r_0(double fz,
                               double pa_z,
                               double r_0_0)
    {
        return r_0_0 * (2.0 * fz * pa_z);

    }

    // SIMD elementary functions for (S|T|D) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_0_xx_s_0(double fx,
                                double pb_xx,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx

                     + pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xx_r_0(double fgb,
                                double fx,
                                double fz,
                                double pb_xx,
                                double r_0_0)
    {
        return r_0_0 * (- fz * fgb

                     + fz * fx

                     + 4.0 * pb_xx * fz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xy_s_0(double pb_xy,
                                double s_0_0)
    {
        return s_0_0 * (pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xy_r_0(double fz,
                                double pb_xy,
                                double r_0_0)
    {
        return r_0_0 * (4.0 * pb_xy * fz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xz_s_0(double pb_xz,
                                double s_0_0)
    {
        return s_0_0 * (pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xz_r_0(double fz,
                                double pb_xz,
                                double r_0_0)
    {
        return r_0_0 * (4.0 * pb_xz * fz);

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
    inline double fvec_0_yy_r_0(double fgb,
                                double fx,
                                double fz,
                                double pb_yy,
                                double r_0_0)
    {
        return r_0_0 * (- fz * fgb

                     + fz * fx

                     + 4.0 * pb_yy * fz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_yz_s_0(double pb_yz,
                                double s_0_0)
    {
        return s_0_0 * (pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_yz_r_0(double fz,
                                double pb_yz,
                                double r_0_0)
    {
        return r_0_0 * (4.0 * pb_yz * fz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_zz_s_0(double fx,
                                double pb_zz,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx

                     + pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_zz_r_0(double fgb,
                                double fx,
                                double fz,
                                double pb_zz,
                                double r_0_0)
    {
        return r_0_0 * (- fz * fgb

                     + fz * fx

                     + 4.0 * pb_zz * fz);

    }

    // SIMD elementary functions for (D|T|S) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_xx_0_s_0(double fx,
                                double pa_xx,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx

                     + pa_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_0_r_0(double fga,
                                double fx,
                                double fz,
                                double pa_xx,
                                double r_0_0)
    {
        return r_0_0 * (- fz * fga

                     + fz * fx

                     + 4.0 * pa_xx * fz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_0_s_0(double pa_xy,
                                double s_0_0)
    {
        return s_0_0 * (pa_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_0_r_0(double fz,
                                double pa_xy,
                                double r_0_0)
    {
        return r_0_0 * (4.0 * pa_xy * fz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_0_s_0(double pa_xz,
                                double s_0_0)
    {
        return s_0_0 * (pa_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_0_r_0(double fz,
                                double pa_xz,
                                double r_0_0)
    {
        return r_0_0 * (4.0 * pa_xz * fz);

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
    inline double fvec_yy_0_r_0(double fga,
                                double fx,
                                double fz,
                                double pa_yy,
                                double r_0_0)
    {
        return r_0_0 * (- fz * fga

                     + fz * fx

                     + 4.0 * pa_yy * fz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_0_s_0(double pa_yz,
                                double s_0_0)
    {
        return s_0_0 * (pa_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_0_r_0(double fz,
                                double pa_yz,
                                double r_0_0)
    {
        return r_0_0 * (4.0 * pa_yz * fz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_0_s_0(double fx,
                                double pa_zz,
                                double s_0_0)
    {
        return s_0_0 * (0.5 * fx

                     + pa_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_0_r_0(double fga,
                                double fx,
                                double fz,
                                double pa_zz,
                                double r_0_0)
    {
        return r_0_0 * (- fz * fga

                     + fz * fx

                     + 4.0 * pa_zz * fz);

    }

    // SIMD elementary functions for (S|T|F) integrals

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
    inline double fvec_0_xxx_r_0(double fgb,
                                 double fx,
                                 double fz,
                                 double pb_x,
                                 double pb_xxx,
                                 double r_0_0)
    {
        return r_0_0 * (- 3.0 * pb_x * fz * fgb

                     + 6.0 * pb_x * fz * fx

                     + 6.0 * pb_xxx * fz);

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
    inline double fvec_0_xxy_r_0(double fgb,
                                 double fx,
                                 double fz,
                                 double pb_xxy,
                                 double pb_y,
                                 double r_0_0)
    {
        return r_0_0 * (- fz * fgb * pb_y

                     + 2.0 * fx * fz * pb_y

                     + 6.0 * pb_xxy * fz);

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
    inline double fvec_0_xxz_r_0(double fgb,
                                 double fx,
                                 double fz,
                                 double pb_xxz,
                                 double pb_z,
                                 double r_0_0)
    {
        return r_0_0 * (- fz * fgb * pb_z

                     + 2.0 * fx * fz * pb_z

                     + 6.0 * pb_xxz * fz);

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
    inline double fvec_0_xyy_r_0(double fgb,
                                 double fx,
                                 double fz,
                                 double pb_x,
                                 double pb_xyy,
                                 double r_0_0)
    {
        return r_0_0 * (- pb_x * fz * fgb

                     + 2.0 * pb_x * fz * fx

                     + 6.0 * pb_xyy * fz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xyz_s_0(double pb_xyz,
                                 double s_0_0)
    {
        return s_0_0 * (pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_0_xyz_r_0(double fz,
                                 double pb_xyz,
                                 double r_0_0)
    {
        return r_0_0 * (6.0 * pb_xyz * fz);

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
    inline double fvec_0_xzz_r_0(double fgb,
                                 double fx,
                                 double fz,
                                 double pb_x,
                                 double pb_xzz,
                                 double r_0_0)
    {
        return r_0_0 * (- pb_x * fz * fgb

                     + 2.0 * pb_x * fz * fx

                     + 6.0 * pb_xzz * fz);

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
    inline double fvec_0_yyy_r_0(double fgb,
                                 double fx,
                                 double fz,
                                 double pb_y,
                                 double pb_yyy,
                                 double r_0_0)
    {
        return r_0_0 * (- 3.0 * pb_y * fz * fgb

                     + 6.0 * pb_y * fz * fx

                     + 6.0 * pb_yyy * fz);

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
    inline double fvec_0_yyz_r_0(double fgb,
                                 double fx,
                                 double fz,
                                 double pb_yyz,
                                 double pb_z,
                                 double r_0_0)
    {
        return r_0_0 * (- fz * fgb * pb_z

                     + 2.0 * fx * fz * pb_z

                     + 6.0 * pb_yyz * fz);

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
    inline double fvec_0_yzz_r_0(double fgb,
                                 double fx,
                                 double fz,
                                 double pb_y,
                                 double pb_yzz,
                                 double r_0_0)
    {
        return r_0_0 * (- pb_y * fz * fgb

                     + 2.0 * pb_y * fz * fx

                     + 6.0 * pb_yzz * fz);

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

    #pragma omp declare simd notinbranch
    inline double fvec_0_zzz_r_0(double fgb,
                                 double fx,
                                 double fz,
                                 double pb_z,
                                 double pb_zzz,
                                 double r_0_0)
    {
        return r_0_0 * (- 3.0 * pb_z * fz * fgb

                     + 6.0 * pb_z * fz * fx

                     + 6.0 * pb_zzz * fz);

    }

    // SIMD elementary functions for (F|T|S) integrals

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
    inline double fvec_xxx_0_r_0(double fga,
                                 double fx,
                                 double fz,
                                 double pa_x,
                                 double pa_xxx,
                                 double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_x * fz * fga

                     + 6.0 * pa_x * fz * fx

                     + 6.0 * pa_xxx * fz);

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
    inline double fvec_xxy_0_r_0(double fga,
                                 double fx,
                                 double fz,
                                 double pa_xxy,
                                 double pa_y,
                                 double r_0_0)
    {
        return r_0_0 * (- fz * fga * pa_y

                     + 2.0 * fx * fz * pa_y

                     + 6.0 * pa_xxy * fz);

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
    inline double fvec_xxz_0_r_0(double fga,
                                 double fx,
                                 double fz,
                                 double pa_xxz,
                                 double pa_z,
                                 double r_0_0)
    {
        return r_0_0 * (- fz * fga * pa_z

                     + 2.0 * fx * fz * pa_z

                     + 6.0 * pa_xxz * fz);

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
    inline double fvec_xyy_0_r_0(double fga,
                                 double fx,
                                 double fz,
                                 double pa_x,
                                 double pa_xyy,
                                 double r_0_0)
    {
        return r_0_0 * (- pa_x * fz * fga

                     + 2.0 * pa_x * fz * fx

                     + 6.0 * pa_xyy * fz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_0_s_0(double pa_xyz,
                                 double s_0_0)
    {
        return s_0_0 * (pa_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_0_r_0(double fz,
                                 double pa_xyz,
                                 double r_0_0)
    {
        return r_0_0 * (6.0 * pa_xyz * fz);

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
    inline double fvec_xzz_0_r_0(double fga,
                                 double fx,
                                 double fz,
                                 double pa_x,
                                 double pa_xzz,
                                 double r_0_0)
    {
        return r_0_0 * (- pa_x * fz * fga

                     + 2.0 * pa_x * fz * fx

                     + 6.0 * pa_xzz * fz);

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
    inline double fvec_yyy_0_r_0(double fga,
                                 double fx,
                                 double fz,
                                 double pa_y,
                                 double pa_yyy,
                                 double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_y * fz * fga

                     + 6.0 * pa_y * fz * fx

                     + 6.0 * pa_yyy * fz);

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
    inline double fvec_yyz_0_r_0(double fga,
                                 double fx,
                                 double fz,
                                 double pa_yyz,
                                 double pa_z,
                                 double r_0_0)
    {
        return r_0_0 * (- fz * fga * pa_z

                     + 2.0 * fx * fz * pa_z

                     + 6.0 * pa_yyz * fz);

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
    inline double fvec_yzz_0_r_0(double fga,
                                 double fx,
                                 double fz,
                                 double pa_y,
                                 double pa_yzz,
                                 double r_0_0)
    {
        return r_0_0 * (- pa_y * fz * fga

                     + 2.0 * pa_y * fz * fx

                     + 6.0 * pa_yzz * fz);

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

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_0_r_0(double fga,
                                 double fx,
                                 double fz,
                                 double pa_z,
                                 double pa_zzz,
                                 double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_z * fz * fga

                     + 6.0 * pa_z * fz * fx

                     + 6.0 * pa_zzz * fz);

    }

    // SIMD elementary functions for (S|T|G) integrals

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
    inline double fvec_0_xxxx_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pb_xx,
                                  double pb_xxxx,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * fx * fz * fgb

                     - 6.0 * pb_xx * fz * fgb

                     + 3.0 * fx * fx * fz

                     + 18.0 * pb_xx * fz * fx

                     + 8.0 * pb_xxxx * fz);

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
    inline double fvec_0_xxxy_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pb_xxxy,
                                  double pb_xy,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pb_xy * fz * fgb

                     + 9.0 * pb_xy * fx * fz

                     + 8.0 * pb_xxxy * fz);

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
    inline double fvec_0_xxxz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pb_xxxz,
                                  double pb_xz,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pb_xz * fz * fgb

                     + 9.0 * pb_xz * fx * fz

                     + 8.0 * pb_xxxz * fz);

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
    inline double fvec_0_xxyy_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pb_xx,
                                  double pb_xxyy,
                                  double pb_yy,
                                  double r_0_0)
    {
        return r_0_0 * (- fx * fz * fgb

                     - pb_xx * fz * fgb

                     + fx * fx * fz

                     - fz * fgb * pb_yy

                     + 3.0 * pb_xx * fz * fx

                     + 3.0 * fx * pb_yy * fz

                     + 8.0 * pb_xxyy * fz);

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
    inline double fvec_0_xxyz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pb_xxyz,
                                  double pb_yz,
                                  double r_0_0)
    {
        return r_0_0 * (- fz * fgb * pb_yz

                     + 3.0 * fx * pb_yz * fz

                     + 8.0 * pb_xxyz * fz);

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
    inline double fvec_0_xxzz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pb_xx,
                                  double pb_xxzz,
                                  double pb_zz,
                                  double r_0_0)
    {
        return r_0_0 * (- fx * fz * fgb

                     - pb_xx * fz * fgb

                     + fx * fx * fz

                     - fz * fgb * pb_zz

                     + 3.0 * pb_xx * fz * fx

                     + 3.0 * fx * pb_zz * fz

                     + 8.0 * pb_xxzz * fz);

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
    inline double fvec_0_xyyy_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pb_xy,
                                  double pb_xyyy,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pb_xy * fz * fgb

                     + 9.0 * pb_xy * fz * fx

                     + 8.0 * pb_xyyy * fz);

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
    inline double fvec_0_xyyz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pb_xyyz,
                                  double pb_xz,
                                  double r_0_0)
    {
        return r_0_0 * (- pb_xz * fz * fgb

                     + 3.0 * pb_xz * fx * fz

                     + 8.0 * pb_xyyz * fz);

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
    inline double fvec_0_xyzz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pb_xy,
                                  double pb_xyzz,
                                  double r_0_0)
    {
        return r_0_0 * (- pb_xy * fz * fgb

                     + 3.0 * pb_xy * fz * fx

                     + 8.0 * pb_xyzz * fz);

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
    inline double fvec_0_xzzz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pb_xz,
                                  double pb_xzzz,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pb_xz * fz * fgb

                     + 9.0 * pb_xz * fz * fx

                     + 8.0 * pb_xzzz * fz);

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
    inline double fvec_0_yyyy_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pb_yy,
                                  double pb_yyyy,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * fx * fz * fgb

                     - 6.0 * pb_yy * fz * fgb

                     + 3.0 * fx * fx * fz

                     + 18.0 * pb_yy * fz * fx

                     + 8.0 * pb_yyyy * fz);

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
    inline double fvec_0_yyyz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pb_yyyz,
                                  double pb_yz,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pb_yz * fz * fgb

                     + 9.0 * pb_yz * fx * fz

                     + 8.0 * pb_yyyz * fz);

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
    inline double fvec_0_yyzz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pb_yy,
                                  double pb_yyzz,
                                  double pb_zz,
                                  double r_0_0)
    {
        return r_0_0 * (- fx * fz * fgb

                     - pb_yy * fz * fgb

                     + fx * fx * fz

                     - fz * fgb * pb_zz

                     + 3.0 * pb_yy * fz * fx

                     + 3.0 * fx * pb_zz * fz

                     + 8.0 * pb_yyzz * fz);

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
    inline double fvec_0_yzzz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pb_yz,
                                  double pb_yzzz,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pb_yz * fz * fgb

                     + 9.0 * pb_yz * fz * fx

                     + 8.0 * pb_yzzz * fz);

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

    #pragma omp declare simd notinbranch
    inline double fvec_0_zzzz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pb_zz,
                                  double pb_zzzz,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * fx * fz * fgb

                     - 6.0 * pb_zz * fz * fgb

                     + 3.0 * fx * fx * fz

                     + 18.0 * pb_zz * fz * fx

                     + 8.0 * pb_zzzz * fz);

    }

    // SIMD elementary functions for (G|T|S) integrals

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
    inline double fvec_xxxx_0_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_xx,
                                  double pa_xxxx,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * fx * fz * fga

                     - 6.0 * pa_xx * fz * fga

                     + 3.0 * fx * fx * fz

                     + 18.0 * pa_xx * fz * fx

                     + 8.0 * pa_xxxx * fz);

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
    inline double fvec_xxxy_0_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_xxxy,
                                  double pa_xy,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xy * fz * fga

                     + 9.0 * pa_xy * fx * fz

                     + 8.0 * pa_xxxy * fz);

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
    inline double fvec_xxxz_0_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_xxxz,
                                  double pa_xz,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xz * fz * fga

                     + 9.0 * pa_xz * fx * fz

                     + 8.0 * pa_xxxz * fz);

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
    inline double fvec_xxyy_0_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_xx,
                                  double pa_xxyy,
                                  double pa_yy,
                                  double r_0_0)
    {
        return r_0_0 * (- fx * fz * fga

                     - pa_xx * fz * fga

                     + fx * fx * fz

                     - fz * fga * pa_yy

                     + 3.0 * pa_xx * fz * fx

                     + 3.0 * fx * pa_yy * fz

                     + 8.0 * pa_xxyy * fz);

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
    inline double fvec_xxyz_0_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_xxyz,
                                  double pa_yz,
                                  double r_0_0)
    {
        return r_0_0 * (- fz * fga * pa_yz

                     + 3.0 * fx * pa_yz * fz

                     + 8.0 * pa_xxyz * fz);

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
    inline double fvec_xxzz_0_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_xx,
                                  double pa_xxzz,
                                  double pa_zz,
                                  double r_0_0)
    {
        return r_0_0 * (- fx * fz * fga

                     - pa_xx * fz * fga

                     + fx * fx * fz

                     - fz * fga * pa_zz

                     + 3.0 * pa_xx * fz * fx

                     + 3.0 * fx * pa_zz * fz

                     + 8.0 * pa_xxzz * fz);

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
    inline double fvec_xyyy_0_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_xy,
                                  double pa_xyyy,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xy * fz * fga

                     + 9.0 * pa_xy * fz * fx

                     + 8.0 * pa_xyyy * fz);

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
    inline double fvec_xyyz_0_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_xyyz,
                                  double pa_xz,
                                  double r_0_0)
    {
        return r_0_0 * (- pa_xz * fz * fga

                     + 3.0 * pa_xz * fx * fz

                     + 8.0 * pa_xyyz * fz);

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
    inline double fvec_xyzz_0_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_xy,
                                  double pa_xyzz,
                                  double r_0_0)
    {
        return r_0_0 * (- pa_xy * fz * fga

                     + 3.0 * pa_xy * fz * fx

                     + 8.0 * pa_xyzz * fz);

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
    inline double fvec_xzzz_0_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_xz,
                                  double pa_xzzz,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xz * fz * fga

                     + 9.0 * pa_xz * fz * fx

                     + 8.0 * pa_xzzz * fz);

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
    inline double fvec_yyyy_0_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_yy,
                                  double pa_yyyy,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * fx * fz * fga

                     - 6.0 * pa_yy * fz * fga

                     + 3.0 * fx * fx * fz

                     + 18.0 * pa_yy * fz * fx

                     + 8.0 * pa_yyyy * fz);

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
    inline double fvec_yyyz_0_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_yyyz,
                                  double pa_yz,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_yz * fz * fga

                     + 9.0 * pa_yz * fx * fz

                     + 8.0 * pa_yyyz * fz);

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
    inline double fvec_yyzz_0_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_yy,
                                  double pa_yyzz,
                                  double pa_zz,
                                  double r_0_0)
    {
        return r_0_0 * (- fx * fz * fga

                     - pa_yy * fz * fga

                     + fx * fx * fz

                     - fz * fga * pa_zz

                     + 3.0 * pa_yy * fz * fx

                     + 3.0 * fx * pa_zz * fz

                     + 8.0 * pa_yyzz * fz);

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
    inline double fvec_yzzz_0_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_yz,
                                  double pa_yzzz,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_yz * fz * fga

                     + 9.0 * pa_yz * fz * fx

                     + 8.0 * pa_yzzz * fz);

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

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_0_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_zz,
                                  double pa_zzzz,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * fx * fz * fga

                     - 6.0 * pa_zz * fz * fga

                     + 3.0 * fx * fx * fz

                     + 18.0 * pa_zz * fz * fx

                     + 8.0 * pa_zzzz * fz);

    }


} // kinrecfunc namespace

